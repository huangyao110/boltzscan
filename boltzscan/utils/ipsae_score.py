from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
from tqdm import tqdm


def calc_ipsae(cif_file, pae_file, ipsae_script=None):
    if ipsae_script is None:
        ipsae_script = Path(__file__).resolve().parent / 'ipsae.py'
    cmd = [sys.executable, str(ipsae_script), str(pae_file), str(cif_file), '10', '10']
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"ipsae script execution failed: {e}") from e


def read_ipsae(ipsae_file):
    df = pd.read_csv(ipsae_file, skipinitialspace=True, sep=' ')
    for col in df.columns:
        if df[col].dtype == 'object':
            df[col] = df[col].map(lambda x: x.strip() if isinstance(x, str) else x)
    return df


@dataclass
class IpsaeScoreSummary:
    output: Path
    predictions_dir: Path
    total_predictions: int
    scored_predictions: int
    matched_rows: int
    warnings: list[str]


@dataclass
class IpsaeValidationSummary:
    output: Path | None
    predictions_dir: Path
    total_predictions: int
    compared_predictions: int
    valid_predictions: int
    invalid_predictions: int
    warnings: list[str]
    result: pd.DataFrame


def score_ipsae_table(res_dir, score_file=None, output=None, id_col=None, force=False, processes=1):
    """Calculate IPSAE files and merge TF-DNA ipSAE/ipTM metrics into a score table."""
    predictions_dir = _resolve_predictions_dir(res_dir)
    scores, warnings, prediction_names = _collect_ipsae_scores(
        predictions_dir,
        force=force,
        processes=processes,
    )
    if not prediction_names:
        raise ValueError(f"No prediction structures found in {predictions_dir}")
    if score_file is None:
        score_df = pd.DataFrame({'model_id': sorted(prediction_names)})
        id_col = 'model_id'
    else:
        score_df = _read_score_table(score_file)
        id_col = _select_id_col(score_df, prediction_names, id_col)

    ids = score_df[id_col].astype(str)
    score_df['ipsae'] = ids.map(lambda name: scores.get(name, {}).get('ipsae'))
    score_df['ipsae_asym_min'] = ids.map(lambda name: scores.get(name, {}).get('ipsae_asym_min'))
    score_df['boltz_iptm'] = ids.map(lambda name: scores.get(name, {}).get('boltz_iptm'))
    score_df['ipsae_iptm'] = ids.map(lambda name: scores.get(name, {}).get('ipsae_iptm'))
    score_df['iptm_diff'] = (score_df['boltz_iptm'] - score_df['ipsae_iptm']).abs()
    score_df['boltz_iptm_global'] = ids.map(lambda name: scores.get(name, {}).get('boltz_iptm_global'))
    score_df['boltz_iptm_asym_min'] = ids.map(lambda name: scores.get(name, {}).get('boltz_iptm_asym_min'))
    score_df['ipsae_iptm_asym_min'] = ids.map(lambda name: scores.get(name, {}).get('ipsae_iptm_asym_min'))
    score_df['pDockQ'] = ids.map(lambda name: scores.get(name, {}).get('pDockQ'))
    # Model-input IDs are opaque hashes in the run-level workflow.  Preserve
    # the biological identifiers from the score table when available instead
    # of trying to parse a TF/target gene from the reusable model ID.
    if score_file is not None:
        if 'tf_name' in score_df.columns:
            score_df['TF'] = score_df['tf_name']
        elif 'TF' not in score_df.columns:
            score_df['TF'] = ids.map(_parse_tf)
        if 'sequence_name' in score_df.columns:
            score_df['TG'] = score_df['sequence_name']
        elif id_col == 'model_id':
            score_df['TG'] = None
        elif 'TG' not in score_df.columns:
            score_df['TG'] = ids.map(_parse_tg)
    score_df = score_df.sort_values('ipsae', ascending=False, na_position='last')
    score_df = _round_score_columns(score_df)

    matched_predictions = int(ids.isin(prediction_names).sum())
    if matched_predictions == 0:
        raise ValueError(
            f"No rows in {score_file} matched prediction directories in {predictions_dir}. "
            "The table must contain IDs matching the prediction directory names."
        )

    if output is None:
        output = predictions_dir / 'ipsae_scored.csv'
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    score_df.to_csv(output, index=False)

    return IpsaeScoreSummary(
        output=output,
        predictions_dir=predictions_dir,
        total_predictions=len(prediction_names),
        scored_predictions=len(scores),
        matched_rows=matched_predictions,
        warnings=warnings,
    )


def validate_ipsae_iptm(res_dir, output=None, tolerance=0.05, force=False, processes=1):
    """Compare IPSAE ipTM_af values against Boltz pair_chains_iptm values."""
    predictions_dir = _resolve_predictions_dir(res_dir)
    rows, warnings = _collect_validation_rows(
        predictions_dir,
        tolerance=tolerance,
        force=force,
        processes=processes,
    )
    result = pd.DataFrame(
        rows,
        columns=[
            'name',
            'boltz_iptm',
            'ipsae_iptm',
            'iptm_diff',
            'boltz_iptm_global',
            'boltz_iptm_asym_min',
            'ipsae_iptm_asym_min',
            'tolerance',
            'valid',
        ],
    )
    compared = int(result['iptm_diff'].notna().sum()) if not result.empty else 0
    valid = int(result['valid'].eq(True).sum()) if not result.empty else 0
    invalid = int(result['valid'].eq(False).sum()) if not result.empty else 0

    output_path = None
    if output is not None:
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_path, index=False)

    return IpsaeValidationSummary(
        output=output_path,
        predictions_dir=predictions_dir,
        total_predictions=len(rows),
        compared_predictions=compared,
        valid_predictions=valid,
        invalid_predictions=invalid,
        warnings=warnings,
        result=result,
    )


def _resolve_predictions_dir(res_dir):
    res_dir = Path(res_dir)
    if not res_dir.exists():
        raise FileNotFoundError(f"Result directory not found: {res_dir}")
    if any(res_dir.glob('*.cif')):
        return res_dir
    if res_dir.name == 'predictions':
        return res_dir
    # The concise run workflow keeps each arm directly under
    # ``predictions/full`` or ``predictions/crop<N>``.  Treat a directory of
    # per-model folders as a predictions root without requiring an extra,
    # confusing ``predictions/predictions`` layer.
    if any(
        child.is_dir() and any(child.glob('*.cif'))
        for child in res_dir.iterdir()
    ):
        return res_dir
    predictions_dir = res_dir / 'predictions'
    if predictions_dir.is_dir():
        return predictions_dir
    raise FileNotFoundError(
        f"Could not find a predictions directory from {res_dir}. "
        "Pass either a Boltz run directory or its predictions/ subdirectory."
    )


def _collect_ipsae_scores(predictions_dir, force=False, processes=1):
    prediction_dirs = sorted(p for p in predictions_dir.iterdir() if p.is_dir())
    flat_cifs = sorted(p for p in predictions_dir.glob('*.cif') if p.is_file())
    scores = {}
    warnings = []
    processes = max(1, int(processes))

    if flat_cifs:
        jobs = flat_cifs
        score_one = _score_flat_prediction
        prediction_names = {path.stem for path in flat_cifs}
    else:
        jobs = prediction_dirs
        score_one = _score_prediction_dir
        prediction_names = {path.name for path in prediction_dirs}

    if processes == 1:
        results = (score_one(job, force) for job in jobs)
        iterator = tqdm(results, total=len(jobs), desc='calc_ipsae')
        for pred_name, metrics, pred_warnings in iterator:
            warnings.extend(pred_warnings)
            if metrics is not None:
                scores[pred_name] = metrics
    else:
        with ThreadPoolExecutor(max_workers=processes) as executor:
            futures = {
                executor.submit(score_one, job, force): job
                for job in jobs
            }
            iterator = tqdm(as_completed(futures), total=len(futures), desc='calc_ipsae')
            for future in iterator:
                pred_name, metrics, pred_warnings = future.result()
                warnings.extend(pred_warnings)
                if metrics is not None:
                    scores[pred_name] = metrics

    return scores, warnings, prediction_names


def _collect_validation_rows(predictions_dir, tolerance, force=False, processes=1):
    prediction_dirs = sorted(p for p in predictions_dir.iterdir() if p.is_dir())
    processes = max(1, int(processes))
    rows = []
    warnings = []

    if processes == 1:
        results = (
            _validate_prediction_dir(pred_dir, tolerance=tolerance, force=force)
            for pred_dir in prediction_dirs
        )
        iterator = tqdm(results, total=len(prediction_dirs), desc='valid')
        for row, pred_warnings in iterator:
            rows.append(row)
            warnings.extend(pred_warnings)
    else:
        with ThreadPoolExecutor(max_workers=processes) as executor:
            futures = {
                executor.submit(
                    _validate_prediction_dir,
                    pred_dir,
                    tolerance=tolerance,
                    force=force,
                ): pred_dir
                for pred_dir in prediction_dirs
            }
            iterator = tqdm(as_completed(futures), total=len(futures), desc='valid')
            for future in iterator:
                row, pred_warnings = future.result()
                rows.append(row)
                warnings.extend(pred_warnings)

    rows.sort(key=lambda row: row['name'])
    return rows, warnings


def _score_prediction_dir(pred_dir, force):
    warnings = []
    cif_file, cif_warning = _select_prediction_file(pred_dir, '*.cif')
    pae_file, pae_warning = _select_prediction_file(pred_dir, 'pae*')
    json_file, json_warning = _select_prediction_file(pred_dir, 'confidence*.json')
    for warning in (cif_warning, pae_warning, json_warning):
        if warning:
            warnings.append(f"{pred_dir.name}: {warning}")
    return _score_prediction_files(
        pred_dir.name,
        cif_file,
        pae_file,
        json_file,
        force,
        warnings,
    )


def _score_flat_prediction(cif_file, force):
    name = cif_file.stem
    pred_dir = cif_file.parent
    warnings = []
    pae_file = pred_dir / f'pae_{name}.npz'
    if not pae_file.exists():
        pae_file = pred_dir / f'pae_{name}.json'
    if not pae_file.exists():
        warnings.append(f"{name}: missing pae_{name}.npz or pae_{name}.json")
        pae_file = None
    json_file = pred_dir / f'confidence_{name}.json'
    if not json_file.exists():
        json_file = None
    return _score_prediction_files(name, cif_file, pae_file, json_file, force, warnings)


def _score_prediction_files(name, cif_file, pae_file, json_file, force, warnings):
    ipsae_file, warnings = _ensure_ipsae_files(
        name,
        cif_file,
        pae_file,
        force=force,
        warnings=warnings,
    )
    if ipsae_file is None:
        return name, None, warnings

    try:
        metrics = _read_ipsae_metrics(ipsae_file)
        if json_file is not None:
            try:
                metrics.update(_read_boltz_iptm_metrics(json_file))
            except Exception as exc:
                warnings.append(f"{name}: could not read Boltz iptm: {exc}")
                metrics['boltz_iptm'] = None
                metrics['boltz_iptm_global'] = None
                metrics['boltz_iptm_asym_min'] = None
        else:
            metrics['boltz_iptm'] = None
            metrics['boltz_iptm_global'] = None
            metrics['boltz_iptm_asym_min'] = None
        return name, metrics, warnings
    except Exception as exc:
        warnings.append(f"{name}: could not read IPSAE metrics: {exc}")
        return name, None, warnings


def _validate_prediction_dir(pred_dir, tolerance, force):
    warnings = []
    row = {
        'name': pred_dir.name,
        'boltz_iptm': None,
        'boltz_iptm_global': None,
        'boltz_iptm_asym_min': None,
        'ipsae_iptm': None,
        'ipsae_iptm_asym_min': None,
        'iptm_diff': None,
        'tolerance': tolerance,
        'valid': None,
    }

    json_file, json_warning = _select_prediction_file(pred_dir, 'confidence*.json')
    if json_warning:
        warnings.append(f"{pred_dir.name}: {json_warning}")
    if json_file is not None:
        try:
            row.update(_read_boltz_iptm_metrics(json_file))
        except Exception as exc:
            warnings.append(f"{pred_dir.name}: could not read Boltz iptm: {exc}")

    ipsae_file, ipsae_warnings = _ensure_ipsae_file(pred_dir, force=force)
    warnings.extend(ipsae_warnings)
    if ipsae_file is not None:
        try:
            metrics = _read_ipsae_metrics(ipsae_file)
            row['ipsae_iptm'] = metrics['ipsae_iptm']
            row['ipsae_iptm_asym_min'] = metrics['ipsae_iptm_asym_min']
        except Exception as exc:
            warnings.append(f"{pred_dir.name}: could not read IPSAE ipTM_af: {exc}")

    if row['boltz_iptm'] is not None and row['ipsae_iptm'] is not None:
        row['iptm_diff'] = abs(row['boltz_iptm'] - row['ipsae_iptm'])
        row['valid'] = row['iptm_diff'] <= tolerance

    return row, warnings


def _ensure_ipsae_file(pred_dir, force=False):
    warnings = []
    cif_file, cif_warning = _select_prediction_file(pred_dir, '*.cif')
    pae_file, pae_warning = _select_prediction_file(pred_dir, 'pae*')
    if cif_warning:
        warnings.append(f"{pred_dir.name}: {cif_warning}")
    if pae_warning:
        warnings.append(f"{pred_dir.name}: {pae_warning}")
    return _ensure_ipsae_files(
        pred_dir.name,
        cif_file,
        pae_file,
        force=force,
        warnings=warnings,
    )


def _ensure_ipsae_files(name, cif_file, pae_file, force, warnings):
    if cif_file is None or pae_file is None:
        return None, warnings

    ipsae_file = _expected_ipsae_file(cif_file)
    if force or not ipsae_file.exists():
        try:
            calc_ipsae(cif_file=cif_file, pae_file=pae_file)
        except Exception as exc:
            warnings.append(f"{name}: IPSAE calculation failed: {exc}")
            return None, warnings

    if not ipsae_file.exists():
        ipsae_file = _find_existing_ipsae_file(cif_file.parent)
    if ipsae_file is None:
        warnings.append(f"{name}: IPSAE output txt not found")
        return None, warnings

    return ipsae_file, warnings


def _read_boltz_iptm_metrics(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    if 'pair_chains_iptm' not in data:
        raise ValueError("missing key: pair_chains_iptm")

    pair_iptm = data['pair_chains_iptm']
    ab = _boltz_pairmax_value(pair_iptm, 'A', 'B')
    ac = _boltz_pairmax_value(pair_iptm, 'A', 'C')
    asym_values = [
        _boltz_pair_value(pair_iptm, 'A', 'B'),
        _boltz_pair_value(pair_iptm, 'B', 'A'),
        _boltz_pair_value(pair_iptm, 'A', 'C'),
        _boltz_pair_value(pair_iptm, 'C', 'A'),
    ]

    metrics = {
        'boltz_iptm': min(ab, ac),
        'boltz_iptm_asym_min': min(asym_values),
        'boltz_iptm_global': None,
    }
    if 'iptm' in data and data['iptm'] is not None:
        metrics['boltz_iptm_global'] = float(data['iptm'])
    return metrics


def _boltz_pairmax_value(pair_iptm, chain1, chain2):
    return max(
        _boltz_pair_value(pair_iptm, chain1, chain2),
        _boltz_pair_value(pair_iptm, chain2, chain1),
    )


def _boltz_pair_value(pair_iptm, chain1, chain2):
    chain_indexes = {'A': '0', 'B': '1', 'C': '2'}
    row_idx = chain_indexes[chain1]
    col_idx = chain_indexes[chain2]
    try:
        row = pair_iptm[row_idx]
    except (KeyError, TypeError):
        row = pair_iptm[int(row_idx)]

    try:
        value = row[col_idx]
    except (KeyError, TypeError):
        value = row[int(col_idx)]
    return float(value)


def _select_prediction_file(pred_dir, pattern):
    candidates = sorted(pred_dir.glob(pattern))
    if pattern == 'pae*':
        candidates = [p for p in candidates if p.suffix in {'.npz', '.json'}]
    if not candidates:
        return None, f"missing {pattern}"
    if len(candidates) == 1:
        return candidates[0], None

    named_matches = [p for p in candidates if pred_dir.name in p.name]
    if len(named_matches) == 1:
        return named_matches[0], None
    return None, f"multiple {pattern} files found"


def _expected_ipsae_file(cif_file):
    return cif_file.with_name(f"{cif_file.stem}_10_10.txt")


def _find_existing_ipsae_file(pred_dir):
    txt_files = sorted(
        p for p in pred_dir.glob('*.txt')
        if 'byres' not in p.stem and p.stem.endswith('_10_10')
    )
    if txt_files:
        return txt_files[0]
    txt_files = sorted(p for p in pred_dir.glob('*.txt') if 'byres' not in p.stem)
    return txt_files[0] if txt_files else None


def _read_ipsae_metrics(ipsae_file):
    df = read_ipsae(ipsae_file)
    required = {'Type', 'Chn1', 'Chn2', 'ipSAE', 'ipTM_af', 'pDockQ'}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"missing columns: {', '.join(sorted(missing))}")

    asym_rows = _tf_dna_asym_rows(df)
    if asym_rows.empty:
        raise ValueError("no asymmetric TF-DNA rows for A-B/B-A/A-C/C-A")

    max_rows = _tf_dna_max_rows(df)
    ipsae = _tf_dna_pairmax_min(df, 'ipSAE')
    ipsae_iptm = _tf_dna_pairmax_min(df, 'ipTM_af')
    ipsae_asym_min = pd.to_numeric(asym_rows['ipSAE'], errors='coerce').min()
    ipsae_iptm_asym_min = pd.to_numeric(asym_rows['ipTM_af'], errors='coerce').min()
    pdockq_rows = max_rows
    pdockq = pd.to_numeric(pdockq_rows['pDockQ'], errors='coerce').max()
    return {
        'ipsae': None if pd.isna(ipsae) else float(ipsae),
        'ipsae_asym_min': None if pd.isna(ipsae_asym_min) else float(ipsae_asym_min),
        'ipsae_iptm': None if pd.isna(ipsae_iptm) else float(ipsae_iptm),
        'ipsae_iptm_asym_min': (
            None if pd.isna(ipsae_iptm_asym_min) else float(ipsae_iptm_asym_min)
        ),
        'pDockQ': None if pd.isna(pdockq) else float(pdockq),
    }


def _tf_dna_asym_rows(df):
    chain1 = df['Chn1'].astype(str)
    chain2 = df['Chn2'].astype(str)
    return df[
        (df['Type'].astype(str) == 'asym')
        & (
            ((chain1 == 'A') & chain2.isin(['B', 'C']))
            | ((chain2 == 'A') & chain1.isin(['B', 'C']))
        )
    ]


def _tf_dna_max_rows(df):
    chain1 = df['Chn1'].astype(str)
    chain2 = df['Chn2'].astype(str)
    return df[
        (df['Type'].astype(str) == 'max')
        & (chain1 == 'A')
        & chain2.isin(['B', 'C'])
    ]


def _tf_dna_pairmax_min(df, value_col):
    values = []
    max_rows = _tf_dna_max_rows(df)
    for dna_chain in ['B', 'C']:
        subset = max_rows[
            (max_rows['Chn1'].astype(str) == 'A')
            & (max_rows['Chn2'].astype(str) == dna_chain)
        ]
        if subset.empty:
            subset = _tf_dna_asym_rows(df)
            subset = subset[
                (
                    (subset['Chn1'].astype(str) == 'A')
                    & (subset['Chn2'].astype(str) == dna_chain)
                )
                | (
                    (subset['Chn1'].astype(str) == dna_chain)
                    & (subset['Chn2'].astype(str) == 'A')
                )
            ]
        if subset.empty:
            raise ValueError(f"missing TF-DNA rows for A-{dna_chain}")

        value = pd.to_numeric(subset[value_col], errors='coerce').max()
        if pd.isna(value):
            raise ValueError(f"{value_col} is not numeric for A-{dna_chain}")
        values.append(float(value))

    return min(values)


def _read_ipsae_iptm(ipsae_file):
    return _read_ipsae_metrics(ipsae_file)['ipsae_iptm']


def _read_score_table(score_file):
    score_file = Path(score_file)
    if not score_file.exists():
        raise FileNotFoundError(f"Score file not found: {score_file}")
    sep = '\t' if score_file.name.endswith(('.tsv', '.tsv.gz')) else ','
    return pd.read_csv(score_file, sep=sep)


def _round_score_columns(score_df):
    for col in [
        'ipsae',
        'ipsae_asym_min',
        'boltz_iptm',
        'ipsae_iptm',
        'iptm_diff',
        'boltz_iptm_global',
        'boltz_iptm_asym_min',
        'ipsae_iptm_asym_min',
        'pDockQ',
    ]:
        if col in score_df.columns:
            score_df[col] = pd.to_numeric(score_df[col], errors='coerce').round(6)
    return score_df


def _select_id_col(score_df, prediction_names, id_col=None):
    if score_df.empty:
        raise ValueError("Score file is empty")
    if id_col is not None:
        if id_col not in score_df.columns:
            raise ValueError(f"ID column not found: {id_col}")
        return id_col

    prediction_names = set(prediction_names)
    preferred = ['model_id', 'boltz_name', 'name', 'id']
    candidates = [c for c in preferred if c in score_df.columns]
    candidates.extend(score_df.columns)

    seen = set()
    best_col = None
    best_overlap = -1
    for col in candidates:
        if col in seen:
            continue
        seen.add(col)
        overlap = score_df[col].astype(str).isin(prediction_names).sum()
        if overlap > best_overlap:
            best_col = col
            best_overlap = overlap

    if best_col is None or best_overlap <= 0:
        raise ValueError(
            "Could not infer the prediction ID column; the table should contain model_id."
        )
    return best_col


def _parse_tf(name):
    parts = str(name).split('-')
    return parts[0] if parts else None


def _parse_tg(name):
    parts = str(name).split('-')
    return parts[-3] if len(parts) >= 3 else None


def print_ipsae_warnings(warnings, max_warnings=20):
    if not warnings:
        return
    print(f"Warnings ({len(warnings)}):", file=sys.stderr)
    for warning in warnings[:max_warnings]:
        print(f"- {warning}", file=sys.stderr)
    if len(warnings) > max_warnings:
        print(f"- ... {len(warnings) - max_warnings} more", file=sys.stderr)
