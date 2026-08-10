"""Build promoter/model candidate layers and TF--DNA prediction YAMLs."""
from dataclasses import dataclass
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
import yaml
from Bio.Seq import Seq
from Bio import SeqIO


def _read_fasta_to_str(path):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, 'fasta')}


def _build_msa_dct(msa_dir, a3m_name='0.a3m'):
    """Walk `msa_dir` and return {subdir_name: <subdir>/<a3m_name>} for every
    subdirectory that actually contains the a3m file."""
    msa_dir = Path(msa_dir)
    dct = {}
    if not msa_dir.exists():
        return dct
    for sub in msa_dir.iterdir():
        if sub.is_dir():
            a3m = sub / a3m_name
            if a3m.exists():
                dct[sub.name] = str(a3m)
    return dct


def _resolve_msa_path(tf_name, msa_paths):
    """Resolve a TF record id against Protenix MSA directory naming variants.

    Older local runs retained `-` and replaced only `.` with `_`. The
    current Protenix splitter sanitises both punctuation characters, e.g.
    `CDS-Solyc05g012020.4.1` -> `CDS_Solyc05g012020_4_1`. Try exact
    names first so an intentionally named directory always wins.
    """
    candidates = (
        tf_name,
        tf_name.replace('.', '_'),
        re.sub(r'[^A-Za-z0-9_]', '_', tf_name),
    )
    for candidate in candidates:
        path = msa_paths.get(candidate)
        if path:
            return path
    return None


def _invert_tf2pwms(tf2pwms):
    """{tf: [pwm, ...]} -> {pwm: [tf, ...]}"""
    inv = {}
    for tf, pwms in tf2pwms.items():
        for pwm in pwms:
            inv.setdefault(str(pwm), []).append(str(tf))
    return inv


@dataclass(frozen=True)
class CropInputSummary:
    """Paths and counts for a coordinate-cropped structural input arm."""

    crop_input_dir: Path
    crop_manifest: Path
    n_model_inputs: int


def canonical_duplex(sequence):
    """Return a strand-invariant representation of a dsDNA sequence.

    This is deliberately a **model-input** operation.  Callers must apply it
    only after promoter-level FIMO de-duplication, because identical DNA on
    different promoters is still a distinct biological candidate edge.
    """
    sequence = str(sequence).upper()
    if not sequence:
        raise ValueError("Cannot canonicalize an empty DNA sequence")
    reverse_complement = str(Seq(sequence).reverse_complement()).upper()
    return min(sequence, reverse_complement)


def _stable_id(prefix, *parts):
    payload = '\0'.join(str(part) for part in parts).encode('utf-8')
    return f"{prefix}_{hashlib.sha256(payload).hexdigest()[:20]}"


def _candidate_id(row):
    """Stable identity for one retained promoter-level FIMO candidate."""
    return _stable_id(
        'candidate',
        row['tf_name'],
        row['motif_id'],
        row['sequence_name'],
        row['start'],
        row['stop'],
        row.get('strand', ''),
        str(row['extracted_motif_seq']).upper(),
    )


def _candidate_sort_columns(df):
    """Return deterministic FIMO priority columns present in ``df``."""
    candidates = ['pvalue', 'score', 'sequence_name', 'start', 'stop', 'strand', 'motif_id']
    return [column for column in candidates if column in df.columns]


def _candidate_sort_ascending(columns):
    return [column != 'score' for column in columns]


def build_promoter_candidates_from_fimo(
    fimo_csv,
    tf_fasta,
    tf2pwms_path,
    promoter_fasta,
    pvalue_thresh=1e-5,
    overlap_thresh=0.0,
    dna_flank=5,
):
    """Build the biological, promoter-level TF--DNA candidate table.

    The de-duplication rule is intentionally local to ``(sequence_name,
    tf_name)``.  Thus, an identical dsDNA bait on two different promoters is
    retained twice: these are two gene-regulatory hypotheses, even when they
    later reuse one structural prediction.
    """
    from boltzscan.fimocistarget.cistarg_impl import filter_fimo_by_seq_and_overlap

    hit = pd.read_csv(fimo_csv)
    required = {'sequence_name', 'motif_id', 'start', 'stop', 'pvalue', 'score'}
    missing = sorted(required.difference(hit.columns))
    if missing:
        raise ValueError(
            f"FIMO table is missing required columns: {', '.join(missing)}"
        )
    print(f'Total hits: {hit.shape[0]}')
    hit['pvalue'] = pd.to_numeric(hit['pvalue'], errors='coerce')
    fimo_df = hit[hit['pvalue'].le(pvalue_thresh)].copy()
    print(f'Hits after filtering with pvalue <= {pvalue_thresh}: {fimo_df.shape[0]}')

    promoter_seqs = _read_fasta_to_str(promoter_fasta)
    pep_seqs = _read_fasta_to_str(tf_fasta)
    with open(tf2pwms_path) as handle:
        tf2pwms = json.load(handle)
    pwm2tf = _invert_tf2pwms(tf2pwms)

    fimo_df['motif_id'] = fimo_df['motif_id'].astype(str)
    fimo_df['promoter_seq'] = fimo_df['sequence_name'].map(promoter_seqs)
    missing_promoters = int(fimo_df['promoter_seq'].isna().sum())
    if missing_promoters:
        print(f'Dropping {missing_promoters} FIMO row(s) with no matching promoter FASTA record')
    fimo_df = fimo_df.dropna(subset=['promoter_seq'])
    fimo_df['tf_name'] = fimo_df['motif_id'].map(pwm2tf)
    fimo_df = fimo_df.explode('tf_name').dropna(subset=['tf_name'])
    fimo_df['tf_name'] = fimo_df['tf_name'].astype(str)
    fimo_df = fimo_df[fimo_df['tf_name'].isin(pep_seqs)].copy()
    print(f'Hits after restricting to TFs in {tf_fasta}: {fimo_df.shape[0]}')
    if fimo_df.empty:
        raise ValueError(
            'No FIMO hit remains after motif-to-TF mapping and TF FASTA restriction. '
            'Check tf2pwms IDs and the p-value threshold.'
        )

    # This function performs **only promoter-level** de-duplication.  Do not
    # introduce a global DNA drop_duplicates here: it would remove valid
    # promoter/gene edges before they can be mapped to a shared model input.
    fimo_df = filter_fimo_by_seq_and_overlap(
        fimo_df,
        overlap_threshold=overlap_thresh,
        extend_range=dna_flank,
    ).copy()
    print(f'Hits after promoter-level redundancy filtering: {fimo_df.shape[0]}')
    if fimo_df.empty:
        raise ValueError('No candidate remains after promoter-level redundancy filtering')

    fimo_df['extracted_motif_seq'] = fimo_df['extracted_motif_seq'].astype(str).str.upper()
    fimo_df['tf_seq'] = fimo_df['tf_name'].map(pep_seqs)
    fimo_df['boltz_name'] = (
        fimo_df['tf_name'] + '-' + fimo_df['motif_id'] + '-'
        + fimo_df['sequence_name'].astype(str) + '-' + fimo_df['start'].astype(str)
        + '-' + fimo_df['stop'].astype(str)
    )
    fimo_df['candidate_id'] = fimo_df.apply(_candidate_id, axis=1)
    sort_columns = _candidate_sort_columns(fimo_df)
    return fimo_df.sort_values(
        sort_columns,
        ascending=_candidate_sort_ascending(sort_columns),
        kind='mergesort',
    ).reset_index(drop=True)


def collapse_model_inputs(promoter_candidates):
    """Create reusable structural jobs without collapsing promoter edges.

    ``promoter_candidates`` must already have passed the promoter-level FIMO
    rule.  Every row is retained in the returned candidate and mapping tables.
    Only exact ``(tf_name, canonical_duplex(extracted_motif_seq))`` groups
    share a ``model_id`` and one prediction/YAML input.
    """
    required = {
        'candidate_id', 'tf_name', 'tf_seq', 'sequence_name', 'motif_id',
        'start', 'stop', 'extracted_motif_seq',
    }
    missing = sorted(required.difference(promoter_candidates.columns))
    if missing:
        raise ValueError(
            'Promoter candidate table is missing columns required for model-input '
            f"collapse: {', '.join(missing)}"
        )

    candidates = promoter_candidates.copy()
    candidates['canonical_dsdna'] = candidates['extracted_motif_seq'].map(canonical_duplex)
    candidates['model_id'] = [
        _stable_id('model', tf_name, dsdna)
        for tf_name, dsdna in zip(candidates['tf_name'], candidates['canonical_dsdna'])
    ]
    sort_columns = _candidate_sort_columns(candidates) + ['candidate_id']
    prioritised = candidates.sort_values(
        sort_columns,
        ascending=_candidate_sort_ascending(sort_columns[:-1]) + [True],
        kind='mergesort',
    )
    representatives = prioritised.drop_duplicates('model_id', keep='first').copy()
    representative_ids = set(representatives['candidate_id'])
    candidates['is_model_representative'] = candidates['candidate_id'].isin(representative_ids)

    group_counts = candidates.groupby('model_id', sort=False).agg(
        n_promoter_candidates=('candidate_id', 'size'),
        n_distinct_promoters=('sequence_name', 'nunique'),
    ).reset_index()
    model_inputs = representatives.loc[:, [
        'model_id', 'tf_name', 'tf_seq', 'canonical_dsdna', 'candidate_id',
    ]].rename(columns={'candidate_id': 'representative_candidate_id'})
    model_inputs = model_inputs.merge(group_counts, on='model_id', how='left', validate='one_to_one')
    model_inputs = model_inputs.sort_values('model_id', kind='mergesort').reset_index(drop=True)
    candidates = candidates.sort_values(
        sort_columns,
        ascending=_candidate_sort_ascending(sort_columns[:-1]) + [True],
        kind='mergesort',
    ).reset_index(drop=True)

    mapping_columns = [
        'candidate_id', 'model_id', 'is_model_representative', 'tf_name',
        'sequence_name', 'motif_id', 'start', 'stop', 'strand', 'pvalue',
        'score', 'extracted_motif_seq', 'canonical_dsdna',
    ]
    mapping_columns = [column for column in mapping_columns if column in candidates.columns]
    promoter_to_model = candidates.loc[:, mapping_columns].copy()
    return candidates, model_inputs, promoter_to_model


def _write_tf_dna_yaml(path, protein_sequence, dsdna_sequence, msa='empty'):
    """Write the fixed A/B/C TF--dsDNA YAML consumed by both predictors."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    dsdna_sequence = canonical_duplex(dsdna_sequence)
    payload = {
        'version': 1,
        'sequences': [
            {'protein': {'id': ['A'], 'sequence': str(protein_sequence), 'msa': str(msa)}},
            {'dna': {'id': ['B'], 'sequence': dsdna_sequence}},
            {'dna': {'id': ['C'], 'sequence': str(Seq(dsdna_sequence).reverse_complement())}},
        ],
    }
    path.write_text(yaml.safe_dump(payload, sort_keys=False))


def write_full_model_inputs(model_inputs, full_input_dir, msa_paths=None):
    """Write one full-length TF--dsDNA YAML for every model-level scan row."""
    required = {'model_id', 'tf_seq', 'canonical_dsdna'}
    missing = sorted(required.difference(model_inputs.columns))
    if missing:
        raise ValueError('Model-level scan table is missing: ' + ', '.join(missing))
    full_input_dir = Path(full_input_dir)
    full_input_dir.mkdir(parents=True, exist_ok=True)
    for row in model_inputs.itertuples(index=False):
        msa = 'empty' if msa_paths is None else msa_paths.get(row.tf_name)
        if msa is None:
            raise FileNotFoundError(f'No MSA supplied for TF {row.tf_name}')
        _write_tf_dna_yaml(
            full_input_dir / f'{row.model_id}.yaml',
            protein_sequence=row.tf_seq,
            dsdna_sequence=row.canonical_dsdna,
            msa=msa,
        )
    return full_input_dir


def write_coordinate_cropped_model_inputs(
    model_inputs,
    intervals_by_tf,
    crop_input_dir,
    crop_manifest_path,
    flank,
    msa_paths=None,
    boundary_source='interface',
    evidence_by_tf=None,
):
    """Write YAML/A3M inputs from one 1-based inclusive interval per TF."""
    if flank < 0:
        raise ValueError('Crop flank must be non-negative')
    if not {'model_id', 'tf_name', 'tf_seq', 'canonical_dsdna'}.issubset(model_inputs.columns):
        raise ValueError('model_inputs must contain model_id, tf_name, tf_seq, canonical_dsdna')

    crop_input_dir = Path(crop_input_dir)
    crop_manifest_path = Path(crop_manifest_path)
    crop_input_dir.mkdir(parents=True, exist_ok=True)
    expected_yamls = {
        f'{model_id}.yaml' for model_id in model_inputs['model_id'].astype(str)
    }
    stale_yamls = [
        path for path in crop_input_dir.glob('*.yaml') if path.name not in expected_yamls
    ]
    if stale_yamls:
        raise FileExistsError(
            f'{crop_input_dir} contains {len(stale_yamls)} stale crop YAML file(s); '
            'use a new run directory for a changed candidate universe'
        )
    required_tfs = sorted(set(model_inputs['tf_name']))
    missing_tfs = [tf_name for tf_name in required_tfs if tf_name not in intervals_by_tf]
    if missing_tfs:
        preview = ', '.join(missing_tfs[:10])
        raise ValueError(
            f'No {boundary_source} interval for {len(missing_tfs)} TF(s): {preview}. '
            'Crop inputs were not written.'
        )

    if msa_paths is not None:
        from boltzscan.msa import crop_a3m_file
        from boltzscan.msa import read_a3m_query

    crop_rules = {}
    for tf_name, group in model_inputs.groupby('tf_name', sort=False):
        protein_sequences = set(group['tf_seq'])
        if len(protein_sequences) != 1:
            raise ValueError(f'TF {tf_name!r} has inconsistent protein sequences in model_inputs')
        protein_sequence = protein_sequences.pop()
        boundary_lo, boundary_hi = map(int, intervals_by_tf[tf_name])
        if boundary_lo < 1 or boundary_hi < boundary_lo or boundary_hi > len(protein_sequence):
            raise ValueError(
                f'Invalid {boundary_source} interval {boundary_lo}-{boundary_hi} '
                f'for {tf_name} length {len(protein_sequence)}'
            )
        crop_lo = max(1, boundary_lo - flank)
        crop_hi = min(len(protein_sequence), boundary_hi + flank)
        crop_msa = None
        if msa_paths is not None:
            source_msa = msa_paths.get(tf_name)
            if source_msa is None:
                raise FileNotFoundError(f'No full-length MSA supplied for TF {tf_name}')
            source_msa = Path(source_msa)
            crop_msa = crop_input_dir / 'msa' / f'{_stable_id("tf", tf_name)}.a3m'
            crop_a3m_file(
                source_msa,
                crop_msa,
                interval=(crop_lo - 1, crop_hi),
            )
            expected_query = protein_sequence[crop_lo - 1:crop_hi].upper()
            if read_a3m_query(crop_msa) != expected_query:
                crop_msa.unlink(missing_ok=True)
                raise ValueError(f'Cropped MSA query mismatch for TF {tf_name}')
            crop_msa = crop_msa.resolve()
        crop_rules[tf_name] = (crop_lo, crop_hi, boundary_lo, boundary_hi, crop_msa)

    manifest_rows = []
    for row in model_inputs.itertuples(index=False):
        crop_lo, crop_hi, boundary_lo, boundary_hi, crop_msa = crop_rules[row.tf_name]
        crop_sequence = row.tf_seq[crop_lo - 1:crop_hi]
        crop_yaml = crop_input_dir / f'{row.model_id}.yaml'
        _write_tf_dna_yaml(
            crop_yaml,
            crop_sequence,
            row.canonical_dsdna,
            msa=crop_msa or 'empty',
        )
        manifest_rows.append({
            'model_id': row.model_id,
            'tf_name': row.tf_name,
            'boundary_source': boundary_source,
            'boundary_start_1based_inclusive': boundary_lo,
            'boundary_stop_1based_inclusive': boundary_hi,
            'boundary_evidence': str((evidence_by_tf or {}).get(row.tf_name, '')),
            'protein_length': len(row.tf_seq),
            'crop_start_1based_inclusive': crop_lo,
            'crop_stop_1based_inclusive': crop_hi,
            'crop_length': len(crop_sequence),
            'crop_flank_aa': flank,
            'crop_input_yaml': str(crop_yaml),
            'crop_msa_path': str(crop_msa or ''),
        })
    manifest = pd.DataFrame(manifest_rows).sort_values('model_id', kind='mergesort')
    crop_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(crop_manifest_path, index=False)
    return CropInputSummary(
        crop_input_dir=crop_input_dir,
        crop_manifest=crop_manifest_path,
        n_model_inputs=len(manifest),
    )


def write_dbd_cropped_model_inputs(
    model_inputs,
    domtbl_path,
    crop_input_dir,
    crop_manifest_path,
    flank,
    msa_paths=None,
):
    """Legacy/internal HMM boundary adapter for coordinate-driven cropping."""
    from boltzscan.pwmmap.dbd import parse_domtbl_dbds

    records_by_tf = {}
    for record in parse_domtbl_dbds(domtbl_path):
        records_by_tf.setdefault(record[0], []).append(record)
    intervals = {
        tf_name: (
            min(record[2] for record in records),
            max(record[3] for record in records),
        )
        for tf_name, records in records_by_tf.items()
    }
    evidence = {
        tf_name: ';'.join(
            f'{pfam}:{start}-{stop}' for _, pfam, start, stop in records
        )
        for tf_name, records in records_by_tf.items()
    }
    return write_coordinate_cropped_model_inputs(
        model_inputs,
        intervals_by_tf=intervals,
        crop_input_dir=crop_input_dir,
        crop_manifest_path=crop_manifest_path,
        flank=flank,
        msa_paths=msa_paths,
        boundary_source='hmmsearch_dbd',
        evidence_by_tf=evidence,
    )


def write_fimo_model_yamls(model_table, msa_dir, domtbl_path, output_dir, crop=None):
    """Write full MSA-backed YAMLs; DBD crop arguments are internal legacy API."""
    from boltzscan.msa import read_a3m_query

    model_table = Path(model_table)
    if not model_table.is_file():
        raise FileNotFoundError(f'Model-level FIMO table not found: {model_table}')
    model_inputs = pd.read_csv(model_table)
    required = {'model_id', 'tf_name', 'tf_seq', 'canonical_dsdna'}
    missing = sorted(required.difference(model_inputs.columns))
    if missing:
        raise ValueError(
            'Model-level FIMO table is missing: ' + ', '.join(missing)
            + '. Use filtered_model_scan_level_res.csv from fimo-scan.'
        )
    if model_inputs.empty:
        raise ValueError(f'Model-level FIMO table is empty: {model_table}')
    if model_inputs[list(required)].isna().any().any():
        raise ValueError('Model-level FIMO table contains empty required values')
    if model_inputs['model_id'].duplicated().any():
        raise ValueError('Model-level FIMO table contains duplicate model_id values')
    if crop is not None:
        if crop < 0:
            raise ValueError('--crop must be a non-negative amino-acid flank')
        if domtbl_path is None:
            raise ValueError('--domtbl is required when --crop is enabled')

    msa_paths = _build_msa_dct(msa_dir)
    tf_msas = {}
    for tf_name, group in model_inputs.groupby('tf_name', sort=False):
        sequences = set(group['tf_seq'].astype(str))
        if len(sequences) != 1:
            raise ValueError(f'TF {tf_name!r} has inconsistent protein sequences')
        protein_sequence = sequences.pop().upper()
        msa_path = _resolve_msa_path(str(tf_name), msa_paths)
        if msa_path is None:
            raise FileNotFoundError(f'MSA not found under {msa_dir} for TF {tf_name}')
        msa_path = Path(msa_path).resolve()
        if read_a3m_query(msa_path) != protein_sequence:
            raise ValueError(f'Protein/MSA query mismatch for TF {tf_name}: {msa_path}')
        tf_msas[tf_name] = msa_path

    output_dir = Path(output_dir)
    full_input_dir = output_dir / 'full'
    expected = {f'{model_id}.yaml' for model_id in model_inputs['model_id'].astype(str)}
    for model_id in model_inputs['model_id'].astype(str):
        if Path(model_id).name != model_id:
            raise ValueError(f'Unsafe model_id in model-level FIMO table: {model_id!r}')
    input_dirs = [full_input_dir]
    if crop is not None:
        input_dirs.append(output_dir / f'crop{crop}')
    for input_dir in input_dirs:
        input_dir.mkdir(parents=True, exist_ok=True)
        stale = [path for path in input_dir.glob('*.yaml') if path.name not in expected]
        if stale:
            raise FileExistsError(
                f'{input_dir} contains {len(stale)} stale YAML file(s); use an empty directory'
            )

    write_full_model_inputs(model_inputs, full_input_dir, msa_paths=tf_msas)
    if crop is None:
        return full_input_dir, None

    crop_input_dir = output_dir / f'crop{crop}'
    crop_summary = write_dbd_cropped_model_inputs(
        model_inputs,
        domtbl_path=domtbl_path,
        crop_input_dir=crop_input_dir,
        crop_manifest_path=crop_input_dir / 'crop_manifest.csv',
        flank=crop,
        msa_paths=tf_msas,
    )
    return full_input_dir, crop_summary
