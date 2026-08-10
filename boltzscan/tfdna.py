"""The concise, auditable TF--DNA workflow used by ``boltzscan run``.

This module is intentionally a thin orchestration layer.  The scientific
rules live with their data transformations:

* promoter-local FIMO de-duplication and model-input collapse live in
  :mod:`boltzscan.utils.boltz_input`;
* reusable structure-interface cropping lives in :mod:`boltzscan.interface`.

Keeping those rules out of task-local notebooks/scripts makes the two output
levels reproducible across species and experiments.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


_STAGES = {'all', 'prepare', 'predict', 'score'}
_RUN_CONFIG = 'run_config.json'


@dataclass(frozen=True)
class TfDnaRunSummary:
    """Locations produced or reused by a run-level TF--DNA workflow."""

    out_dir: Path
    candidates_dir: Path
    full_input_dir: Path
    crop_input_dir: Path | None
    predictions_dir: Path
    results_dir: Path


def _selected_stages(stage):
    if stage not in _STAGES:
        raise ValueError(f"Unknown run stage: {stage}")
    return ('prepare', 'predict', 'score') if stage == 'all' else (stage,)


def _read_run_config(out_dir):
    path = Path(out_dir) / _RUN_CONFIG
    if not path.is_file():
        raise FileNotFoundError(
            f"Run configuration not found: {path}. Run prepare first."
        )
    try:
        config = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid run configuration: {path}") from exc
    if config.get('schema_version') != 1:
        raise ValueError(f"Unsupported run configuration schema: {path}")
    return config


def _write_run_config(out_dir, *, model, crop, seed):
    path = Path(out_dir) / _RUN_CONFIG
    staged = path.with_suffix('.json.tmp')
    staged.write_text(json.dumps({
        'schema_version': 1,
        'model': model,
        'crop': crop,
        'seed': seed,
    }, indent=2) + '\n')
    staged.replace(path)
    return path


def _resolve_run_settings(out_dir, *, prepare, reuse, model, crop, seed):
    config_path = Path(out_dir) / _RUN_CONFIG
    use_existing = not prepare or (reuse and config_path.is_file())
    if not use_existing:
        return model or 'esmfold2', crop, seed

    config = _read_run_config(out_dir)
    requested = {'model': model, 'crop': crop, 'seed': seed}
    for name, value in requested.items():
        if value is not None and value != config.get(name):
            raise ValueError(
                f"--{name}={value} conflicts with {config_path.name}: "
                f"{config.get(name)}"
            )
    return config['model'], config.get('crop'), config.get('seed')


def _normalise_excluded_species(values):
    result = []
    for value in values or ():
        result.extend(part.strip() for part in str(value).split(',') if part.strip())
    return tuple(dict.fromkeys(result))


def _refuse_existing(path, label, reuse):
    if path.exists() and any(path.iterdir()) and not reuse:
        raise FileExistsError(
            f"{label} already exists at {path}. Pass --resume to use its existing outputs "
            "or choose a new --run directory."
        )


def _resolve_pwm_inputs(
    *,
    tf_fasta,
    out_dir,
    refs,
    exclude_species,
    collapse_clusters,
    cpu,
    reuse,
    pfam=None,
):
    """Map the run TFs and return its motif, mapping, and DBD paths."""
    from boltzscan.pwmmap.mapper import map_species

    out_dir = Path(out_dir)
    pwm_dir = out_dir / 'pwm'
    _refuse_existing(pwm_dir, 'PWM mapping directory', reuse)
    refs_for_mapping = Path(refs)
    if exclude_species:
        from boltzscan.pwmmap.leaveout import create_reference_subset

        subset_dir = out_dir / 'refs_loso'
        if subset_dir.exists():
            if not reuse:
                raise FileExistsError(
                    f"LOSO reference subset already exists at {subset_dir}. Pass --resume "
                    "to use this immutable subset or choose a new --run directory."
                )
        else:
            summary = create_reference_subset(
                refs_for_mapping,
                subset_dir,
                exclude_species=list(exclude_species),
            )
            print(
                f"LOSO reference subset: retained {summary.n_retained_dbd_rows}/"
                f"{summary.n_input_dbd_rows} DBD records after excluding "
                f"{', '.join(exclude_species)}"
            )
        refs_for_mapping = subset_dir

    mapping_json = pwm_dir / 'tf2pwms.json'
    if mapping_json.exists() and reuse:
        return pwm_dir / 'meme', mapping_json, pwm_dir / 'pfam.domtbl'

    # Cluster representatives are optional; when requested after a LOSO
    # subset, cluster only the filtered store so target-species motifs cannot
    # influence representative choice.
    if collapse_clusters and exclude_species:
        from boltzscan.pwmmap.cluster import cluster_reference_motifs

        cluster_reference_motifs(refs_dir=refs_for_mapping, cpu=cpu)
    summary = map_species(
        species_fasta=tf_fasta,
        out_dir=pwm_dir,
        refs_dir=refs_for_mapping,
        domtbl=None,
        pfam=pfam,
        cpu=cpu,
        collapse_clusters=collapse_clusters,
    )
    print(
        f"PWM transfer: {summary.n_mapped}/{summary.n_species_tfs} TFs mapped; "
        f"{summary.n_motifs} scan-ready motif(s)"
    )
    return pwm_dir / 'meme', mapping_json, pwm_dir / 'pfam.domtbl'


def _run_prepare_stage(
    *,
    tf_fasta,
    promoters,
    out_dir,
    refs,
    exclude_species,
    collapse_clusters,
    pvalue,
    overlap_thresh,
    dna_flank,
    crop,
    model,
    reuse,
    pfam=None,
):
    """Map PWMs, scan promoters, and materialize model inputs."""
    from boltzscan.fimocistarget.runner import run_fimo_scan
    from boltzscan.utils.boltz_input import (
        write_fimo_model_yamls,
        write_full_model_inputs,
    )

    out_dir = Path(out_dir)
    from boltzscan.predict.runners import available_cpu_count
    from boltzscan.runlog import CommandLog

    cpu = available_cpu_count()

    with CommandLog(
        out_dir / 'map-pwm.log',
        ['boltzscan', 'run', 'step=pwm-map'],
        SimpleNamespace(
            refs=refs,
            pfam=pfam,
            exclude_species=list(exclude_species),
            pwm_clustered=collapse_clusters,
            cpu=cpu,
            reuse=reuse,
        ),
    ):
        motif_dir, tf2pwms_path, _resolved_domtbl = _resolve_pwm_inputs(
            tf_fasta=tf_fasta,
            out_dir=out_dir,
            refs=refs,
            exclude_species=exclude_species,
            collapse_clusters=collapse_clusters,
            cpu=cpu,
            reuse=reuse,
            pfam=pfam,
        )

    scan_dir = out_dir / 'scan'
    if not reuse:
        _refuse_existing(scan_dir, 'FIMO scan directory', reuse)
    with CommandLog(
        out_dir / 'fimo-scan.log',
        ['boltzscan', 'run', 'step=fimo-scan'],
        SimpleNamespace(
            promoters=promoters,
            motif_dir=motif_dir,
            tf2pwms=tf2pwms_path,
            tf_fasta=tf_fasta,
            pvalue=pvalue,
            overlap_thresh=overlap_thresh,
            dna_flank=dna_flank,
            pwm_clustered=collapse_clusters,
            reuse=reuse,
        ),
    ):
        scan_summary = run_fimo_scan(
            target_fasta=promoters,
            motif_dir=motif_dir,
            tf2pwms_path=tf2pwms_path,
            tf_fasta=tf_fasta,
            output_dir=scan_dir,
            pvalue_thresh=pvalue,
            overlap_thresh=overlap_thresh,
            dna_flank=dna_flank,
            no_pwm_cluster=not collapse_clusters,
            reuse=reuse,
        )

    from boltzscan.predict.runners import engine_for_model

    msa_dir = out_dir / 'msa'
    requires_msa = engine_for_model(model) == 'boltz'
    if requires_msa:
        from boltzscan.msa import DEFAULT_MSA_SERVER_URL, run_msa

        with CommandLog(
            out_dir / 'msa.log',
            ['boltzscan', 'run', 'step=msa'],
            SimpleNamespace(
                fasta=tf_fasta,
                output=msa_dir,
                model=model,
                server_url=DEFAULT_MSA_SERVER_URL,
            ),
        ):
            msa_summary = run_msa(
                tf_fasta,
                msa_dir,
                server_url=DEFAULT_MSA_SERVER_URL,
            )
            print(
                f"MSA: requested={msa_summary.requested}, "
                f"generated={msa_summary.completed}, reused={msa_summary.skipped}"
            )

    input_arm = 'interface' if crop is not None else 'full'
    full_input_dir = out_dir / 'inputs' / input_arm
    with CommandLog(
        out_dir / 'fimo2yaml.log',
        ['boltzscan', 'run', 'step=model-inputs'],
        SimpleNamespace(
            model_table=scan_summary.model_level,
            crop=crop,
            msa='required' if requires_msa else 'empty',
        ),
    ):
        model_inputs = pd.read_csv(scan_summary.model_level)
        if crop is not None:
            from boltzscan.interface import write_interface_localizer_inputs

            write_interface_localizer_inputs(
                model_inputs,
                full_input_dir,
                msa_dir=msa_dir if requires_msa else None,
            )
            crop_summary = None
        elif requires_msa:
            full_input_dir, crop_summary = write_fimo_model_yamls(
                scan_summary.model_level,
                msa_dir,
                None,
                out_dir / 'inputs',
                crop=None,
            )
        else:
            write_full_model_inputs(model_inputs, full_input_dir)
            crop_summary = None
        message = (
            f"Model inputs: {scan_summary.n_promoter_candidates} promoter candidates -> "
            f"{scan_summary.n_model_inputs} reusable TF+dsDNA predictions"
        )
        if crop is not None:
            message += f"; {len(list(full_input_dir.glob('*.yaml')))} full interface localizer(s)"
        print(message)


def _predict_arm(*, input_dir, run_dir, model, seed):
    from boltzscan.predict.runners import run_prediction

    if not input_dir.is_dir() or not any(input_dir.glob('*.yaml')):
        raise FileNotFoundError(f"No prepared YAML inputs found for prediction: {input_dir}")
    return run_prediction(
        input_dir=input_dir,
        out_dir=run_dir,
        model=model,
        seed=seed,
    )


def _run_predict_stage(*, out_dir, crop, model, seed):
    from boltzscan.predict.wash import wash_prediction_outputs
    from boltzscan.predict.runners import native_prediction_dir
    from boltzscan.runlog import CommandLog

    out_dir = Path(out_dir)
    wash_mode = 'soft'
    if crop is None:
        arms = prediction_arms = ['full']
    else:
        from boltzscan.interface import run_interface_stage

        with CommandLog(
            out_dir / 'find-interface.log',
            ['boltzscan', 'run', 'step=find-interface'],
            SimpleNamespace(model=model, flank=crop, cutoff=5.0, seed=seed),
        ):
            interface = run_interface_stage(
                out_dir,
                flank=crop,
                cutoff=5.0,
                model=model,
                seed=seed,
            )
            print(
                f"Interfaces: {interface.n_tfs} TF(s) -> "
                f"{interface.n_model_inputs} crop inputs; {interface.boundaries}"
            )
        arms = ['interface', f'crop{crop}']
        prediction_arms = [f'crop{crop}']
    for arm in prediction_arms:
        input_dir = out_dir / 'inputs' / arm
        with CommandLog(
            out_dir / f'predict-{arm}.log',
            ['boltzscan', 'run', f'step=predict-{arm}'],
            SimpleNamespace(
                model=model,
                input_dir=input_dir,
                configuration=(
                    f'boltzscan_{model}' if model in {'boltz1_ode', 'boltz2_ode'}
                    else 'native_default'
                ),
            ),
        ):
            result = _predict_arm(
                input_dir=input_dir,
                run_dir=out_dir,
                model=model,
                seed=seed,
            )
            if result != 0:
                raise RuntimeError(f"{model} prediction failed for {arm} with exit code {result}")
            native_output = native_prediction_dir(input_dir, out_dir, model)
            print(f"{model} native prediction complete: {arm} -> {native_output}")

    with CommandLog(
        out_dir / 'wash.log',
        ['boltzscan', 'run', 'step=prediction-wash'],
        SimpleNamespace(model=model, mode=wash_mode),
    ):
        washed = wash_prediction_outputs(
            out_dir,
            model=model,
            mode=wash_mode,
            arms=arms,
        )
        print(f"Prediction wash ({wash_mode}) -> {washed.method_root}")


def _score_prediction_root(out_dir, model):
    from boltzscan.predict.runners import prediction_root_name

    return Path(out_dir) / prediction_root_name(model)


def score_tf_dna_run(*, out_dir, model=None, processes=1):
    """Score every available structural arm and project scores to promoter rows."""
    from boltzscan.predict.runners import PREDICTION_MODELS, prediction_root_name

    out_dir = Path(out_dir)
    if not out_dir.is_dir():
        raise FileNotFoundError(f'Run directory not found: {out_dir}')
    if model is None:
        models = [
            name for name in PREDICTION_MODELS
            if (out_dir / prediction_root_name(name)).is_dir()
        ]
        if len(models) != 1:
            found = ', '.join(models) if models else 'none'
            raise ValueError(
                f'Could not select one prediction model in {out_dir} (found: {found}). '
                'Run `boltzscan wash` first or pass --model.'
            )
        model = models[0]

    prediction_root = _score_prediction_root(out_dir, model)
    if not prediction_root.is_dir():
        raise FileNotFoundError(f'Prediction directory not found: {prediction_root}')
    crop_arms = sorted(
        (
            child.name for child in prediction_root.iterdir()
            if child.is_dir() and child.name.startswith('crop') and child.name[4:].isdigit()
        ),
        key=lambda name: int(name[4:]),
    )
    arms = (['full'] if (prediction_root / 'full').is_dir() else []) + crop_arms
    if not arms:
        raise FileNotFoundError(
            f'No full or crop<N> prediction arm found in {prediction_root}'
        )
    print(
        f'ipSAE protocol: model={model}, arms={"/".join(arms)}, '
        f'PAE cutoff=10, distance cutoff=10 A, CPU workers={processes}'
    )
    _run_score_stage_impl(
        out_dir=out_dir,
        crop=int(crop_arms[0][4:]) if len(crop_arms) == 1 else None,
        processes=processes,
        model=model,
        arms=arms,
    )
    return model, tuple(arms)


def _run_score_stage_impl(*, out_dir, crop, processes, model='esmfold2', arms=None):
    from boltzscan.utils.ipsae_score import print_ipsae_warnings, score_ipsae_table

    out_dir = Path(out_dir)
    prediction_root = _score_prediction_root(out_dir, model)
    scan_dir = out_dir / 'scan'
    model_inputs_path = scan_dir / 'filtered_model_scan_level_res.csv'
    mapping_path = scan_dir / 'filtered_promoter_scan_level_res.csv'
    if not model_inputs_path.is_file() or not mapping_path.is_file():
        raise FileNotFoundError(
            'Model-input mapping tables are missing. Run `boltzscan run --stage prepare` first.'
        )
    results_dir = out_dir / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(mapping_path)
    if arms is None:
        arms = ('full',) if crop is None else (f'crop{crop}',)
    for arm in arms:
        prediction_dir = prediction_root / arm
        model_score_path = results_dir / f'{arm}_model_input_ipsae.csv.gz'
        summary = score_ipsae_table(
            res_dir=prediction_dir,
            score_file=model_inputs_path,
            output=model_score_path,
            id_col='model_id',
            processes=processes,
        )
        print_ipsae_warnings(summary.warnings)
        model_scores = pd.read_csv(model_score_path)
        score_columns = [
            'model_id', 'ipsae', 'ipsae_asym_min', 'boltz_iptm', 'ipsae_iptm',
            'iptm_diff', 'boltz_iptm_global', 'boltz_iptm_asym_min',
            'ipsae_iptm_asym_min', 'pDockQ',
        ]
        score_columns = [column for column in score_columns if column in model_scores.columns]
        promoter_scores = mapping.merge(
            model_scores.loc[:, score_columns],
            on='model_id',
            how='left',
            validate='many_to_one',
        )
        promoter_scores['TF'] = promoter_scores['tf_name']
        promoter_scores['TG'] = promoter_scores['sequence_name']
        promoter_scores = promoter_scores.sort_values(
            ['ipsae', 'pvalue', 'candidate_id'],
            ascending=[False, True, True],
            na_position='last',
            kind='mergesort',
        )
        promoter_score_path = results_dir / f'{arm}_promoter_level_ipsae.csv.gz'
        promoter_scores.to_csv(promoter_score_path, index=False, compression='gzip')
        print(
            f"ipSAE: {summary.scored_predictions}/{summary.total_predictions} model inputs scored; "
            f"{len(promoter_scores)} promoter candidates written -> {promoter_score_path}"
        )

    promoter_candidates_path = scan_dir / 'filtered_promoter_scan_level_res.csv'
    crop_arms = [
        arm for arm in arms
        if arm.startswith('crop') and arm[4:].isdigit()
    ]
    if promoter_candidates_path.is_file() and 'full' in arms and len(crop_arms) == 1:
        from boltzscan.ranking import write_three_arm_score_tables

        ranking = write_three_arm_score_tables(
            out_dir,
            crop=int(crop_arms[0][4:]),
            model=model,
        )
        print(
            f"Three-arm rankings: {ranking.n_candidates} promoter candidates, "
            f"{ranking.n_model_inputs} model inputs -> {ranking.comparison}"
        )


def _run_score_stage(*, out_dir, crop, processes, model='esmfold2'):
    from boltzscan.runlog import CommandLog

    out_dir = Path(out_dir)
    with CommandLog(
        out_dir / 'ipsae.log',
        ['boltzscan', 'run', 'step=score'],
        SimpleNamespace(model=model, crop=crop, processes=processes),
    ):
        return score_tf_dna_run(
            out_dir=out_dir,
            processes=processes,
            model=model,
        )


def run_tf_dna_workflow(
    *,
    tf_fasta,
    promoters,
    out_dir,
    refs='data/pwms/_refs',
    pfam=None,
    exclude_species=(),
    collapse_clusters=True,
    pvalue=1e-5,
    overlap_thresh=0.0,
    dna_flank=5,
    crop=None,
    model=None,
    seed=None,
    stage='all',
    reuse=False,
):
    """Run requested parts of the concise TF--DNA workflow.

    ``prepare`` always includes PWM transfer + FIMO because those steps define
    the candidate universe.  ``predict`` and ``score`` intentionally consume
    existing inputs/results when requested alone, making it safe to inspect
    the candidate tables before submitting GPU work.
    """
    out_dir = Path(out_dir)
    stages = _selected_stages(stage)
    out_dir.mkdir(parents=True, exist_ok=True)
    model, crop, seed = _resolve_run_settings(
        out_dir,
        prepare='prepare' in stages,
        reuse=reuse,
        model=model,
        crop=crop,
        seed=seed,
    )
    if 'prepare' in stages:
        if tf_fasta is None:
            raise ValueError('--tf is required when prepare is selected')
        if promoters is None:
            raise ValueError('--promoters is required when prepare is selected')
        tf_fasta = Path(tf_fasta)
        promoters = Path(promoters)
        if not tf_fasta.is_file():
            raise FileNotFoundError(f"TF FASTA not found: {tf_fasta}")
        if not promoters.is_file():
            raise FileNotFoundError(f"Promoter FASTA not found: {promoters}")
    if crop is not None and crop < 0:
        raise ValueError('--crop must be a non-negative amino-acid flank')
    if dna_flank < 0:
        raise ValueError('--dna-flank must be non-negative')
    from boltzscan.predict.runners import PREDICTION_MODELS

    if model not in PREDICTION_MODELS:
        raise ValueError(f"Unknown prediction model: {model}")
    exclude_species = _normalise_excluded_species(exclude_species)

    if 'prepare' in stages:
        _run_prepare_stage(
            tf_fasta=tf_fasta,
            promoters=promoters,
            out_dir=out_dir,
            refs=refs,
            pfam=pfam,
            exclude_species=exclude_species,
            collapse_clusters=collapse_clusters,
            pvalue=pvalue,
            overlap_thresh=overlap_thresh,
            dna_flank=dna_flank,
            crop=crop,
            model=model,
            reuse=reuse,
        )
        config_path = _write_run_config(
            out_dir, model=model, crop=crop, seed=seed
        )
        print(f"Run configuration: {config_path}")
    if 'predict' in stages:
        _run_predict_stage(
            out_dir=out_dir,
            crop=crop,
            model=model,
            seed=seed,
        )
    if 'score' in stages:
        from boltzscan.predict.runners import available_cpu_count

        _run_score_stage(
            out_dir=out_dir,
            crop=crop,
            processes=available_cpu_count(),
            model=model,
        )

    from boltzscan.predict.runners import prediction_root_name

    return TfDnaRunSummary(
        out_dir=out_dir,
        candidates_dir=out_dir / 'scan',
        full_input_dir=out_dir / 'inputs' / ('interface' if crop is not None else 'full'),
        crop_input_dir=(out_dir / 'inputs' / f'crop{crop}') if crop is not None else None,
        predictions_dir=out_dir / prediction_root_name(model),
        results_dir=out_dir / 'results',
    )
