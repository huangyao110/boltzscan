"""Create directly comparable FIMO/full/crop promoter-level score tables."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd


_STRUCTURE_METRICS = (
    'ipsae',
    'ipsae_asym_min',
    'boltz_iptm',
    'ipsae_iptm',
    'iptm_diff',
    'boltz_iptm_global',
    'boltz_iptm_asym_min',
    'ipsae_iptm_asym_min',
    'pDockQ',
)


@dataclass(frozen=True)
class ThreeArmRankingSummary:
    """Files and row counts for one three-arm comparison."""

    fimo: Path
    full: Path
    crop: Path
    comparison: Path
    metadata: Path
    n_candidates: int
    n_model_inputs: int
    n_target_promoters: int


def _require_table(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f'Missing {label}: {path}')
    table = pd.read_csv(path)
    if table.empty:
        raise ValueError(f'{label} is empty: {path}')
    if 'candidate_id' not in table.columns:
        raise ValueError(f'{label} has no candidate_id column: {path}')
    if table['candidate_id'].duplicated().any():
        raise ValueError(f'{label} contains duplicate candidate_id values: {path}')
    return table


def _sequential_rank(table: pd.DataFrame, columns, ascending) -> pd.Series:
    ordered = table.sort_values(
        columns,
        ascending=ascending,
        na_position='last',
        kind='mergesort',
    )
    ranks = pd.Series(range(1, len(ordered) + 1), index=ordered.index, dtype='int64')
    return ranks.reindex(table.index)


def _target_site_rank(table: pd.DataFrame, rank_column: str) -> pd.Series:
    ordered = table.sort_values(rank_column, kind='mergesort')
    ranks = ordered.groupby('sequence_name', sort=False).cumcount().add(1)
    ranks.index = ordered.index
    return ranks.reindex(table.index).astype('int64')


def _prefixed_metrics(table: pd.DataFrame, prefix: str) -> pd.DataFrame:
    columns = ['candidate_id'] + [
        metric for metric in _STRUCTURE_METRICS if metric in table.columns
    ]
    return table.loc[:, columns].rename(columns={
        metric: f'{prefix}_{metric}' for metric in columns if metric != 'candidate_id'
    })


def _validate_candidate_universe(
    candidates: pd.DataFrame,
    full: pd.DataFrame,
    crop: pd.DataFrame,
) -> None:
    expected = set(candidates['candidate_id'].astype(str))
    for label, table in (('full', full), ('crop', crop)):
        observed = set(table['candidate_id'].astype(str))
        if observed != expected:
            raise ValueError(
                f'{label} promoter score candidate universe differs from FIMO: '
                f'{len(expected - observed)} missing, {len(observed - expected)} extra'
            )


def write_three_arm_score_tables(
    run_dir: str | Path,
    *,
    crop: int = 20,
    model: str = 'esmfold2',
) -> ThreeArmRankingSummary:
    """Write FIMO, FIMO+full, and FIMO+crop promoter rankings.

    FIMO defines the candidate gate for all three arms.  The sequence-only arm
    is ranked by p-value (then PWM score), while the two shape-readout arms are
    ranked by their respective ipSAE values with FIMO statistics used only as
    deterministic tie-breakers.  No uncalibrated weighted composite is used.
    """

    run_dir = Path(run_dir)
    candidate_path = run_dir / 'scan' / 'filtered_promoter_scan_level_res.csv'
    if not candidate_path.is_file():
        candidate_path = run_dir / 'candidates' / 'promoter_candidates.csv.gz'
    candidates = _require_table(
        candidate_path,
        'promoter candidate table',
    )
    full = _require_table(
        run_dir / 'results' / 'full_promoter_level_ipsae.csv.gz',
        'full promoter score table',
    )
    crop_arm = f'crop{crop}'
    crop_scores = _require_table(
        run_dir / 'results' / f'{crop_arm}_promoter_level_ipsae.csv.gz',
        f'{crop_arm} promoter score table',
    )
    _validate_candidate_universe(candidates, full, crop_scores)

    required = {
        'candidate_id', 'model_id', 'tf_name', 'sequence_name', 'motif_id',
        'start', 'stop', 'strand', 'score', 'pvalue', 'extracted_motif_seq',
        'canonical_dsdna',
    }
    missing = sorted(required.difference(candidates.columns))
    if missing:
        raise ValueError('Promoter candidate table is missing: ' + ', '.join(missing))

    base_columns = [
        'candidate_id', 'model_id', 'tf_name', 'sequence_name', 'motif_id',
        'start', 'stop', 'strand', 'extracted_motif_seq', 'canonical_dsdna',
        'is_model_representative', 'score', 'pvalue', 'qvalue', 'real_start',
        'real_stop', 'length', 'boltz_name',
    ]
    base_columns = [column for column in base_columns if column in candidates.columns]
    comparison = candidates.loc[:, base_columns].copy().rename(columns={
        'score': 'fimo_pwm_score',
        'pvalue': 'fimo_pvalue',
        'qvalue': 'fimo_qvalue',
    })
    pvalues = pd.to_numeric(comparison['fimo_pvalue'], errors='coerce')
    if pvalues.isna().any() or pvalues.le(0).any():
        raise ValueError('FIMO pvalue must be finite and > 0 for every candidate')
    comparison['fimo_neglog10_pvalue'] = -np.log10(pvalues)

    comparison = comparison.merge(
        _prefixed_metrics(full, 'full'),
        on='candidate_id',
        how='left',
        validate='one_to_one',
    ).merge(
        _prefixed_metrics(crop_scores, crop_arm),
        on='candidate_id',
        how='left',
        validate='one_to_one',
    )
    biological_tie_breakers = [
        'sequence_name', 'start', 'stop', 'strand', 'motif_id', 'candidate_id',
    ]
    comparison['fimo_rank'] = _sequential_rank(
        comparison,
        ['fimo_pvalue', 'fimo_pwm_score'] + biological_tie_breakers,
        [True, False] + [True] * len(biological_tie_breakers),
    )
    comparison['full_rank'] = _sequential_rank(
        comparison,
        ['full_ipsae', 'fimo_pvalue', 'fimo_pwm_score'] + biological_tie_breakers,
        [False, True, False] + [True] * len(biological_tie_breakers),
    )
    crop_ipsae = f'{crop_arm}_ipsae'
    crop_rank = f'{crop_arm}_rank'
    comparison[crop_rank] = _sequential_rank(
        comparison,
        [crop_ipsae, 'fimo_pvalue', 'fimo_pwm_score'] + biological_tie_breakers,
        [False, True, False] + [True] * len(biological_tie_breakers),
    )
    comparison['fimo_target_site_rank'] = _target_site_rank(comparison, 'fimo_rank')
    comparison['full_target_site_rank'] = _target_site_rank(comparison, 'full_rank')
    comparison[f'{crop_arm}_target_site_rank'] = _target_site_rank(
        comparison, crop_rank
    )
    comparison['full_vs_fimo_rank_gain'] = (
        comparison['fimo_rank'] - comparison['full_rank']
    )
    comparison[f'{crop_arm}_vs_fimo_rank_gain'] = (
        comparison['fimo_rank'] - comparison[crop_rank]
    )
    comparison[f'full_minus_{crop_arm}_ipsae'] = (
        comparison['full_ipsae'] - comparison[crop_ipsae]
    )
    comparison['best_structure_ipsae'] = comparison[
        ['full_ipsae', crop_ipsae]
    ].max(axis=1)
    comparison['best_structure_arm'] = np.select(
        [
            comparison['full_ipsae'].gt(comparison[crop_ipsae]),
            comparison['full_ipsae'].lt(comparison[crop_ipsae]),
        ],
        ['full', crop_arm],
        default='tie',
    )
    comparison['TF'] = comparison['tf_name']
    comparison['TG'] = comparison['sequence_name']

    results_dir = run_dir / 'results'
    fimo_path = results_dir / 'fimo_only_scored.csv'
    full_path = results_dir / f'fimo_plus_{model}_scored.csv'
    crop_path = results_dir / f'fimo_plus_{crop_arm}_{model}_scored.csv'
    comparison_path = results_dir / 'three_arm_scored_comparison.csv'
    metadata_path = results_dir / 'three_arm_scoring_summary.json'

    fimo_columns = [
        'fimo_rank', 'fimo_target_site_rank', 'fimo_neglog10_pvalue',
    ] + base_columns
    fimo_columns = [
        {'score': 'fimo_pwm_score', 'pvalue': 'fimo_pvalue', 'qvalue': 'fimo_qvalue'}.get(
            column, column
        )
        for column in fimo_columns
    ]
    fimo_table = comparison.loc[:, list(dict.fromkeys(fimo_columns))].sort_values(
        'fimo_rank', kind='mergesort'
    )
    fimo_table.insert(1, 'ranking_score', fimo_table['fimo_neglog10_pvalue'])

    common = [column for column in fimo_table.columns if column != 'ranking_score']
    full_metrics = [f'full_{metric}' for metric in _STRUCTURE_METRICS]
    full_columns = [
        'full_rank', 'full_target_site_rank', 'full_ipsae', 'fimo_rank',
        'fimo_neglog10_pvalue',
    ] + common + full_metrics + ['full_vs_fimo_rank_gain']
    full_columns = list(dict.fromkeys(
        column for column in full_columns if column in comparison.columns
    ))
    full_table = comparison.loc[:, full_columns].sort_values(
        'full_rank', kind='mergesort'
    )
    full_table.insert(1, 'ranking_score', full_table['full_ipsae'])

    crop_metrics = [f'{crop_arm}_{metric}' for metric in _STRUCTURE_METRICS]
    crop_columns = [
        crop_rank, f'{crop_arm}_target_site_rank', crop_ipsae, 'fimo_rank',
        'fimo_neglog10_pvalue',
    ] + common + crop_metrics + [f'{crop_arm}_vs_fimo_rank_gain']
    crop_columns = list(dict.fromkeys(
        column for column in crop_columns if column in comparison.columns
    ))
    crop_table = comparison.loc[:, crop_columns].sort_values(
        crop_rank, kind='mergesort'
    )
    crop_table.insert(1, 'ranking_score', crop_table[crop_ipsae])

    fimo_table.to_csv(fimo_path, index=False)
    full_table.to_csv(full_path, index=False)
    crop_table.to_csv(crop_path, index=False)
    comparison.sort_values('fimo_rank', kind='mergesort').to_csv(
        comparison_path, index=False
    )

    summary = {
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'candidate_level': 'promoter-level; identical dsDNA on different promoters retained',
        'candidate_gate': 'the same retained FIMO candidate universe for all three arms',
        'n_candidates': len(comparison),
        'n_model_inputs': int(comparison['model_id'].nunique()),
        'n_target_promoters': int(comparison['sequence_name'].nunique()),
        'n_full_scores': int(comparison['full_ipsae'].notna().sum()),
        f'n_{crop_arm}_scores': int(comparison[crop_ipsae].notna().sum()),
        'ranking': {
            'fimo_only': (
                'pvalue ascending, PWM score descending, then promoter/coordinate/strand/motif'
            ),
            f'fimo_plus_{model}': (
                'full ipSAE descending; FIMO statistics then promoter/site fields are tie-breakers'
            ),
            f'fimo_plus_{crop_arm}_{model}': (
                f'{crop_arm} ipSAE descending; FIMO statistics then promoter/site fields are tie-breakers'
            ),
        },
        'weighted_composite_score': None,
        'rank_gain_definition': 'FIMO rank minus structure-arm rank; positive means promoted by shape readout',
        'rank_spearman': {
            'fimo_vs_full': float(comparison['fimo_rank'].corr(
                comparison['full_rank'], method='spearman'
            )),
            f'fimo_vs_{crop_arm}': float(comparison['fimo_rank'].corr(
                comparison[crop_rank], method='spearman'
            )),
            f'full_vs_{crop_arm}': float(comparison['full_rank'].corr(
                comparison[crop_rank], method='spearman'
            )),
        },
        'files': {
            'fimo_only': str(fimo_path),
            f'fimo_plus_{model}': str(full_path),
            f'fimo_plus_{crop_arm}_{model}': str(crop_path),
            'comparison': str(comparison_path),
        },
    }
    metadata_path.write_text(json.dumps(summary, indent=2) + '\n')
    return ThreeArmRankingSummary(
        fimo=fimo_path,
        full=full_path,
        crop=crop_path,
        comparison=comparison_path,
        metadata=metadata_path,
        n_candidates=len(comparison),
        n_model_inputs=int(comparison['model_id'].nunique()),
        n_target_promoters=int(comparison['sequence_name'].nunique()),
    )
