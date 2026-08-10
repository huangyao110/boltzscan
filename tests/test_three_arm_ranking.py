import json

import pandas as pd

from boltzscan.ranking import write_three_arm_score_tables


def test_three_arm_rankings_keep_promoter_edges_and_use_arm_specific_ipsae(tmp_path):
    run = tmp_path / 'run'
    (run / 'candidates').mkdir(parents=True)
    (run / 'results').mkdir()
    candidates = pd.DataFrame([
        {
            'candidate_id': 'c1', 'model_id': 'm1', 'tf_name': 'TF1',
            'sequence_name': 'gene1', 'motif_id': 'motif', 'start': 1, 'stop': 4,
            'strand': '+', 'score': 12.0, 'pvalue': 1e-8, 'qvalue': 0.1,
            'extracted_motif_seq': 'GCAT', 'canonical_dsdna': 'ATGC',
            'is_model_representative': True,
        },
        {
            'candidate_id': 'c2', 'model_id': 'm2', 'tf_name': 'TF1',
            'sequence_name': 'gene2', 'motif_id': 'motif', 'start': 5, 'stop': 8,
            'strand': '+', 'score': 10.0, 'pvalue': 1e-6, 'qvalue': 0.2,
            'extracted_motif_seq': 'AAAA', 'canonical_dsdna': 'AAAA',
            'is_model_representative': True,
        },
        {
            # Same physical duplex/model as c1 on another promoter: retain it.
            'candidate_id': 'c3', 'model_id': 'm1', 'tf_name': 'TF1',
            'sequence_name': 'gene3', 'motif_id': 'motif', 'start': 9, 'stop': 12,
            'strand': '-', 'score': 11.0, 'pvalue': 1e-7, 'qvalue': 0.15,
            'extracted_motif_seq': 'ATGC', 'canonical_dsdna': 'ATGC',
            'is_model_representative': False,
        },
    ])
    candidates.to_csv(
        run / 'candidates' / 'promoter_candidates.csv.gz',
        index=False,
        compression='gzip',
    )
    base = candidates.loc[:, [
        'candidate_id', 'model_id', 'tf_name', 'sequence_name', 'motif_id',
        'start', 'stop', 'strand', 'score', 'pvalue', 'extracted_motif_seq',
        'canonical_dsdna',
    ]]
    full = base.assign(ipsae=[0.2, 0.8, 0.2])
    crop = base.assign(ipsae=[0.9, 0.3, 0.9])
    full.to_csv(
        run / 'results' / 'full_promoter_level_ipsae.csv.gz',
        index=False,
        compression='gzip',
    )
    crop.to_csv(
        run / 'results' / 'crop20_promoter_level_ipsae.csv.gz',
        index=False,
        compression='gzip',
    )

    summary = write_three_arm_score_tables(run, crop=20, model='boltz2')
    comparison = pd.read_csv(summary.comparison).set_index('candidate_id')
    assert len(comparison) == 3
    assert comparison.model_id.nunique() == 2
    assert comparison.loc['c1', 'fimo_rank'] == 1
    assert comparison.loc['c2', 'full_rank'] == 1
    assert comparison.loc['c1', 'crop20_rank'] == 1
    assert comparison.loc['c3', 'crop20_rank'] == 2
    assert comparison.loc['c1', 'best_structure_arm'] == 'crop20'
    metadata = json.loads(summary.metadata.read_text())
    assert metadata['weighted_composite_score'] is None
    assert metadata['n_candidates'] == 3
    assert summary.full.name == 'fimo_plus_boltz2_scored.csv'
    assert summary.crop.name == 'fimo_plus_crop20_boltz2_scored.csv'
