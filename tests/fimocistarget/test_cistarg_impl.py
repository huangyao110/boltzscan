import pandas as pd

from boltzscan.fimocistarget.cistarg_impl import filter_fimo_by_seq_and_overlap


def _hits(rows):
    return pd.DataFrame(rows)


def test_fimo_coordinates_are_1_based_inclusive():
    result = filter_fimo_by_seq_and_overlap(_hits([{
        'sequence_name': 'p1', 'tf_name': 'tf', 'promoter_seq': 'AACCGGTTAACC',
        'start': 3, 'stop': 5, 'strand': '+', 'score': 1.0,
        'pvalue': 1e-6, 'motif_id': 'm1',
    }]), extend_range=2)

    hit = result.iloc[0]
    assert hit.extracted_motif_seq == 'AACCGGT'
    assert (hit.real_start, hit.real_stop, hit.length) == (0, 7, 7)


def test_overlap_uses_core_and_prefers_lower_pvalue():
    result = filter_fimo_by_seq_and_overlap(_hits([
        {
            'sequence_name': 'p1', 'tf_name': 'tf',
            'promoter_seq': 'AACCTGAGTTCGATCGGATACCTGAGCACTTGA',
            'start': 12, 'stop': 14, 'strand': '+', 'score': 20.0,
            'pvalue': 1e-4, 'motif_id': 'higher_score',
        },
        {
            'sequence_name': 'p1', 'tf_name': 'tf',
            'promoter_seq': 'AACCTGAGTTCGATCGGATACCTGAGCACTTGA',
            'start': 13, 'stop': 15, 'strand': '-', 'score': 10.0,
            'pvalue': 1e-6, 'motif_id': 'lower_pvalue',
        },
        {
            'sequence_name': 'p1', 'tf_name': 'tf',
            'promoter_seq': 'AACCTGAGTTCGATCGGATACCTGAGCACTTGA',
            'start': 16, 'stop': 17, 'strand': '+', 'score': 1.0,
            'pvalue': 1e-5, 'motif_id': 'adjacent_core',
        },
    ]), extend_range=5)

    assert set(result.motif_id) == {'lower_pvalue', 'adjacent_core'}


def test_equal_hits_have_a_deterministic_strand_choice():
    rows = [
        {
            'sequence_name': 'p1', 'tf_name': 'tf', 'promoter_seq': 'AACCGGTTAACC',
            'start': 6, 'stop': 4, 'strand': '-', 'score': 1.0,
            'pvalue': 1e-6, 'motif_id': 'm1',
        },
        {
            'sequence_name': 'p1', 'tf_name': 'tf', 'promoter_seq': 'AACCGGTTAACC',
            'start': 4, 'stop': 6, 'strand': '+', 'score': 1.0,
            'pvalue': 1e-6, 'motif_id': 'm1',
        },
    ]

    forward = filter_fimo_by_seq_and_overlap(_hits(rows), extend_range=0)
    reversed_rows = filter_fimo_by_seq_and_overlap(_hits(list(reversed(rows))), extend_range=0)
    assert forward.iloc[0].strand == reversed_rows.iloc[0].strand == '+'
