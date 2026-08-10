import gzip

import pytest

from boltzscan.pwmmap.pfam import (
    extract_runtime_pfam,
    packaged_runtime_pfam_paths,
    required_pfam_accessions,
    resolve_runtime_pfam,
    runtime_pfam_paths,
    validate_runtime_pfam,
)


def _full_pfam_text(accessions):
    return ''.join(
        f'HMMER3/f\nNAME  test_{accession}\nACC   {accession}.42\nDATE  2026-01-01\n//\n'
        for accession in accessions
    )


def test_extracts_exact_runtime_profile_set_and_resolves_from_refs(tmp_path):
    required = sorted(required_pfam_accessions())
    source = tmp_path / 'Pfam-A.hmm.gz'
    with gzip.open(source, 'wt') as handle:
        handle.write(_full_pfam_text([*required, 'PF99999']))

    summary = extract_runtime_pfam(source, tmp_path / 'refs' / 'pfam')

    assert summary.n_profiles == len(required) == 70
    assert summary.size_bytes < source.stat().st_size * 20
    assert resolve_runtime_pfam(tmp_path / 'refs') == summary.hmm.resolve()
    assert validate_runtime_pfam(*runtime_pfam_paths(tmp_path / 'refs')) == summary


def test_extraction_rejects_a_missing_required_profile(tmp_path):
    required = sorted(required_pfam_accessions())
    source = tmp_path / 'incomplete.hmm'
    source.write_text(_full_pfam_text(required[:-1]))

    with pytest.raises(ValueError, match='missing 1 required profiles'):
        extract_runtime_pfam(source, tmp_path / 'pfam')


def test_packaged_runtime_library_is_complete():
    summary = validate_runtime_pfam(*packaged_runtime_pfam_paths())

    assert summary.n_profiles == 70
    assert summary.size_bytes < 4 * 1024 * 1024
