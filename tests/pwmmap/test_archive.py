import json
from io import BytesIO

import pytest

from boltzscan.pwmmap import archive as archive_module
from boltzscan.pwmmap.archive import (
    DEFAULT_REFERENCE_RELEASE_URL,
    _download_url,
    _preferred_release_source,
    install_reference_store,
    pack_reference_store,
    validate_runtime_store,
)
from boltzscan.pwmmap.pfam import extract_runtime_pfam, required_pfam_accessions


def _runtime_store(tmp_path):
    refs = tmp_path / '_refs'
    (refs / 'motif_store' / 'txt').mkdir(parents=True)
    (refs / 'motif_store' / 'meme').mkdir(parents=True)
    (refs / 'build_manifest.json').write_text(json.dumps({'n_refs': 1}))
    (refs / 'ref_dbd.fasta').write_text('>ref__PF00001__0\nAAAA\n')
    (refs / 'ref_index.tsv').write_text(
        'ref_id\tsource\tspecies\tfamily\tpfam_acc\tdbd_seq_id\tmotif_ids\n'
        'ref\ttest\tSpecies\tFamily\tPF00001\tref__PF00001__0\tM1\n'
    )
    (refs / 'ref_proteins.fasta').write_text('>ref\nAAAA\n')
    (refs / 'motif_clusters.tsv').write_text(
        'motif_id\tfamily\trepresentative_id\nM1\tFamily\tM1\n'
    )
    (refs / 'motif_store' / 'txt' / 'M1.txt').write_text('Pos\tA\tC\tG\tT\n1\t1\t0\t0\t0\n')
    (refs / 'motif_store' / 'meme' / 'M1.meme').write_text('MEME version 5\n\nMOTIF M1\n')
    full_pfam = tmp_path / 'Pfam-A.hmm'
    full_pfam.write_text(''.join(
        f'HMMER3/f\nNAME  test_{accession}\nACC   {accession}.1\n//\n'
        for accession in sorted(required_pfam_accessions())
    ))
    extract_runtime_pfam(full_pfam, refs / 'pfam')
    return refs


def test_pack_and_install_runtime_reference_store(tmp_path):
    source = _runtime_store(tmp_path)
    archive = tmp_path / 'release' / 'boltzscan_pwm_refs.tar.gz'

    packed = pack_reference_store(source, archive)

    assert packed.archive == archive
    assert packed.checksum_file.read_text().startswith(packed.sha256)
    assert packed.n_dbd_rows == packed.n_txt_motifs == packed.n_meme_motifs == 1
    assert packed.n_pfam_profiles == len(required_pfam_accessions())
    assert (source / 'release_manifest.json').is_file()

    progress = []
    installed = install_reference_store(
        archive,
        tmp_path / 'installed_refs',
        packed.sha256,
        progress=progress.append,
    )
    assert installed.archive_sha256 == packed.sha256
    assert validate_runtime_store(installed.refs_dir)['n_dbd_rows'] == 1
    assert (installed.refs_dir / 'motif_store' / 'meme' / 'M1.meme').is_file()
    assert (installed.refs_dir / 'pfam' / 'plant_tf_pfam.hmm').is_file()
    assert installed.n_pfam_profiles == len(required_pfam_accessions())
    assert any(message.startswith('Archive copied:') for message in progress)
    assert 'SHA256 verified.' in progress
    assert 'Extracting PWM reference archive...' in progress
    assert 'Reference store validated; installing atomically...' in progress


def test_install_rejects_wrong_checksum_without_touching_destination(tmp_path):
    source = _runtime_store(tmp_path)
    packed = pack_reference_store(source, tmp_path / 'release.tar.gz')
    destination = tmp_path / 'destination'

    with pytest.raises(ValueError, match='SHA256 mismatch'):
        install_reference_store(packed.archive, destination, '0' * 64)

    assert not destination.exists()


def test_default_release_prefers_a_matching_source_checkout_archive(
    tmp_path, monkeypatch,
):
    source = _runtime_store(tmp_path)
    packed = pack_reference_store(source, tmp_path / 'local-release.tar.gz')
    monkeypatch.setattr(
        archive_module,
        'DEFAULT_LOCAL_REFERENCE_RELEASE',
        packed.archive,
    )
    progress = []

    selected = _preferred_release_source(
        DEFAULT_REFERENCE_RELEASE_URL,
        packed.sha256,
        progress=progress.append,
    )

    assert selected == packed.archive
    assert progress == [f'Using verified source-checkout PWM release: {packed.archive}']


def test_google_drive_share_url_is_converted_to_direct_download():
    share_url = 'https://drive.google.com/file/d/FILE_ID/view?usp=drive_link'

    assert _download_url(share_url) == (
        'https://drive.usercontent.google.com/download'
        '?id=FILE_ID&export=download&confirm=t'
    )
    assert _download_url('/tmp/local-release.tar.gz') == '/tmp/local-release.tar.gz'


def test_remote_download_reports_connection_and_progress(tmp_path, monkeypatch):
    payload = b'x' * 2048

    class Response(BytesIO):
        headers = {'Content-Length': str(len(payload))}

    monkeypatch.setattr(
        archive_module,
        'urlopen',
        lambda request, timeout: Response(payload),
    )
    progress = []
    output = tmp_path / 'download.tar.gz'

    archive_module._download(
        'https://example.org/release.tar.gz',
        output,
        progress=progress.append,
    )

    assert output.read_bytes() == payload
    assert progress[0].startswith('Connecting to example.org')
    assert any(message.startswith('Download:  10%') for message in progress)
    assert any(message.startswith('Download: 100%') for message in progress)
    assert progress[-1] == 'Download complete: 2.0 KiB'


def test_remote_download_times_out_with_an_offline_install_hint(tmp_path, monkeypatch):
    calls = []

    def timeout(request, *, timeout):
        calls.append(timeout)
        raise TimeoutError('timed out')

    monkeypatch.setattr(archive_module, 'urlopen', timeout)

    with pytest.raises(RuntimeError, match='install-pwm-refs --url FILE'):
        archive_module._download(
            'https://example.org/release.tar.gz',
            tmp_path / 'download.tar.gz',
        )

    assert calls == [archive_module.DOWNLOAD_TIMEOUT_SECONDS]
