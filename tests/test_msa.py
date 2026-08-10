import io
import tarfile

import pytest

from boltzscan import msa
from boltzscan.cli import _build_parser
from boltzscan.msa import crop_a3m_file, read_a3m


def _result_archive(a3m_text):
    buffer = io.BytesIO()
    payload = a3m_text.encode()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        member = tarfile.TarInfo("0.a3m")
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))
    return buffer.getvalue()


class _Response:
    def __init__(self, *, json_data=None, content=b""):
        self._json_data = json_data
        self._content = content

    def raise_for_status(self):
        return None

    def json(self):
        return self._json_data

    def iter_content(self, chunk_size):
        del chunk_size
        yield self._content


class _Session:
    archive = b""
    submissions = []

    def __init__(self):
        self.headers = {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def request(self, method, url, **kwargs):
        self.submissions.append((method, url, kwargs))
        return _Response(json_data={"id": "job-1", "status": "COMPLETE"})

    def get(self, url, **kwargs):
        return _Response(content=self.archive)


def test_run_msa_uses_server_and_keeps_only_final_a3m(tmp_path, monkeypatch):
    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(">TF.1 description\nACDEFG\n")
    _Session.archive = _result_archive(">query_0\nACDEFG\n>hit\nAC-EFG\n")
    _Session.submissions = []
    monkeypatch.setattr(msa.requests, "Session", _Session)

    summary = msa.run_msa(fasta, tmp_path / "msa", server_url="https://msa.test/api")

    output = tmp_path / "msa" / "TF_1" / "0.a3m"
    assert output.read_text().startswith(">query_0\nACDEFG\n")
    assert summary.requested == 1
    assert summary.completed == 1
    assert summary.skipped == 0
    assert list((tmp_path / "msa").rglob("*.tar.gz")) == []
    assert _Session.submissions[0][0:2] == (
        "POST", "https://msa.test/api/ticket/msa"
    )
    assert _Session.submissions[0][2]["data"]["q"] == ">query_0\nACDEFG\n"


def test_run_msa_reuses_matching_existing_a3m(tmp_path, monkeypatch):
    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(">TF-1\nACDEFG\n")
    output = tmp_path / "msa" / "TF_1" / "0.a3m"
    output.parent.mkdir(parents=True)
    output.write_text(">query_0\nACDEFG\n")

    def unexpected_session():
        raise AssertionError("server should not be called")

    monkeypatch.setattr(msa.requests, "Session", unexpected_session)
    summary = msa.run_msa(fasta, tmp_path / "msa")

    assert summary.completed == 0
    assert summary.skipped == 1


def test_run_msa_rejects_existing_msa_for_different_protein(tmp_path):
    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(">TF-1\nACDEFG\n")
    output = tmp_path / "msa" / "TF_1" / "0.a3m"
    output.parent.mkdir(parents=True)
    output.write_text(">query_0\nAAAAAA\n")

    with pytest.raises(RuntimeError, match="Existing MSA query does not match TF-1"):
        msa.run_msa(fasta, tmp_path / "msa")


def test_run_msa_rejects_illegal_protein_characters(tmp_path):
    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(">TF1\nACD*EFG\n")

    with pytest.raises(ValueError, match=r"Illegal protein character.*\*"):
        msa.run_msa(fasta, tmp_path / "msa")


def test_msa_cli_uses_public_flags():
    args = _build_parser().parse_args(
        ["msa", "--fasta", "proteins.fa", "--output", "msa", "--jobs", "2"]
    )
    assert args.fasta_file == "proteins.fa"
    assert args.output_dir == "msa"
    assert args.jobs == 2


def test_crop_a3m_uses_query_residues_not_lowercase_insertions(tmp_path):
    source = tmp_path / '0.a3m'
    source.write_text(
        '>query\nAB-CD-EF\n'
        '>hit\nAqB-CrD-sEtF\n'
        '>all_gap\n--------\n'
        '>last_record\nXB-YD-ZF\n'
    )

    output = crop_a3m_file(source, tmp_path / 'crop_0.a3m', interval=(1, 5))

    records = read_a3m(output)
    assert [record.id for record in records] == ['query', 'hit', 'last_record']
    assert [str(record.seq) for record in records] == ['B-CD-E', 'B-CrD-sE', 'B-YD-Z']
    assert msa.read_a3m_query(output) == 'BCDE'
