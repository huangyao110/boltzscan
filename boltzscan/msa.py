"""Minimal client for the Protenix-compatible remote MSA service.

This module intentionally implements only the HTTP workflow needed by
BoltzScan: submit one protein sequence, wait for completion, and retain the
returned ``0.a3m``.  It does not import or execute Protenix.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
import os
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
import time

from Bio import SeqIO
import requests

from boltzscan import __version__


DEFAULT_MSA_SERVER_URL = os.environ.get(
    "BOLTZSCAN_MSA_SERVER_URL", "https://protenix-server.com/api/msa"
)
_PROTEIN_ALPHABET = frozenset("ARNDCEQGHILKMFPSTWYVXJBZOU")
_PENDING_STATUSES = frozenset({"UNKNOWN", "RUNNING", "PENDING"})
_SUBMIT_AGAIN_STATUSES = frozenset({"UNKNOWN", "RATELIMIT"})


@dataclass(frozen=True)
class MsaRunSummary:
    requested: int
    completed: int
    skipped: int
    outputs: tuple[Path, ...]


@dataclass(frozen=True)
class _Query:
    record_id: str
    sequence: str
    output: Path


def _safe_record_id(record_id: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9]", "_", record_id)
    if not safe:
        raise ValueError(f"FASTA record has no usable identifier: {record_id!r}")
    return safe


def _read_queries(fasta_files, output_dir: Path) -> list[_Query]:
    queries = []
    seen = {}
    for fasta_file in fasta_files:
        fasta_path = Path(fasta_file)
        if not fasta_path.is_file():
            raise FileNotFoundError(f"Protein FASTA not found: {fasta_path}")
        records = list(SeqIO.parse(str(fasta_path), "fasta"))
        if not records:
            raise ValueError(f"No FASTA records found: {fasta_path}")
        for record in records:
            record_id = str(record.id)
            safe_id = _safe_record_id(record_id)
            if safe_id in seen:
                raise ValueError(
                    f"FASTA identifiers {seen[safe_id]!r} and {record_id!r} "
                    f"both map to output directory {safe_id!r}"
                )
            seen[safe_id] = record_id
            sequence = str(record.seq).upper()
            invalid = sorted(set(sequence).difference(_PROTEIN_ALPHABET))
            if not sequence:
                raise ValueError(f"Empty protein sequence: {record_id}")
            if invalid:
                raise ValueError(
                    f"Illegal protein character(s) {invalid} in FASTA record {record_id}"
                )
            queries.append(
                _Query(record_id, sequence, output_dir / safe_id / "0.a3m")
            )
    return queries


def read_a3m_query(path: Path) -> str:
    """Return the ungapped query sequence from the first A3M record."""
    sequence = []
    in_first_record = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if in_first_record:
                    break
                in_first_record = True
                continue
            if in_first_record:
                sequence.append(line)
    if not in_first_record or not sequence:
        raise ValueError(f"MSA contains no query sequence: {path}")
    return "".join(
        char for char in "".join(sequence) if char != "-" and not char.islower()
    ).upper()


def read_a3m(path):
    """Read A3M records while preserving lowercase insertion characters."""
    return list(SeqIO.parse(path, "fasta"))


def _a3m_groups(sequence):
    """Group each alignment character with lowercase insertions before it."""
    groups = []
    insertions = ""
    for char in sequence:
        if char.islower():
            insertions += char
        else:
            groups.append((insertions, char))
            insertions = ""
    return groups


def crop_a3m_file(a3m_file, out_path, interval=None):
    """Crop A3M rows to a 0-based, half-open query-residue interval."""
    records = read_a3m(a3m_file)
    if not records:
        raise ValueError(f"A3M contains no sequences: {a3m_file}")

    query_groups = _a3m_groups(str(records[0].seq))
    residue_columns = [
        index for index, (_, char) in enumerate(query_groups) if char != "-"
    ]
    if not residue_columns:
        raise ValueError(f"A3M query contains no residues: {a3m_file}")
    start, stop = interval if interval is not None else (0, len(residue_columns))
    if not 0 <= start < stop <= len(residue_columns):
        raise ValueError(
            f"Invalid A3M crop interval [{start}, {stop}) for query length "
            f"{len(residue_columns)}"
        )
    first_column = residue_columns[start]
    last_column = residue_columns[stop - 1]

    cropped_records = []
    for record in records:
        groups = _a3m_groups(str(record.seq))
        if len(groups) != len(query_groups):
            raise ValueError(
                f"A3M alignment width mismatch for {record.id}: "
                f"{len(groups)} != {len(query_groups)}"
            )
        selected = groups[first_column:last_column + 1]
        crop_seq = selected[0][1] + "".join(
            insertions + char for insertions, char in selected[1:]
        )
        if record is not records[0] and set(
            char for char in crop_seq if not char.islower()
        ) <= {"-"}:
            continue
        cropped_records.append((record.description, crop_seq))

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    staged_path = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            prefix=f".{out_path.name}.",
            dir=out_path.parent,
            delete=False,
            encoding="utf-8",
        ) as staged:
            staged_path = Path(staged.name)
            for description, sequence in cropped_records:
                staged.write(f">{description}\n{sequence}\n")
        staged_path.replace(out_path)
    finally:
        if staged_path is not None:
            staged_path.unlink(missing_ok=True)
    return out_path


def _request_json(session, method, url, **kwargs):
    last_error = None
    for attempt in range(5):
        try:
            response = session.request(method, url, **kwargs)
            response.raise_for_status()
            result = response.json()
            if not isinstance(result, dict):
                raise ValueError("response is not a JSON object")
            return result
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            if attempt < 4:
                time.sleep(min(2 ** attempt, 10))
    raise RuntimeError(f"MSA server request failed: {last_error}") from None


def _submit(session, server_url: str, sequence: str) -> tuple[str, dict]:
    payload = {"q": f">query_0\n{sequence}\n", "mode": "env", "email": ""}
    for attempt in range(6):
        result = _request_json(
            session,
            "POST",
            f"{server_url}/ticket/msa",
            data=payload,
            timeout=(10, 120),
        )
        status = str(result.get("status", "ERROR")).upper()
        if status not in _SUBMIT_AGAIN_STATUSES:
            if status == "MAINTENANCE":
                raise RuntimeError("MSA server is under maintenance; please try again later")
            if status == "ERROR":
                raise RuntimeError(
                    "MSA server rejected the job; check the protein sequence or try later"
                )
            job_id = result.get("id")
            if not job_id:
                raise RuntimeError(f"MSA server returned no job id (status={status})")
            return str(job_id), result
        if attempt < 5:
            time.sleep(10)
    raise RuntimeError(f"MSA server remained unavailable (status={status})")


def _wait_for_completion(session, server_url: str, job_id: str, result: dict):
    deadline = time.monotonic() + 2 * 60 * 60
    status = str(result.get("status", "UNKNOWN")).upper()
    while status in _PENDING_STATUSES:
        if time.monotonic() >= deadline:
            raise RuntimeError(f"MSA server job timed out after 2 hours: {job_id}")
        time.sleep(10)
        result = _request_json(
            session,
            "GET",
            f"{server_url}/ticket/{job_id}",
            timeout=(10, 120),
        )
        status = str(result.get("status", "ERROR")).upper()
    if status == "MAINTENANCE":
        raise RuntimeError("MSA server is under maintenance; please try again later")
    if status != "COMPLETE":
        raise RuntimeError(f"MSA server job failed: {job_id} (status={status})")


def _download_a3m(session, server_url: str, job_id: str, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    archive_path = None
    staged_path = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=".msa_", suffix=".tar.gz", dir=output.parent, delete=False
        ) as archive:
            archive_path = Path(archive.name)
            response = session.get(
                f"{server_url}/result/download/{job_id}",
                stream=True,
                timeout=(10, 300),
            )
            response.raise_for_status()
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    archive.write(chunk)

        with tarfile.open(archive_path, "r:gz") as tar:
            members = [
                member for member in tar.getmembers()
                if member.isfile() and Path(member.name).name == "0.a3m"
            ]
            if len(members) != 1:
                raise RuntimeError(
                    f"MSA result archive contains {len(members)} files named 0.a3m"
                )
            source = tar.extractfile(members[0])
            if source is None:
                raise RuntimeError("Cannot read 0.a3m from MSA result archive")
            with tempfile.NamedTemporaryFile(
                prefix=".0_", suffix=".a3m", dir=output.parent, delete=False
            ) as staged:
                staged_path = Path(staged.name)
                shutil.copyfileobj(source, staged)
        staged_path.replace(output)
    finally:
        if archive_path is not None:
            archive_path.unlink(missing_ok=True)
        if staged_path is not None:
            staged_path.unlink(missing_ok=True)


def _fetch_query(query: _Query, server_url: str) -> tuple[Path, bool]:
    if query.output.is_file():
        if read_a3m_query(query.output) != query.sequence:
            raise ValueError(
                f"Existing MSA query does not match {query.record_id}: {query.output}"
            )
        print(f"MSA {query.record_id}: reuse {query.output}")
        return query.output, True

    headers = {"User-Agent": f"boltzscan/{__version__}"}
    with requests.Session() as session:
        session.headers.update(headers)
        job_id, result = _submit(session, server_url, query.sequence)
        print(f"MSA {query.record_id}: submitted job {job_id}")
        _wait_for_completion(session, server_url, job_id, result)
        _download_a3m(session, server_url, job_id, query.output)

    if read_a3m_query(query.output) != query.sequence:
        query.output.unlink(missing_ok=True)
        raise ValueError(f"Downloaded MSA query does not match {query.record_id}")
    print(f"MSA {query.record_id}: wrote {query.output}")
    return query.output, False


def run_msa_files(fasta_files, output_dir, jobs=1, server_url=DEFAULT_MSA_SERVER_URL):
    """Generate one ``<record_id>/0.a3m`` per FASTA record."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs = int(jobs)
    if jobs < 1:
        raise ValueError("MSA jobs must be at least 1")
    server_url = str(server_url).rstrip("/")
    if not server_url.startswith(("http://", "https://")):
        raise ValueError(f"Invalid MSA server URL: {server_url}")

    queries = _read_queries(fasta_files, output_dir)
    results = []
    errors = []
    with ThreadPoolExecutor(max_workers=min(jobs, len(queries))) as executor:
        futures = {
            executor.submit(_fetch_query, query, server_url): query
            for query in queries
        }
        for future in as_completed(futures):
            query = futures[future]
            try:
                results.append(future.result())
            except Exception as exc:
                errors.append(f"{query.record_id}: {exc}")
    if errors:
        details = "; ".join(errors[:5])
        if len(errors) > 5:
            details += f"; ... {len(errors) - 5} more"
        raise RuntimeError(f"MSA failed for {len(errors)}/{len(queries)} proteins: {details}")

    outputs = tuple(sorted(path for path, _ in results))
    skipped = sum(was_skipped for _, was_skipped in results)
    return MsaRunSummary(
        requested=len(queries),
        completed=len(results) - skipped,
        skipped=skipped,
        outputs=outputs,
    )


def run_msa(fasta_file, output_dir, jobs=1, server_url=DEFAULT_MSA_SERVER_URL):
    return run_msa_files([fasta_file], output_dir, jobs=jobs, server_url=server_url)
