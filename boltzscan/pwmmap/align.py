"""blastp-based DBD %ID: fast all-vs-all with domain-level, coverage-gated identity."""
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from boltzscan.toolchain import resolve_executable

OUTFMT = "6 qseqid sseqid pident length nident qlen slen"


@dataclass
class Hit:
    query_tf: str
    query_pfam: str
    ref_dbd_id: str
    ref_pfam: str
    pct_id: float


def resolve_blast_bins(blastp=None, makeblastdb=None):
    bp = resolve_executable("blastp", explicit=blastp)
    mk = resolve_executable("makeblastdb", explicit=makeblastdb)
    if not bp or not mk:
        raise FileNotFoundError(
            "BLAST+ executables not found; run `boltzscan doctor --fix`"
        )
    return bp, mk


def _pfam_of(seq_id):
    parts = seq_id.split("__")
    return parts[1] if len(parts) >= 2 else ""


def _tf_of(seq_id):
    return seq_id.split("__")[0]


def make_blast_db(ref_dbd_fasta, db_path, makeblastdb):
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([makeblastdb, "-in", str(ref_dbd_fasta), "-dbtype", "prot",
                    "-out", str(db_path)], check=True)
    return db_path


def parse_blast_rows(rows, min_cov):
    best = {}
    for r in rows:
        q, s = r[0], r[1]
        length, nident, qlen, slen = int(r[3]), int(r[4]), int(r[5]), int(r[6])
        qpf, spf = _pfam_of(q), _pfam_of(s)
        if qpf != spf or not qpf:
            continue
        if length < min_cov * qlen or length < min_cov * slen:
            continue
        pid = nident / length if length else 0.0
        key = (q, s)
        if key not in best or pid > best[key].pct_id:
            best[key] = Hit(_tf_of(q), qpf, s, spf, pid)
    return list(best.values())


def blast_dbd_pct_id(query_fasta, db_path, blastp, cpu=8, min_cov=0.8):
    with tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=False) as tmp:
        out = tmp.name
    try:
        subprocess.run([blastp, "-query", str(query_fasta), "-db", str(db_path),
                        "-outfmt", OUTFMT, "-num_threads", str(cpu),
                        "-max_target_seqs", "2000", "-out", out], check=True)
        rows = [ln.split("\t") for ln in Path(out).read_text().splitlines() if ln.strip()]
    finally:
        Path(out).unlink(missing_ok=True)
    return parse_blast_rows(rows, min_cov)
