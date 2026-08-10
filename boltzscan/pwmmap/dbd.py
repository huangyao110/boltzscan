"""Extract DNA-binding-domain subsequences from proteins using Pfam-A hits.

Both the reference proteins (Stage A) and the species TFs (Stage B) go through
this so the two sides use identical domain boundaries -> comparable %ID.
"""
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from boltzscan.utils.find_tf import DBD_ACCS
from boltzscan.pwmmap.pfam import packaged_runtime_pfam_paths
from boltzscan.pwmmap.thresholds import family_for_pfam
from boltzscan.toolchain import resolve_executable

DEFAULT_PFAM = str(packaged_runtime_pfam_paths()[0])


@dataclass
class DbdRecord:
    tf_id: str
    pfam_acc: str
    family: str | None
    start: int   # 1-based inclusive (envelope)
    end: int
    seq: str


def _base_acc(a):
    return a.split(".")[0]


def parse_domtbl_dbds(domtbl_path):
    """Return [(tf_id, pfam_acc, env_from, env_to)] for DBD Pfam hits only."""
    out = []
    with open(domtbl_path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            c = ln.split()
            acc = _base_acc(c[4])
            if acc not in DBD_ACCS:
                continue
            out.append((c[0], acc, int(c[19]), int(c[20])))
    return out


def dbd_crop_interval(records, tf_id, protein_length, flank):
    """Return the 1-based inclusive union-of-DBDs crop interval for one TF.

    ``records`` may be the complete output of :func:`parse_domtbl_dbds` or a
    pre-filtered subset.  A crop is meaningful only when a DNA-binding-domain
    hit exists, so absence of a hit is an error rather than a silent full-length
    fallback.  Multiple DBDs are represented by their bounding union.
    """
    if protein_length <= 0:
        raise ValueError(f"Protein length must be positive for {tf_id!r}")
    if flank < 0:
        raise ValueError("Crop flank must be non-negative")
    hits = [record for record in records if record[0] == tf_id]
    if not hits:
        raise ValueError(f"No DNA-binding-domain Pfam hit found for {tf_id!r}")

    dbd_start = min(int(record[2]) for record in hits)
    dbd_stop = max(int(record[3]) for record in hits)
    if dbd_stop < 1 or dbd_start > protein_length:
        raise ValueError(
            f"DBD coordinates {dbd_start}-{dbd_stop} for {tf_id!r} do not overlap "
            f"its protein length ({protein_length})"
        )
    dbd_start = max(1, dbd_start)
    dbd_stop = min(protein_length, dbd_stop)
    return max(1, dbd_start - flank), min(protein_length, dbd_stop + flank)


def _run_hmmsearch(proteins, pfam, domtbl, cpu):
    domtbl.parent.mkdir(parents=True, exist_ok=True)
    hmmsearch = resolve_executable("hmmsearch")
    if hmmsearch is None:
        raise FileNotFoundError(
            "HMMER hmmsearch not found; run `boltzscan doctor --fix`"
        )
    cmd = [hmmsearch, "--cut_ga", "--cpu", str(cpu),
           "--domtblout", str(domtbl), "-o", "/dev/null", str(pfam), str(proteins)]
    print("[dbd] " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def extract_dbds(proteins_fasta, domtbl=None, pfam=DEFAULT_PFAM, cpu=8, work_dir=None):
    proteins_fasta = Path(proteins_fasta)
    if domtbl is None:
        work_dir = Path(work_dir) if work_dir else proteins_fasta.parent
        domtbl = work_dir / "pfam.domtbl"
        if not domtbl.exists():
            _run_hmmsearch(proteins_fasta, Path(pfam), domtbl, cpu)
    rows = parse_domtbl_dbds(domtbl)
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(str(proteins_fasta), "fasta")}
    recs = []
    for tf_id, acc, a, b in rows:
        s = seqs.get(tf_id)
        if not s:
            continue
        recs.append(DbdRecord(tf_id, acc, family_for_pfam(acc), a, b, s[a - 1:b]))
    return recs


def write_dbd_fasta(records, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    written_ids = []
    with open(path, "w") as fh:
        for i, r in enumerate(records):
            seq_id = f"{r.tf_id}__{r.pfam_acc}__{i}"
            fh.write(f">{seq_id}\n{r.seq}\n")
            written_ids.append(seq_id)
    return written_ids
