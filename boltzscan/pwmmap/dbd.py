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
from boltzscan.pwmmap.thresholds import family_for_pfam

DEFAULT_PFAM = "data/pfam/Pfam-A.hmm"


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


def _run_hmmsearch(proteins, pfam, domtbl, cpu):
    domtbl.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["hmmsearch", "--cut_ga", "--cpu", str(cpu),
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
    with open(path, "w") as fh:
        for i, r in enumerate(records):
            fh.write(f">{r.tf_id}__{r.pfam_acc}__{i}\n{r.seq}\n")
