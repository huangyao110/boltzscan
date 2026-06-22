# Cross-species PWM mapping via DBD homology — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reusable two-stage pipeline that transfers cisBP + JASPAR PWMs onto any species' TF set via DBD sequence identity, producing the `data/pwms/OB_rose_pwms` layout (`tf2pwms.json` + `txt/` + `meme/`); first target `data/pwms/cm_pwms/` for 野菊.

**Architecture:** Stage A (`build-pwm-refs`) downloads cisBP + JASPAR once, resolves reference protein sequences via UniProt, extracts reference DBDs with our Pfam-A HMMs, and writes a cached reference store under `data/pwms/_refs/`. Stage B (`map-pwm`) extracts a species' DBDs, blastp's them against the cached reference DBDs, keeps hits clearing a family-specific %ID threshold (same Pfam only), and emits the per-species PWM database. Swapping species re-runs Stage B only.

**Tech Stack:** Python 3.11, Biopython, pandas, numpy, requests, HMMER (`hmmsearch`), NCBI BLAST+ (`blastp`/`makeblastdb` 2.17 from `…/envs/swt/bin`), pytest.

## Global Constraints

- Python target 3.11; do not use 3.13-only syntax (vendored boltz pins `<3.13`).
- New code lives in `boltzscan/pwmmap/`; CLI logic stays thin in `boltzscan/cli.py` (subparser + one-line handler delegating to the module), per repo convention.
- Default Pfam DB: `data/pfam/Pfam-A.hmm` (already hmmpress-ed). Reuse `boltzscan.utils.find_tf.DBD_ACCS` for the plant-TF DBD Pfam set; do not duplicate it.
- BLAST binaries: resolve `blastp`/`makeblastdb` from PATH, else `/home/zlab/miniconda3/envs/swt/bin/`; overridable via `--blastp` / `--makeblastdb`.
- cisBP entire-dataset URLs (verified, build 3_10):
  - `https://cisbp.ccbr.utoronto.ca/data/3_10/DataFiles/Bulk_downloads/EntireDataset/TF_Information.zip`
  - `https://cisbp.ccbr.utoronto.ca/data/3_10/DataFiles/Bulk_downloads/EntireDataset/PWMs.zip`
- JASPAR REST: `https://jaspar.elixir.no/api/v1/matrix/?collection=CORE&page_size=500` (paginate via `next`); detail per `matrix_id` gives `pfm`, `uniprot_ids`, `family`.
- UniProt: sequences `https://rest.uniprot.org/uniprotkb/<acc>.fasta`; ID-mapping API `https://rest.uniprot.org/idmapping/`.
- Commits: the user has a standing "no unsolicited commits" preference. Each task below ends with a commit step (TDD hygiene); when executing, **stage and ask** before committing, or batch commits at the user's request.
- Tests run with `python -m pytest`. Network-touching tests are marked `@pytest.mark.network` and excluded from the default unit run (`-m "not network"`).

---

## File Structure

```
boltzscan/pwmmap/
  __init__.py
  models.py          # RefTf dataclass (shared by sources + refs)
  thresholds.py      # pfam_acc->family, family->%ID cutoff
  dbd.py             # DBD extraction from proteins+domtbl (or run hmmsearch)
  pwmio.py           # JASPAR pfm -> cisBP-style txt + meme; cisBP txt passthrough
  sources/
    __init__.py
    cisbp.py         # download entire dataset, parse high-quality refs
    jaspar.py        # REST download, pfm conversion
    uniprot.py       # DBID/uniprot_id -> protein sequence
  refs.py            # Stage A driver -> data/pwms/_refs/
  align.py           # makeblastdb + blastp + coverage-gated DBD %ID
  mapper.py          # Stage B driver -> data/pwms/<sp>_pwms/
boltzscan/cli.py     # +build-pwm-refs, +map-pwm
tests/pwmmap/
  __init__.py
  conftest.py        # tiny fixtures (fake TF_Information, jaspar page, blast rows)
  test_thresholds.py
  test_dbd.py
  test_pwmio.py
  test_cisbp.py
  test_jaspar.py
  test_uniprot.py
  test_refs.py
  test_align.py
  test_mapper.py
  test_cli.py
```

---

### Task 1: Package scaffold + thresholds table

**Files:**
- Create: `boltzscan/pwmmap/__init__.py` (empty), `boltzscan/pwmmap/sources/__init__.py` (empty)
- Create: `boltzscan/pwmmap/thresholds.py`
- Create: `tests/pwmmap/__init__.py` (empty), `tests/pwmmap/test_thresholds.py`

**Interfaces:**
- Produces: `family_for_pfam(pfam_acc: str) -> str | None`; `cutoff_for(pfam_acc: str, mode: str = "family", global_thr: float = 0.70) -> float`; constants `DEFAULT_CUTOFF = 0.70`, dicts `PFAM_FAMILY`, `FAMILY_CUTOFF`.

- [ ] **Step 1: Write the failing test**

```python
# tests/pwmmap/test_thresholds.py
from boltzscan.pwmmap import thresholds as T

def test_pfam_maps_to_family():
    assert T.family_for_pfam("PF00010") == "bHLH"
    assert T.family_for_pfam("PF00046") == "Homeobox"
    assert T.family_for_pfam("PF99999") is None

def test_family_cutoff_used_in_family_mode():
    # bHLH has a calibrated cutoff distinct from the default
    assert T.cutoff_for("PF00010", mode="family") == T.FAMILY_CUTOFF["bHLH"]

def test_unknown_family_falls_back_to_default():
    assert T.cutoff_for("PF99999", mode="family") == T.DEFAULT_CUTOFF

def test_global_mode_ignores_family():
    assert T.cutoff_for("PF00010", mode="global", global_thr=0.55) == 0.55
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/pwmmap/test_thresholds.py -q`
Expected: FAIL (module `boltzscan.pwmmap.thresholds` not found).

- [ ] **Step 3: Write minimal implementation**

```python
# boltzscan/pwmmap/thresholds.py
"""Pfam DBD accession -> TF family, and family -> DBD %ID transfer cutoff.

Cutoffs follow the cisBP / Weirauch et al. 2014 family-specific DBD identity
thresholds for motif inference (fraction identity over the aligned domain).
Families without a calibrated value use DEFAULT_CUTOFF.
"""

# Canonical Pfam-acc -> family. Seeded from boltzscan.utils.find_tf.DBD_ACCS so
# both modules agree on which Pfam is which family.
PFAM_FAMILY = {
    "PF00046": "Homeobox", "PF05920": "Homeobox", "PF00249": "MYB",
    "PF13837": "MYB", "PF00010": "bHLH", "PF00847": "AP2", "PF02365": "NAC",
    "PF02362": "B3", "PF00096": "C2H2", "PF13894": "C2H2", "PF13912": "C2H2",
    "PF13465": "C2H2", "PF03106": "WRKY", "PF03514": "GRAS", "PF03101": "FAR1",
    "PF00642": "C3H", "PF00170": "bZIP", "PF07716": "bZIP", "PF00319": "MADS",
    "PF03195": "LBD", "PF02701": "Dof", "PF03634": "TCP", "PF00447": "HSF",
    "PF03110": "SBP", "PF00320": "GATA", "PF02319": "E2F", "PF00808": "NF-Y",
    "PF02045": "NF-Y", "PF08879": "GRF", "PF04770": "ZF-HD", "PF02042": "RWP-RK",
    "PF00643": "BBX", "PF06203": "CCT", "PF03859": "CAMTA", "PF04690": "YABBY",
    "PF04504": "GeBP", "PF03638": "CPP", "PF04873": "EIL", "PF06943": "LSD",
    "PF08536": "Whirly", "PF01422": "NF-X1", "PF04689": "S1Fa", "PF01698": "LFY",
    "PF17538": "LFY", "PF05142": "SRS", "PF05687": "BES1", "PF06217": "BBR-BPC",
    "PF28235": "VOZ",
}

DEFAULT_CUTOFF = 0.70

# Family-specific DBD %ID cutoffs (fraction). Values from cisBP/Weirauch 2014;
# families not listed inherit DEFAULT_CUTOFF.
FAMILY_CUTOFF = {
    "Homeobox": 0.70, "bHLH": 0.69, "bZIP": 0.70, "MYB": 0.62, "AP2": 0.59,
    "NAC": 0.62, "B3": 0.60, "C2H2": 0.78, "WRKY": 0.66, "Dof": 0.70,
    "GATA": 0.70, "MADS": 0.70, "HSF": 0.67, "NF-Y": 0.70, "E2F": 0.70,
    "SBP": 0.70, "TCP": 0.66, "GRAS": 0.60, "ZF-HD": 0.70, "LBD": 0.62,
}

def family_for_pfam(pfam_acc):
    return PFAM_FAMILY.get(pfam_acc)

def cutoff_for(pfam_acc, mode="family", global_thr=0.70):
    if mode == "global":
        return global_thr
    fam = family_for_pfam(pfam_acc)
    return FAMILY_CUTOFF.get(fam, DEFAULT_CUTOFF)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/pwmmap/test_thresholds.py -q`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add boltzscan/pwmmap/__init__.py boltzscan/pwmmap/sources/__init__.py \
        boltzscan/pwmmap/thresholds.py tests/pwmmap/__init__.py tests/pwmmap/test_thresholds.py
git commit -m "feat(pwmmap): package scaffold + DBD family/cutoff table"
```

---

### Task 2: DBD extraction from proteins + domtbl

**Files:**
- Create: `boltzscan/pwmmap/dbd.py`
- Create: `tests/pwmmap/test_dbd.py`, `tests/pwmmap/conftest.py`

**Interfaces:**
- Consumes: `boltzscan.utils.find_tf.DBD_ACCS`, `thresholds.family_for_pfam`.
- Produces: dataclass `DbdRecord(tf_id: str, pfam_acc: str, family: str | None, start: int, end: int, seq: str)`; `parse_domtbl_dbds(domtbl_path) -> list[tuple[str,str,int,int]]` returning `(tf_id, pfam_acc, env_from, env_to)` for rows whose Pfam acc is in `DBD_ACCS`; `extract_dbds(proteins_fasta, domtbl, pfam="data/pfam/Pfam-A.hmm", cpu=8, work_dir=None) -> list[DbdRecord]` (runs `hmmsearch --cut_ga` into `work_dir/pfam.domtbl` if `domtbl` is None); `write_dbd_fasta(records, path)` writing IDs `"{tf_id}__{pfam_acc}__{i}"`.

- [ ] **Step 1: Write the failing test**

```python
# tests/pwmmap/conftest.py
import pytest

@pytest.fixture
def tiny_domtbl(tmp_path):
    # domtblout columns: target(0) ... query_name(3) query_acc(4) ... dom_score(13)
    # ... env_from(19) env_to(20). PF00010 is a DBD (bHLH); PF99999 is not.
    rows = [
        "geneA - 120 HLH PF00010.30 54 1e-20 70 0 1 1 1e-22 1e-18 65 0 1 54 11 64 10 65 0.95 -",
        "geneB - 90 SomethingElse PF99999.1 40 1 5 0 1 1 1 1 4 0 1 40 5 44 4 45 0.9 -",
    ]
    p = tmp_path / "pfam.domtbl"
    p.write_text("# comment\n" + "\n".join(" ".join(r.split()) for r in rows) + "\n")
    return p

@pytest.fixture
def tiny_proteins(tmp_path):
    p = tmp_path / "prot.fasta"
    p.write_text(">geneA\n" + "M"*9 + "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRS" + "G"*60 + "\n"
                 ">geneB\nMKKKKKKKKK\n")
    return p
```

```python
# tests/pwmmap/test_dbd.py
from boltzscan.pwmmap import dbd

def test_parse_domtbl_keeps_only_dbd_pfams(tiny_domtbl):
    rows = dbd.parse_domtbl_dbds(tiny_domtbl)
    assert rows == [("geneA", "PF00010", 10, 65)]  # PF99999 dropped

def test_extract_dbds_slices_envelope(tiny_domtbl, tiny_proteins):
    recs = dbd.extract_dbds(tiny_proteins, domtbl=tiny_domtbl)
    assert len(recs) == 1
    r = recs[0]
    assert r.tf_id == "geneA" and r.pfam_acc == "PF00010" and r.family == "bHLH"
    assert len(r.seq) == 65 - 10 + 1   # 1-based inclusive envelope

def test_write_dbd_fasta_id_format(tiny_domtbl, tiny_proteins, tmp_path):
    recs = dbd.extract_dbds(tiny_proteins, domtbl=tiny_domtbl)
    out = tmp_path / "dbd.fasta"
    dbd.write_dbd_fasta(recs, out)
    text = out.read_text()
    assert text.startswith(">geneA__PF00010__0\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/pwmmap/test_dbd.py -q`
Expected: FAIL (module `boltzscan.pwmmap.dbd` not found).

- [ ] **Step 3: Write minimal implementation**

```python
# boltzscan/pwmmap/dbd.py
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/pwmmap/test_dbd.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add boltzscan/pwmmap/dbd.py tests/pwmmap/conftest.py tests/pwmmap/test_dbd.py
git commit -m "feat(pwmmap): DBD extraction from proteins + Pfam domtbl"
```

---

### Task 3: PWM I/O — JASPAR PFM -> cisBP-style txt + MEME

**Files:**
- Create: `boltzscan/pwmmap/pwmio.py`
- Create: `tests/pwmmap/test_pwmio.py`

**Interfaces:**
- Consumes: `boltzscan.utils.io_utils.txt_to_meme` (existing PWM-txt -> MEME converter).
- Produces: `pfm_to_txt(pfm: dict[str, list[float]], path) -> Path` (writes header `Pos\tA\tC\tG\tT` then per-position **probabilities**); `write_txt_and_meme(pfm: dict, motif_id: str, txt_dir, meme_dir) -> tuple[Path, Path]`; `copy_cisbp_pwm(src_txt, motif_id, txt_dir, meme_dir) -> tuple[Path, Path]` (cisBP M*.txt is already in `Pos A C G T` form: copy to txt, convert to meme).

- [ ] **Step 1: Write the failing test**

```python
# tests/pwmmap/test_pwmio.py
from boltzscan.pwmmap import pwmio

def test_pfm_to_txt_normalizes_counts(tmp_path):
    pfm = {"A": [10, 0], "C": [0, 0], "G": [0, 10], "T": [0, 0]}
    p = pwmio.pfm_to_txt(pfm, tmp_path / "MA0001.1.txt")
    lines = p.read_text().strip().splitlines()
    assert lines[0].split() == ["Pos", "A", "C", "G", "T"]
    # position 1 is all-A -> A prob 1.0; position 2 all-G -> G prob 1.0
    assert lines[1].split()[1] == "1.0"
    assert lines[2].split()[3] == "1.0"

def test_write_txt_and_meme_emits_both(tmp_path):
    pfm = {"A": [5, 1], "C": [1, 1], "G": [1, 1], "T": [1, 5]}
    txt, meme = pwmio.write_txt_and_meme(pfm, "MA0001.1", tmp_path/"txt", tmp_path/"meme")
    assert txt.exists() and meme.exists()
    assert "MOTIF MA0001.1" in meme.read_text()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/pwmmap/test_pwmio.py -q`
Expected: FAIL (module not found).

- [ ] **Step 3: Write minimal implementation**

```python
# boltzscan/pwmmap/pwmio.py
"""Convert PWMs between JASPAR PFM (counts), cisBP-style txt (probabilities),
and MEME, reusing the existing txt->MEME converter."""
import shutil
from pathlib import Path

from boltzscan.utils.io_utils import txt_to_meme


def pfm_to_txt(pfm, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["A", "C", "G", "T"]
    width = len(pfm["A"])
    with open(path, "w") as fh:
        fh.write("Pos\tA\tC\tG\tT\n")
        for i in range(width):
            counts = [float(pfm[c][i]) for c in cols]
            tot = sum(counts) or 1.0
            probs = [x / tot for x in counts]
            fh.write(f"{i+1}\t" + "\t".join(repr(round(p, 6)) for p in probs) + "\n")
    return path


def write_txt_and_meme(pfm, motif_id, txt_dir, meme_dir):
    txt_dir, meme_dir = Path(txt_dir), Path(meme_dir)
    txt = pfm_to_txt(pfm, txt_dir / f"{motif_id}.txt")
    meme_dir.mkdir(parents=True, exist_ok=True)
    txt_to_meme(input_path=str(txt), output_path=str(meme_dir / f"{motif_id}.meme"),
                force=True)
    return txt, meme_dir / f"{motif_id}.meme"


def copy_cisbp_pwm(src_txt, motif_id, txt_dir, meme_dir):
    txt_dir, meme_dir = Path(txt_dir), Path(meme_dir)
    txt_dir.mkdir(parents=True, exist_ok=True)
    meme_dir.mkdir(parents=True, exist_ok=True)
    dst_txt = txt_dir / f"{motif_id}.txt"
    shutil.copyfile(src_txt, dst_txt)
    txt_to_meme(input_path=str(dst_txt), output_path=str(meme_dir / f"{motif_id}.meme"),
                force=True)
    return dst_txt, meme_dir / f"{motif_id}.meme"
```

> Note: `txt_to_meme` signature is in `boltzscan/utils/io_utils.py`; confirm its kwargs (`input_path`, `output_path`, `nsites`, `background`, `force`) at implementation time and adapt the two calls above if names differ. Empty/degenerate cisBP txt files (the MEME-format `rose_*_pwm.txt` variant seen in OB_rose) must be detected and skipped with a logged warning rather than crashing the converter.

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/pwmmap/test_pwmio.py -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add boltzscan/pwmmap/pwmio.py tests/pwmmap/test_pwmio.py
git commit -m "feat(pwmmap): PFM/txt/MEME conversion helpers"
```

---

### Task 4: cisBP source — download entire dataset + parse high-quality refs

**Files:**
- Create: `boltzscan/pwmmap/models.py`, `boltzscan/pwmmap/sources/cisbp.py`
- Create: `tests/pwmmap/test_cisbp.py`

**Interfaces:**
- Produces: dataclass `RefTf(ref_id: str, source: str, species: str, family: str, dbid: str, motif_ids: list[str], uniprot_ids: list[str])` (in `models.py`); `download_cisbp(dest, refresh=False) -> tuple[Path, Path]` returning `(tf_information_txt, pwms_dir)` (downloads + unzips the two verified URLs, skips if present); `parse_cisbp_refs(tf_info_path) -> list[RefTf]` keeping rows with a real `Motif_ID` (not `.`/empty) and `TF_Status == 'D'`, grouped by `DBID`.

- [ ] **Step 1: Write the failing test**

```python
# tests/pwmmap/test_cisbp.py
from boltzscan.pwmmap.sources import cisbp

def _fake_tf_info(tmp_path):
    header = ("TF_ID\tFamily_ID\tTSource_ID\tMotif_ID\tMSource_ID\tDBID\tTF_Name\t"
              "TF_Species\tTF_Status\tFamily_Name\tDBDs\tDBD_Count\n")
    rows = [
        "T1\tF1\tS\tM001_3.00\tMS\tAT1G01250\tERF1\tArabidopsis_thaliana\tD\tAP2\tAP2\t1",
        "T2\tF1\tS\tM002_3.00\tMS\tAT1G01260\tERF2\tArabidopsis_thaliana\tD\tAP2\tAP2\t1",
        "T3\tF1\tS\t.\tMS\tAT9G99999\tX\tArabidopsis_thaliana\tD\tAP2\tAP2\t1",     # no motif
        "T4\tF1\tS\tM003_3.00\tMS\tSORGHUM01\tY\tSorghum_bicolor\tI\tAP2\tAP2\t1", # inferred
        "T1\tF1\tS\tM010_3.00\tMS\tAT1G01250\tERF1\tArabidopsis_thaliana\tD\tAP2\tAP2\t1", # 2nd motif for T1
    ]
    p = tmp_path / "TF_Information.txt"
    p.write_text(header + "\n".join(rows) + "\n")
    return p

def test_parse_keeps_determined_with_motif_and_groups(tmp_path):
    refs = {r.dbid: r for r in cisbp.parse_cisbp_refs(_fake_tf_info(tmp_path))}
    assert set(refs) == {"AT1G01250", "AT1G01260"}     # T3 (no motif) & T4 (inferred) dropped
    assert sorted(refs["AT1G01250"].motif_ids) == ["M001_3.00", "M010_3.00"]
    assert refs["AT1G01250"].source == "cisbp"
    assert refs["AT1G01260"].family == "AP2"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/pwmmap/test_cisbp.py -q`
Expected: FAIL (module not found).

- [ ] **Step 3: Write minimal implementation**

```python
# boltzscan/pwmmap/models.py
from dataclasses import dataclass, field

@dataclass
class RefTf:
    ref_id: str
    source: str            # "cisbp" | "jaspar"
    species: str
    family: str
    dbid: str
    motif_ids: list = field(default_factory=list)
    uniprot_ids: list = field(default_factory=list)
```

```python
# boltzscan/pwmmap/sources/cisbp.py
"""Download cisBP entire dataset (TF_Information + PWMs) and parse high-quality refs."""
import io
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

from boltzscan.pwmmap.models import RefTf

BASE = "https://cisbp.ccbr.utoronto.ca/data/3_10/DataFiles/Bulk_downloads/EntireDataset"
TF_INFO_URL = f"{BASE}/TF_Information.zip"
PWMS_URL = f"{BASE}/PWMs.zip"


def _download_zip(url, dest_dir):
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    r = requests.get(url, stream=True, timeout=600)
    r.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        z.extractall(dest_dir)
    return dest_dir


def download_cisbp(dest, refresh=False):
    dest = Path(dest)
    tf_info = next(dest.glob("**/TF_Information.txt"), None)
    pwms_dir = next((p for p in dest.glob("**/pwms_all_motifs") if p.is_dir()), None) \
        or next((p for p in dest.glob("**/pwms*") if p.is_dir()), None)
    if refresh or tf_info is None:
        _download_zip(TF_INFO_URL, dest / "tf_information")
        tf_info = next((dest / "tf_information").glob("**/TF_Information.txt"))
    if refresh or pwms_dir is None:
        _download_zip(PWMS_URL, dest / "pwms")
        pwms_dir = dest / "pwms"
    return tf_info, pwms_dir


def parse_cisbp_refs(tf_info_path):
    df = pd.read_csv(tf_info_path, sep="\t", dtype=str).fillna("")
    df = df[(df["TF_Status"] == "D") & (df["Motif_ID"] != "") & (df["Motif_ID"] != ".")]
    by_dbid = defaultdict(lambda: {"motifs": set(), "family": "", "species": ""})
    for _, row in df.iterrows():
        d = by_dbid[row["DBID"]]
        d["motifs"].add(row["Motif_ID"])
        d["family"] = row["Family_Name"]
        d["species"] = row["TF_Species"]
    return [
        RefTf(ref_id=f"cisbp:{dbid}", source="cisbp", species=v["species"],
              family=v["family"], dbid=dbid, motif_ids=sorted(v["motifs"]))
        for dbid, v in by_dbid.items()
    ]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/pwmmap/test_cisbp.py -q`
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
git add boltzscan/pwmmap/models.py boltzscan/pwmmap/sources/cisbp.py tests/pwmmap/test_cisbp.py
git commit -m "feat(pwmmap): cisBP entire-dataset download + high-quality ref parsing"
```

---

### Task 5: JASPAR source + UniProt sequence resolver

**Files:**
- Create: `boltzscan/pwmmap/sources/jaspar.py`, `boltzscan/pwmmap/sources/uniprot.py`
- Create: `tests/pwmmap/test_jaspar.py`, `tests/pwmmap/test_uniprot.py`

**Interfaces:**
- Produces (jaspar): `iter_jaspar_matrices(collection="CORE", page_size=500) -> Iterator[dict]` (yields detail dicts); `jaspar_refs_and_pwms(dest, txt_dir, meme_dir, refresh=False) -> list[RefTf]` (downloads, writes `<matrix_id>.txt`/`.meme` via `pwmio.write_txt_and_meme`, returns RefTf with `motif_ids=[matrix_id]`, `uniprot_ids`, family).
- Produces (uniprot): `fetch_fasta(acc) -> str | None`; `map_ids_to_uniprot(ids) -> dict[str,str]` (UniProt ID-mapping); `resolve_sequences(refs, cache_path) -> dict[ref_id, str]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/pwmmap/test_jaspar.py
from boltzscan.pwmmap.sources import jaspar
from boltzscan.pwmmap import pwmio

def test_detail_to_ref_and_pwm(tmp_path, monkeypatch):
    detail = {"matrix_id": "MA0570.1", "name": "ABF1", "family": ["bZIP"],
              "species": [{"name": "Arabidopsis thaliana"}], "uniprot_ids": ["P00001"],
              "pfm": {"A": [5,1], "C": [1,1], "G": [1,1], "T": [1,5]}}
    monkeypatch.setattr(jaspar, "iter_jaspar_matrices", lambda **k: iter([detail]))
    refs = jaspar.jaspar_refs_and_pwms(tmp_path/"raw", tmp_path/"txt", tmp_path/"meme")
    assert len(refs) == 1
    r = refs[0]
    assert r.source == "jaspar" and r.motif_ids == ["MA0570.1"]
    assert r.uniprot_ids == ["P00001"]
    assert (tmp_path/"txt"/"MA0570.1.txt").exists()
    assert (tmp_path/"meme"/"MA0570.1.meme").exists()
```

```python
# tests/pwmmap/test_uniprot.py
from boltzscan.pwmmap.sources import uniprot

def test_parse_fasta_body():
    fa = ">sp|P00001|X\nMKTA\nYYYY\n"
    assert uniprot._seq_from_fasta(fa) == "MKTAYYYY"

def test_resolve_uses_uniprot_id_directly(monkeypatch, tmp_path):
    from boltzscan.pwmmap.models import RefTf
    monkeypatch.setattr(uniprot, "fetch_fasta", lambda acc: ">x\nMKK\n" if acc == "P9" else None)
    refs = [RefTf("jaspar:MA1", "jaspar", "sp", "bZIP", "MA1", ["MA1"], ["P9"])]
    seqs = uniprot.resolve_sequences(refs, cache_path=tmp_path/"cache.json")
    assert seqs["jaspar:MA1"] == "MKK"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/pwmmap/test_jaspar.py tests/pwmmap/test_uniprot.py -q`
Expected: FAIL (modules not found).

- [ ] **Step 3: Write minimal implementation**

```python
# boltzscan/pwmmap/sources/jaspar.py
"""Download JASPAR CORE matrices via REST and convert PFMs to txt+meme."""
import json
from pathlib import Path

import requests

from boltzscan.pwmmap.models import RefTf
from boltzscan.pwmmap import pwmio

API = "https://jaspar.elixir.no/api/v1/matrix/"


def iter_jaspar_matrices(collection="CORE", page_size=500):
    url = f"{API}?collection={collection}&page_size={page_size}"
    while url:
        r = requests.get(url, headers={"Accept": "application/json"}, timeout=120)
        r.raise_for_status()
        page = r.json()
        for stub in page["results"]:
            d = requests.get(stub["url"], headers={"Accept": "application/json"},
                             timeout=120)
            d.raise_for_status()
            yield d.json()
        url = page.get("next")


def jaspar_refs_and_pwms(dest, txt_dir, meme_dir, refresh=False, **iter_kwargs):
    dest = Path(dest); dest.mkdir(parents=True, exist_ok=True)
    refs = []
    for det in iter_jaspar_matrices(**iter_kwargs):
        mid = det["matrix_id"]
        pfm = det.get("pfm")
        if not pfm:
            continue
        (dest / f"{mid}.json").write_text(json.dumps(det))
        pwmio.write_txt_and_meme(pfm, mid, txt_dir, meme_dir)
        species = "; ".join(s.get("name", "") for s in det.get("species", [])) or "?"
        fam = (det.get("family") or ["?"])[0]
        refs.append(RefTf(ref_id=f"jaspar:{mid}", source="jaspar", species=species,
                          family=fam, dbid=mid, motif_ids=[mid],
                          uniprot_ids=det.get("uniprot_ids") or []))
    return refs
```

```python
# boltzscan/pwmmap/sources/uniprot.py
"""Resolve reference TF protein sequences from UniProt."""
import json
import time
from pathlib import Path

import requests

FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"
IDMAP = "https://rest.uniprot.org/idmapping"


def _seq_from_fasta(text):
    if not text or not text.startswith(">"):
        return None
    return "".join(text.splitlines()[1:])


def fetch_fasta(acc):
    try:
        r = requests.get(FASTA.format(acc), timeout=60)
        if r.status_code == 200 and r.text.startswith(">"):
            return _seq_from_fasta(r.text)
    except requests.RequestException:
        return None
    return None


def map_ids_to_uniprot(ids, from_db="Ensembl_Genomes", to_db="UniProtKB"):
    """Batch map source IDs -> UniProt accessions; returns {src_id: uniprot_acc}."""
    if not ids:
        return {}
    sub = requests.post(f"{IDMAP}/run",
                        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
                        timeout=120)
    sub.raise_for_status()
    job = sub.json()["jobId"]
    for _ in range(60):
        st = requests.get(f"{IDMAP}/status/{job}", timeout=120).json()
        if st.get("results") is not None or st.get("jobStatus") == "FINISHED":
            break
        time.sleep(3)
    res = requests.get(f"{IDMAP}/results/{job}?size=500", timeout=120).json()
    return {x["from"]: x["to"]["primaryAccession"] if isinstance(x["to"], dict)
            else x["to"] for x in res.get("results", [])}


def resolve_sequences(refs, cache_path):
    cache_path = Path(cache_path)
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    out = {}
    pending_cisbp = []
    for ref in refs:
        if ref.ref_id in cache:
            if cache[ref.ref_id]:
                out[ref.ref_id] = cache[ref.ref_id]
            continue
        acc = ref.uniprot_ids[0] if ref.uniprot_ids else None
        if acc:
            seq = fetch_fasta(acc)
            cache[ref.ref_id] = seq or ""
            if seq:
                out[ref.ref_id] = seq
        else:
            pending_cisbp.append(ref)
    # cisBP refs without a direct UniProt id: batch-map their DBID
    mapped = map_ids_to_uniprot([r.dbid for r in pending_cisbp]) if pending_cisbp else {}
    for r in pending_cisbp:
        acc = mapped.get(r.dbid)
        seq = fetch_fasta(acc) if acc else None
        cache[r.ref_id] = seq or ""
        if seq:
            out[r.ref_id] = seq
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache))
    return out
```

> Note: UniProt ID-mapping `from_db` depends on the cisBP source namespace (Ensembl plant gene IDs, RefSeq, etc.). At implementation, sample cisBP `DBID`s per `TfSource_Name` and pick the right `from` db(s); loop over a few candidate `from_db`s and keep the first that resolves. This is the coverage risk flagged in the spec — log resolved/total.

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/pwmmap/test_jaspar.py tests/pwmmap/test_uniprot.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add boltzscan/pwmmap/sources/jaspar.py boltzscan/pwmmap/sources/uniprot.py \
        tests/pwmmap/test_jaspar.py tests/pwmmap/test_uniprot.py
git commit -m "feat(pwmmap): JASPAR REST source + UniProt sequence resolver"
```

---

### Task 6: Stage A driver — build the cached reference store

**Files:**
- Create: `boltzscan/pwmmap/refs.py`
- Create: `tests/pwmmap/test_refs.py`

**Interfaces:**
- Consumes: `sources.cisbp`, `sources.jaspar`, `sources.uniprot`, `dbd`, `pwmio`.
- Produces: dataclass `RefStore(root: Path, ref_dbd_fasta: Path, ref_index_tsv: Path, motif_txt_dir: Path, motif_meme_dir: Path)`; `build_reference_db(refs_dir="data/pwms/_refs", pfam="data/pfam/Pfam-A.hmm", cpu=8, refresh=False, include_cisbp=True, include_jaspar=True) -> RefStore`. Writes `ref_proteins.fasta`, `ref_dbd.fasta` (IDs `cisbp:DBID__PFxxxxx__i` / `jaspar:MA..__PFxxxxx__i`), `ref_index.tsv` (`ref_id, source, species, family, pfam_acc, dbd_seq_id, motif_ids`), `motif_store/{txt,meme}/`, `build_manifest.json`. Idempotent (skips present artifacts unless `refresh`).
- `load_ref_store(refs_dir) -> RefStore`; `load_ref_index(refs_dir) -> dict[ref_id, dict]`.

- [ ] **Step 1: Write the failing test** (network bypassed via monkeypatch)

```python
# tests/pwmmap/test_refs.py
from boltzscan.pwmmap import refs as R
from boltzscan.pwmmap.models import RefTf

def test_build_reference_db_offline(tmp_path, monkeypatch):
    # Two refs sharing a (real) Pfam DBD so hmmsearch finds them.
    refs = [RefTf("cisbp:G1","cisbp","Arabidopsis","bHLH","G1",["M001_3.00"]),
            RefTf("jaspar:MA1.1","jaspar","Arabidopsis","bHLH","MA1.1",["MA1.1"],["P9"])]
    # a real bHLH DBD sequence (~60 aa) padded; both refs get the same for the test
    bhlh = ("MKRAHHNALERRRRDHIKDSFSSLRDSVPSLQGEKASRAQILDKATEYIQYMRRKNHTHQQDIDDLKRQNALLEQQVRAL")
    monkeypatch.setattr(R, "_gather_cisbp", lambda d, refresh: ([refs[0]], {"M001_3.00": None}))
    monkeypatch.setattr(R, "_gather_jaspar", lambda d, t, m, refresh: [refs[1]])
    monkeypatch.setattr(R, "_resolve_seqs",
                        lambda rr, c: {"cisbp:G1": bhlh, "jaspar:MA1.1": bhlh})
    store = R.build_reference_db(tmp_path/"_refs", include_cisbp=True, include_jaspar=True)
    idx = R.load_ref_index(tmp_path/"_refs")
    assert store.ref_dbd_fasta.exists()
    # both refs produced at least one DBD record tagged with a Pfam acc
    assert any(v["pfam_acc"].startswith("PF") for v in idx.values())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/pwmmap/test_refs.py -q`
Expected: FAIL (module not found).

- [ ] **Step 3: Write minimal implementation**

```python
# boltzscan/pwmmap/refs.py
"""Stage A: build the reusable cisBP+JASPAR reference DBD store (run once)."""
import json
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from boltzscan.pwmmap import dbd, pwmio
from boltzscan.pwmmap.sources import cisbp, jaspar, uniprot


@dataclass
class RefStore:
    root: Path
    ref_dbd_fasta: Path
    ref_index_tsv: Path
    motif_txt_dir: Path
    motif_meme_dir: Path


def _gather_cisbp(dest, refresh):
    tf_info, pwms_dir = cisbp.download_cisbp(dest, refresh=refresh)
    return cisbp.parse_cisbp_refs(tf_info), pwms_dir


def _gather_jaspar(dest, txt_dir, meme_dir, refresh):
    return jaspar.jaspar_refs_and_pwms(dest, txt_dir, meme_dir, refresh=refresh)


def _resolve_seqs(refs, cache):
    return uniprot.resolve_sequences(refs, cache)


def load_ref_store(refs_dir):
    root = Path(refs_dir)
    return RefStore(root, root/"ref_dbd.fasta", root/"ref_index.tsv",
                    root/"motif_store"/"txt", root/"motif_store"/"meme")


def load_ref_index(refs_dir):
    out = {}
    with open(Path(refs_dir)/"ref_index.tsv") as f:
        header = f.readline().rstrip("\n").split("\t")
        for ln in f:
            vals = ln.rstrip("\n").split("\t")
            row = dict(zip(header, vals))
            out[row["dbd_seq_id"]] = row
    return out


def build_reference_db(refs_dir="data/pwms/_refs", pfam="data/pfam/Pfam-A.hmm",
                       cpu=8, refresh=False, include_cisbp=True, include_jaspar=True):
    root = Path(refs_dir); root.mkdir(parents=True, exist_ok=True)
    txt_dir = root/"motif_store"/"txt"; meme_dir = root/"motif_store"/"meme"
    txt_dir.mkdir(parents=True, exist_ok=True); meme_dir.mkdir(parents=True, exist_ok=True)

    refs = []
    if include_cisbp:
        cis_refs, pwms_dir = _gather_cisbp(root/"cisbp", refresh)
        refs += cis_refs
        # copy needed cisBP motif txt -> motif_store, convert to meme
        wanted = {m for r in cis_refs for m in r.motif_ids}
        for m in wanted:
            src = Path(pwms_dir)/f"{m}.txt"
            if src.exists():
                try:
                    pwmio.copy_cisbp_pwm(src, m, txt_dir, meme_dir)
                except Exception as e:  # degenerate cisBP txt
                    print(f"[refs] skip cisBP {m}: {e}")
    if include_jaspar:
        refs += _gather_jaspar(root/"jaspar", txt_dir, meme_dir, refresh)

    seqs = _resolve_seqs(refs, root/"uniprot_cache.json")

    # write reference protein fasta (only refs with a sequence)
    prot_fa = root/"ref_proteins.fasta"
    by_id = {}
    with open(prot_fa, "w") as fh:
        for ref in refs:
            s = seqs.get(ref.ref_id)
            if not s:
                continue
            fh.write(f">{ref.ref_id}\n{s}\n")
            by_id[ref.ref_id] = ref

    # extract reference DBDs (hmmsearch over ref proteins), write ref_dbd.fasta + index
    recs = dbd.extract_dbds(prot_fa, domtbl=None, pfam=pfam, cpu=cpu, work_dir=root)
    dbd.write_dbd_fasta(recs, root/"ref_dbd.fasta")
    with open(root/"ref_index.tsv", "w") as fh:
        fh.write("ref_id\tsource\tspecies\tfamily\tpfam_acc\tdbd_seq_id\tmotif_ids\n")
        for i, r in enumerate(recs):
            ref = by_id.get(r.tf_id)
            if not ref:
                continue
            sid = f"{r.tf_id}__{r.pfam_acc}__{i}"
            fh.write(f"{ref.ref_id}\t{ref.source}\t{ref.species}\t{ref.family}\t"
                     f"{r.pfam_acc}\t{sid}\t{';'.join(ref.motif_ids)}\n")

    (root/"build_manifest.json").write_text(json.dumps({
        "n_refs": len(refs), "n_with_seq": len(by_id), "n_dbd": len(recs),
        "include_cisbp": include_cisbp, "include_jaspar": include_jaspar,
    }, indent=2))
    print(f"[refs] {len(refs)} refs, {len(by_id)} with seq, {len(recs)} DBDs -> {root}")
    return load_ref_store(root)
```

> Note: `write_dbd_fasta` uses a global index `i`; `ref_index.tsv` must use the **same** index when forming `dbd_seq_id`. The implementation above reuses the enumerate index `i` over the same `recs` list, so the IDs match `write_dbd_fasta`. Keep these two in lockstep (or have `write_dbd_fasta` return the list of written IDs and consume that).

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/pwmmap/test_refs.py -q`
Expected: PASS (1 passed). (Requires `hmmsearch` + `data/pfam/Pfam-A.hmm`; the bHLH fixture will hit PF00010.)

- [ ] **Step 5: Commit**

```bash
git add boltzscan/pwmmap/refs.py tests/pwmmap/test_refs.py
git commit -m "feat(pwmmap): Stage A reference-store builder"
```

---

### Task 7: blastp DBD %ID engine

**Files:**
- Create: `boltzscan/pwmmap/align.py`
- Create: `tests/pwmmap/test_align.py`

**Interfaces:**
- Produces: dataclass `Hit(query_tf: str, query_pfam: str, ref_dbd_id: str, ref_pfam: str, pct_id: float)`; `resolve_blast_bins(blastp=None, makeblastdb=None) -> tuple[str,str]`; `make_blast_db(ref_dbd_fasta, db_path, makeblastdb) -> Path`; `parse_blast_rows(rows: list[list[str]], min_cov: float) -> list[Hit]` (pure, testable: `rows` are outfmt-6 fields `qseqid sseqid pident length nident qlen slen`); `blast_dbd_pct_id(query_fasta, db_path, blastp, cpu=8, min_cov=0.8) -> list[Hit]`. Pfam acc is parsed from the `__PF…__` segment of each id; cross-Pfam pairs are dropped.

- [ ] **Step 1: Write the failing test**

```python
# tests/pwmmap/test_align.py
from boltzscan.pwmmap import align

def test_parse_blast_rejects_low_coverage_and_cross_pfam():
    rows = [
        # full-length, same pfam, 50/55 identical -> kept, pct = 50/55
        ["qA__PF00010__0", "rB__PF00010__3", "90.9", "55", "50", "60", "58"],
        # short HSP (len 20 of qlen 60) -> coverage <0.8 -> dropped
        ["qA__PF00010__0", "rC__PF00010__1", "100.0", "20", "20", "60", "58"],
        # cross-pfam -> dropped
        ["qA__PF00010__0", "rD__PF00046__2", "95.0", "55", "52", "60", "58"],
    ]
    hits = align.parse_blast_rows(rows, min_cov=0.8)
    assert len(hits) == 1
    h = hits[0]
    assert h.query_tf == "qA" and h.ref_dbd_id == "rB__PF00010__3"
    assert h.query_pfam == "PF00010" and round(h.pct_id, 4) == round(50/55, 4)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/pwmmap/test_align.py -q`
Expected: FAIL (module not found).

- [ ] **Step 3: Write minimal implementation**

```python
# boltzscan/pwmmap/align.py
"""blastp-based DBD %ID: fast all-vs-all with domain-level, coverage-gated identity."""
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

_SWT = "/home/zlab/miniconda3/envs/swt/bin"
OUTFMT = "6 qseqid sseqid pident length nident qlen slen"


@dataclass
class Hit:
    query_tf: str
    query_pfam: str
    ref_dbd_id: str
    ref_pfam: str
    pct_id: float


def resolve_blast_bins(blastp=None, makeblastdb=None):
    bp = blastp or shutil.which("blastp") or f"{_SWT}/blastp"
    mk = makeblastdb or shutil.which("makeblastdb") or f"{_SWT}/makeblastdb"
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
    subprocess.run([blastp, "-query", str(query_fasta), "-db", str(db_path),
                    "-outfmt", OUTFMT, "-num_threads", str(cpu),
                    "-max_target_seqs", "2000", "-out", out], check=True)
    rows = [ln.split("\t") for ln in Path(out).read_text().splitlines() if ln.strip()]
    Path(out).unlink(missing_ok=True)
    return parse_blast_rows(rows, min_cov)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/pwmmap/test_align.py -q`
Expected: PASS (1 passed).

- [ ] **Step 5: Add a live blastp integration test (marked, optional)**

```python
# append to tests/pwmmap/test_align.py
import pytest
@pytest.mark.network  # not network, but external binary; excluded from default run
def test_real_blast_roundtrip(tmp_path):
    ref = tmp_path/"ref.fasta"; ref.write_text(">rX__PF00010__0\nMKRAHHNALERRRRDHIKDSFSSLRDSVPSLQ\n")
    qry = tmp_path/"q.fasta"; qry.write_text(">qX__PF00010__0\nMKRAHHNALERRRRDHIKDSFSSLRDSVPSLQ\n")
    bp, mk = align.resolve_blast_bins()
    db = align.make_blast_db(ref, tmp_path/"db", mk)
    hits = align.blast_dbd_pct_id(qry, db, bp, min_cov=0.8)
    assert hits and hits[0].pct_id == 1.0
```

- [ ] **Step 6: Commit**

```bash
git add boltzscan/pwmmap/align.py tests/pwmmap/test_align.py
git commit -m "feat(pwmmap): blastp coverage-gated DBD %ID engine"
```

---

### Task 8: Stage B driver — map a species to the reference store

**Files:**
- Create: `boltzscan/pwmmap/mapper.py`
- Create: `tests/pwmmap/test_mapper.py`

**Interfaces:**
- Consumes: `dbd`, `refs.load_ref_store/load_ref_index`, `align`, `thresholds`, `pwmio`.
- Produces: dataclass `MapSummary(out_dir: Path, n_species_tfs: int, n_mapped: int, n_motifs: int)`; `map_species(species_fasta, out_dir, refs_dir="data/pwms/_refs", domtbl=None, threshold_mode="family", threshold=0.70, min_cov=0.8, blastp=None, makeblastdb=None, pfam="data/pfam/Pfam-A.hmm", cpu=8) -> MapSummary`. Writes `<out_dir>/tf2pwms.json`, `txt/`, `meme/`, `map_report.tsv`.

- [ ] **Step 1: Write the failing test** (uses a tiny prebuilt store; bypass blast with monkeypatch)

```python
# tests/pwmmap/test_mapper.py
import json
from boltzscan.pwmmap import mapper, align
from boltzscan.pwmmap.dbd import DbdRecord

def _make_store(tmp_path):
    root = tmp_path/"_refs"; (root/"motif_store"/"txt").mkdir(parents=True)
    (root/"motif_store"/"meme").mkdir(parents=True)
    (root/"ref_dbd.fasta").write_text(">cisbp:G1__PF00010__0\nMKRAHH\n")
    (root/"ref_index.tsv").write_text(
        "ref_id\tsource\tspecies\tfamily\tpfam_acc\tdbd_seq_id\tmotif_ids\n"
        "cisbp:G1\tcisbp\tAth\tbHLH\tPF00010\tcisbp:G1__PF00010__0\tM001_3.00\n")
    for ext, d in (("txt","txt"),("meme","meme")):
        (root/"motif_store"/d/f"M001_3.00.{ext}").write_text("x")
    return root

def test_map_species_transfers_motif_above_threshold(tmp_path, monkeypatch):
    root = _make_store(tmp_path)
    sp = tmp_path/"cm.fasta"; sp.write_text(">CM1\nMKRAHHWXYZ\n")
    # species DBD extraction -> one bHLH DBD on CM1
    monkeypatch.setattr(mapper.dbd, "extract_dbds",
        lambda *a, **k: [DbdRecord("CM1","PF00010","bHLH",1,6,"MKRAHH")])
    # blast -> CM1 matches ref G1 DBD at 0.95 (above bHLH cutoff 0.69)
    monkeypatch.setattr(mapper, "_blast",
        lambda *a, **k: [align.Hit("CM1","PF00010","cisbp:G1__PF00010__0","PF00010",0.95)])
    s = mapper.map_species(sp, tmp_path/"cm_pwms", refs_dir=root, threshold_mode="family")
    j = json.loads((tmp_path/"cm_pwms"/"tf2pwms.json").read_text())
    assert j["CM1"] == ["M001_3.00"]
    assert (tmp_path/"cm_pwms"/"txt"/"M001_3.00.txt").exists()
    assert s.n_mapped == 1

def test_below_threshold_not_transferred(tmp_path, monkeypatch):
    root = _make_store(tmp_path)
    sp = tmp_path/"cm.fasta"; sp.write_text(">CM1\nMKRAHHWXYZ\n")
    monkeypatch.setattr(mapper.dbd, "extract_dbds",
        lambda *a, **k: [DbdRecord("CM1","PF00010","bHLH",1,6,"MKRAHH")])
    monkeypatch.setattr(mapper, "_blast",
        lambda *a, **k: [align.Hit("CM1","PF00010","cisbp:G1__PF00010__0","PF00010",0.40)])
    s = mapper.map_species(sp, tmp_path/"cm_pwms", refs_dir=root, threshold_mode="family")
    assert s.n_mapped == 0
    assert json.loads((tmp_path/"cm_pwms"/"tf2pwms.json").read_text()) == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/pwmmap/test_mapper.py -q`
Expected: FAIL (module not found).

- [ ] **Step 3: Write minimal implementation**

```python
# boltzscan/pwmmap/mapper.py
"""Stage B: map a species' TFs to reference motifs via DBD %ID (run per species)."""
import json
import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from boltzscan.pwmmap import dbd, align
from boltzscan.pwmmap.refs import load_ref_store, load_ref_index
from boltzscan.pwmmap.thresholds import cutoff_for


@dataclass
class MapSummary:
    out_dir: Path
    n_species_tfs: int
    n_mapped: int
    n_motifs: int


def _blast(query_fasta, store, blastp, makeblastdb, cpu, min_cov, work_dir):
    db = align.make_blast_db(store.ref_dbd_fasta, Path(work_dir)/"refdb", makeblastdb)
    return align.blast_dbd_pct_id(query_fasta, db, blastp, cpu=cpu, min_cov=min_cov)


def map_species(species_fasta, out_dir, refs_dir="data/pwms/_refs", domtbl=None,
                threshold_mode="family", threshold=0.70, min_cov=0.8,
                blastp=None, makeblastdb=None, pfam="data/pfam/Pfam-A.hmm", cpu=8):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    txt_dir = out_dir/"txt"; meme_dir = out_dir/"meme"
    txt_dir.mkdir(exist_ok=True); meme_dir.mkdir(exist_ok=True)
    store = load_ref_store(refs_dir)
    index = load_ref_index(refs_dir)            # dbd_seq_id -> row(family, motif_ids,...)
    bp, mk = align.resolve_blast_bins(blastp, makeblastdb)

    recs = dbd.extract_dbds(species_fasta, domtbl=domtbl, pfam=pfam, cpu=cpu, work_dir=out_dir)
    sp_dbd_fa = out_dir/"species_dbd.fasta"
    dbd.write_dbd_fasta(recs, sp_dbd_fa)
    n_species_tfs = len({r.tf_id for r in recs})

    hits = _blast(sp_dbd_fa, store, bp, mk, cpu, min_cov, out_dir)

    tf2motifs = defaultdict(set)
    report = [("species_tf", "pfam_acc", "ref_id", "source", "species", "pct_id", "motif_id")]
    for h in hits:
        row = index.get(h.ref_dbd_id)
        if not row:
            continue
        cut = cutoff_for(h.query_pfam, mode=threshold_mode, global_thr=threshold)
        if h.pct_id < cut:
            continue
        for m in row["motif_ids"].split(";"):
            if not m:
                continue
            tf2motifs[h.query_tf].add(m)
            report.append((h.query_tf, h.query_pfam, row["ref_id"], row["source"],
                           row["species"], f"{h.pct_id:.3f}", m))

    # write tf2pwms.json + copy matched motif files from the store
    tf2pwms = {tf: sorted(ms) for tf, ms in tf2motifs.items()}
    (out_dir/"tf2pwms.json").write_text(json.dumps(tf2pwms, indent=2))
    needed = {m for ms in tf2pwms.values() for m in ms}
    for m in needed:
        for d, dst in ((store.motif_txt_dir, txt_dir), (store.motif_meme_dir, meme_dir)):
            for ext in (".txt", ".meme"):
                src = Path(d)/f"{m}{ext}"
                if src.exists():
                    shutil.copyfile(src, dst/f"{m}{ext}")
    with open(out_dir/"map_report.tsv", "w") as fh:
        for row in report:
            fh.write("\t".join(map(str, row)) + "\n")

    return MapSummary(out_dir, n_species_tfs, len(tf2pwms), len(needed))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/pwmmap/test_mapper.py -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add boltzscan/pwmmap/mapper.py tests/pwmmap/test_mapper.py
git commit -m "feat(pwmmap): Stage B species->reference DBD mapping"
```

---

### Task 9: CLI wiring + docs

**Files:**
- Modify: `boltzscan/cli.py` (add two subparsers near `find-tf`/`cistarg`; two handlers; two dispatch entries)
- Modify: `CLAUDE.md` (Common commands section)
- Create: `tests/pwmmap/test_cli.py`

**Interfaces:**
- Consumes: `refs.build_reference_db`, `mapper.map_species`.
- Produces: CLI commands `build-pwm-refs` and `map-pwm` (flags exactly as in the spec's CLI section).

- [ ] **Step 1: Write the failing test**

```python
# tests/pwmmap/test_cli.py
from boltzscan.cli import _build_parser

def test_map_pwm_parses():
    p = _build_parser()
    ns = p.parse_args(["map-pwm", "-f", "cm.fasta", "-o", "data/pwms/cm_pwms",
                       "--domtbl", "x/pfam.domtbl", "--threshold-mode", "family"])
    assert ns.command == "map-pwm" and ns.proteins == "cm.fasta"

def test_build_pwm_refs_parses():
    p = _build_parser()
    ns = p.parse_args(["build-pwm-refs", "--refs", "data/pwms/_refs", "-c", "16"])
    assert ns.command == "build-pwm-refs" and ns.cpu == 16
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/pwmmap/test_cli.py -q`
Expected: FAIL (unknown command `map-pwm`).

- [ ] **Step 3: Add subparsers (in `_build_parser`, after the `find-tf` block)**

```python
    # build-pwm-refs
    p = sub_parsers.add_parser('build-pwm-refs',
        help='Stage A: download cisBP+JASPAR, extract reference DBDs, cache the reference store')
    p.add_argument('--refs', default='data/pwms/_refs', help='Reference store dir (default: data/pwms/_refs)')
    p.add_argument('-m', '--pfam', default='data/pfam/Pfam-A.hmm')
    p.add_argument('-c', '--cpu', type=int, default=8)
    p.add_argument('--no-cisbp', action='store_true')
    p.add_argument('--no-jaspar', action='store_true')
    p.add_argument('--refresh', action='store_true', help='Re-download even if present')

    # map-pwm
    p = sub_parsers.add_parser('map-pwm',
        help='Stage B: map a species TF FASTA to reference motifs via DBD %ID (no network)')
    p.add_argument('-f', '--proteins', required=True, help='Species TF protein FASTA')
    p.add_argument('-o', '--output', required=True, help='Output dir, e.g. data/pwms/cm_pwms')
    p.add_argument('--refs', default='data/pwms/_refs')
    p.add_argument('--domtbl', default=None, help='Reuse a find-tf pfam.domtbl (else run hmmsearch)')
    p.add_argument('--threshold-mode', default='family', choices=['family', 'global'])
    p.add_argument('--threshold', type=float, default=0.70, help='Global %ID cutoff (global mode)')
    p.add_argument('--min-cov', type=float, default=0.8, help='Min DBD coverage for a blast hit')
    p.add_argument('--blastp', default=None)
    p.add_argument('--makeblastdb', default=None)
    p.add_argument('-m', '--pfam', default='data/pfam/Pfam-A.hmm')
    p.add_argument('-c', '--cpu', type=int, default=8)
```

- [ ] **Step 4: Add handlers + dispatch entries**

```python
def _cmd_build_pwm_refs(args):
    from boltzscan.pwmmap.refs import build_reference_db
    store = build_reference_db(refs_dir=args.refs, pfam=args.pfam, cpu=args.cpu,
                               refresh=args.refresh,
                               include_cisbp=not args.no_cisbp,
                               include_jaspar=not args.no_jaspar)
    print(f"Reference store ready at {store.root}")


def _cmd_map_pwm(args):
    from boltzscan.pwmmap.mapper import map_species
    s = map_species(species_fasta=args.proteins, out_dir=args.output, refs_dir=args.refs,
                    domtbl=args.domtbl, threshold_mode=args.threshold_mode,
                    threshold=args.threshold, min_cov=args.min_cov,
                    blastp=args.blastp, makeblastdb=args.makeblastdb,
                    pfam=args.pfam, cpu=args.cpu)
    print(f"Wrote {s.out_dir}/tf2pwms.json "
          f"({s.n_mapped}/{s.n_species_tfs} TFs mapped, {s.n_motifs} motifs)")
```

Add to `_DISPATCH`:

```python
    'build-pwm-refs': _cmd_build_pwm_refs,
    'map-pwm': _cmd_map_pwm,
```

- [ ] **Step 5: Run test + help to verify**

Run: `python -m pytest tests/pwmmap/test_cli.py -q && python bscan.py map-pwm --help`
Expected: tests PASS; help prints the map-pwm usage.

- [ ] **Step 6: Document in CLAUDE.md**

Add under "Common commands" (after the `find-tf` block):

```bash
# Cross-species PWM DB via DBD homology (two stages)
# Stage A (once): download cisBP + JASPAR, cache reference DBD store
python bscan.py build-pwm-refs --refs data/pwms/_refs -c 16
# Stage B (per species): map a TF FASTA to reference motifs -> OB_rose-style DB
python bscan.py map-pwm -f <species_tfs.fasta> -o data/pwms/<sp>_pwms \
    --domtbl <tf_out>/pfam.domtbl --refs data/pwms/_refs
```

- [ ] **Step 7: Commit**

```bash
git add boltzscan/cli.py CLAUDE.md tests/pwmmap/test_cli.py
git commit -m "feat(pwmmap): wire build-pwm-refs + map-pwm CLI; docs"
```

---

### Task 10: End-to-end validation on 野菊 (manual, gated)

**Files:** none (operational run + spot-check)

This task is the real-data validation; it depends on Stage A succeeding (network + UniProt coverage), so run it interactively and inspect before trusting output.

- [ ] **Step 1: Build the reference store (long; network).** Run:
`python bscan.py build-pwm-refs --refs data/pwms/_refs -c 16`
Expected: `build_manifest.json` shows non-trivial `n_with_seq` and `n_dbd`; inspect resolved/total coverage. If cisBP UniProt mapping coverage is poor, revisit `uniprot.map_ids_to_uniprot` `from_db` (see Task 5 note) before proceeding.

- [ ] **Step 2: Map a 50-TF 野菊 subset first.** Run:
```bash
head -100 tasks/zyf/野菊基因组信息/tf_out/tf_proteins.fasta > /tmp/cm50.fasta
python bscan.py map-pwm -f /tmp/cm50.fasta -o /tmp/cm50_pwms --refs data/pwms/_refs
```
Expected: `tf2pwms.json` non-empty for conserved-family TFs (bHLH/MYB/Homeobox); `map_report.tsv` shows sane %IDs ≥ family cutoff.

- [ ] **Step 3: Full 野菊 run.** Run:
```bash
python bscan.py map-pwm -f tasks/zyf/野菊基因组信息/tf_out/tf_proteins.fasta \
    -o data/pwms/cm_pwms --domtbl tasks/zyf/野菊基因组信息/tf_out/pfam.domtbl --refs data/pwms/_refs
```
Expected: `data/pwms/cm_pwms/{tf2pwms.json, txt/, meme/, map_report.tsv}` in OB_rose shape. Spot-check: average motifs/TF and per-family hit rate are in a plausible range vs `data/pwms/OB_rose_pwms/tf2pwms.json` (969 TFs). Confirm a known family (e.g. a 野菊 WRKY) maps only to WRKY motifs.

- [ ] **Step 4 (no commit — data artifacts).** Per repo policy, do not stage `data/pwms/_refs` or `data/pwms/cm_pwms` (large/derived). Report counts to the user.

---

## Self-Review

**Spec coverage:**
- Two-stage reusable architecture → Tasks 6 (A) + 8 (B), CLI Task 9. ✓
- cisBP entire-dataset download + high-quality filter → Task 4. ✓
- JASPAR REST + PFM conversion → Tasks 3, 5. ✓
- UniProt sequence resolution (both sources) → Task 5. ✓
- DBD extraction, shared both sides, same Pfam HMMs → Task 2 (used by 6 + 8). ✓
- blastp standard+fast, coverage-gated DBD %ID, same-Pfam only → Task 7. ✓
- Family-specific thresholds, cross-kingdom safety via Pfam keying → Tasks 1 + 8. ✓
- OB_rose output shape (tf2pwms.json + txt/ + meme/) + audit report → Task 8. ✓
- CLI + docs → Task 9. End-to-end on 野菊 → Task 10. ✓

**Placeholder scan:** No "TBD"/"add error handling"-style steps; the two `> Note:` blocks call out genuine implementation-time decisions (txt_to_meme kwargs; UniProt `from_db`) with concrete guidance, not placeholders.

**Type consistency:** `RefTf` (models.py) shared by sources + refs; `DbdRecord` (dbd.py) used by refs + mapper; `Hit` (align.py) used by mapper; DBD-id format `tf__PFxxxxx__i` produced by `write_dbd_fasta` and parsed by `align._pfam_of/_tf_of` and `refs` index — consistent. `cutoff_for(pfam_acc, mode, global_thr)` signature identical in Task 1 def and Task 8 call.

**Known risk carried from spec:** cisBP→UniProt coverage (Task 5 note + Task 10 Step 1 gate). This is validated before the full run rather than assumed.
