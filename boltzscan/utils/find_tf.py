"""Identify plant transcription factors from a protein FASTA (PlantTFDB-style).

Pipeline (sequence-only, no folding / GPU):
  1. hmmsearch --cut_ga of every protein against Pfam-A  -> domtblout
  2. parse domtblout into per-protein {pfam_acc: (count, best_domain_score)}
  3. apply an ordered PlantTFDB-style rule table to assign a TF family
     (combinatorial rules first, single-domain families next, promiscuous
      zinc-fingers last so they never steal a protein that belongs to a
      co-domain family).

Outputs (in --out-dir):
  pfam.domtbl          raw hmmsearch result (reused unless --force)
  tf_annotation.tsv    one row per TF: gene_id, family, dbd_pfam, evidence, score
  tf_proteins.fasta    the TF protein sequences (handy as boltzscan --tf-fasta)
  + a per-family count summary on stderr.

The DBD->Pfam rule set extends the curated table in scratch/rose_pfam_parse.py.
Several PlantTFDB splits are NOT cleanly separable with stock Pfam-A
(G2-like vs MYB_related, NF-YB vs NF-YC, WOX/HD-ZIP-I/II inside HB-other);
those are labelled at the resolvable granularity and noted in CAVEATS below.

boltzscan env. Usage:
  python boltzscan/utils/find_tf.py \
      --proteins "tasks/zyf/野菊基因组信息/protein_coreset.fa" \
      --pfam data/pfam/Pfam-A.hmm \
      --out-dir "tasks/zyf/野菊基因组信息/tf_out" --cpu 8
"""
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from boltzscan.toolchain import resolve_executable

from Bio import SeqIO

DEFAULT_PFAM = "data/pfam/Pfam-A.hmm"

# ---------------------------------------------------------------------------
# Pfam accessions used by the rules (base accession, no .version)
# ---------------------------------------------------------------------------
# DNA-binding domains (the "is this a TF" evidence set; from rose_pfam_parse.py)
B3, AP2, AUXRESP = "PF02362", "PF00847", "PF06507"
SRF, KBOX = "PF00319", "PF01486"
HOMEO, HOMEOKN = "PF00046", "PF05920"
KNOX1, KNOX2, ELK, PHD, START = "PF03790", "PF03791", "PF03789", "PF00628", "PF01852"
MYB, MYB4, RESPREG = "PF00249", "PF13837", "PF00072"
BBOX, CCT = "PF00643", "PF06203"
WRC, QLQ = "PF08879", "PF08880"
NFYA, NFYBC = "PF02045", "PF00808"

# All Pfam accs that, on their own, count as a plant-TF DNA-binding domain.
DBD_ACCS = {
    HOMEO, HOMEOKN, MYB, MYB4, "PF00010", AP2, "PF02365", B3,
    "PF00096", "PF13894", "PF13912", "PF13465", "PF03106", "PF03514",
    "PF03101", "PF00642", "PF00170", "PF07716", SRF, "PF03195", "PF02701",
    "PF03634", "PF00447", "PF03110", "PF00320", "PF02319", NFYBC, NFYA,
    WRC, "PF04770", "PF02042", BBOX, CCT, "PF03859", "PF04690", "PF04504",
    "PF03638", "PF04873", "PF06943", "PF08536", "PF01422", "PF04689",
    "PF01698", "PF17538", "PF05142", "PF05687", "PF06217", "PF28235",
}

# ---------------------------------------------------------------------------
# PlantTFDB-style family rules, applied top-to-bottom; FIRST match wins.
# Each rule: (family, all_of, none_of)
#   all_of : list of requirements; each item is either
#              "PFxxxxx"            -> present (>=1 copy)
#              ("PFxxxxx", n)       -> present with >= n copies
#              ["PFa","PFb", ...]   -> any one of these present
#   none_of: list of Pfam accs that must be ABSENT
# ---------------------------------------------------------------------------
RULES = [
    # --- AP2/ERF & B3 superfamilies (order matters: combos before singles) ---
    ("RAV",          [AP2, B3],                 []),
    ("AP2",          [(AP2, 2)],                [B3]),
    ("ERF",          [AP2],                     [B3]),
    ("ARF",          [B3, AUXRESP],             []),
    ("B3",           [B3],                      [AP2, AUXRESP]),
    # --- MADS ---
    ("MIKC_MADS",    [SRF, KBOX],               []),
    ("M-type_MADS",  [SRF],                     [KBOX]),
    # --- Homeobox superfamily ---
    ("TALE",         [[HOMEOKN, KNOX1, KNOX2, ELK]], []),  # KNOX/BELL
    ("HB-PHD",       [HOMEO, PHD],              []),
    ("HD-ZIP",       [HOMEO, START],            []),        # HD-ZIP III/IV
    ("HB-other",     [HOMEO],                   []),        # incl. HD-ZIP I/II, WOX
    # --- MYB superfamily ---
    ("ARR-B",        [RESPREG, [MYB, MYB4]],    []),
    ("MYB",          [(MYB, 2)],                []),         # >=2 Myb repeats (R2R3 etc.)
    ("MYB_related",  [[MYB, MYB4]],             []),         # incl. G2-like, Trihelix*
    # --- other combinatorial ---
    ("GRF",          [WRC],                     []),         # WRC (+QLQ)
    ("CO-like",      [BBOX, CCT],               []),
    ("DBB",          [BBOX],                    [CCT]),
    ("NF-YA",        [NFYA],                    []),
    ("NF-YB/NF-YC",  [NFYBC],                   []),         # not split by Pfam
    # --- single clean DBD families ---
    ("bHLH",         ["PF00010"],               []),
    ("bZIP",         [["PF00170", "PF07716"]],  []),
    ("NAC",          ["PF02365"],               []),
    ("WRKY",         ["PF03106"],               []),
    ("Dof",          ["PF02701"],               []),
    ("GATA",         ["PF00320"],               []),
    ("GRAS",         ["PF03514"],               []),
    ("TCP",          ["PF03634"],               []),
    ("HSF",          ["PF00447"],               []),
    ("SBP",          ["PF03110"],               []),
    ("FAR1",         ["PF03101"],               []),
    ("LBD",          ["PF03195"],               []),
    ("E2F/DP",       ["PF02319"],               []),
    ("ZF-HD",        ["PF04770"],               []),
    ("Nin-like",     ["PF02042"],               []),
    ("CAMTA",        ["PF03859"],               []),
    ("YABBY",        ["PF04690"],               []),
    ("GeBP",         ["PF04504"],               []),
    ("CPP",          ["PF03638"],               []),
    ("EIL",          ["PF04873"],               []),
    ("LSD",          ["PF06943"],               []),
    ("Whirly",       ["PF08536"],               []),
    ("NF-X1",        ["PF01422"],               []),
    ("S1Fa-like",    ["PF04689"],               []),
    ("LFY",          [["PF01698", "PF17538"]],  []),
    ("SRS",          ["PF05142"],               []),
    ("BES1",         ["PF05687"],               []),
    ("BBR-BPC",      ["PF06217"],               []),
    ("VOZ",          ["PF28235"],               []),
    # --- promiscuous zinc fingers LAST (common in non-TF proteins) ---
    ("C3H",          ["PF00642"],               []),
    ("C2H2",         [["PF00096", "PF13894", "PF13912", "PF13465"]], []),
]

# ---------------------------------------------------------------------------
# Transposon / transposable-element domains. A protein hitting ANY of these is
# a TE-derived ORF, not a curated TF, and is dropped EVEN IF it carries a Myb /
# SANT-like DBD (the classic Harbinger/Helitron false-positive). The ONE
# exception is FAR1: that family's DBD is itself domesticated from a Mutator
# transposase, so FAR1 calls are kept. MULE (PF10551) is deliberately absent
# here because it is a genuine part of FAR1/FHY3.
# ---------------------------------------------------------------------------
TE_FORBIDDEN = {
    "PF10545",  # MADF_DNA_bdg   (PIF/Harbinger transposase)
    "PF14214",  # Helitron_like_N
    "PF04827",  # Plant_tran     (plant transposase)
    "PF04937",  # DUF659         (transposon)
    "PF00078",  # RVT_1          reverse transcriptase
    "PF07727",  # RVT_2
    "PF13456",  # RVT_3
    "PF00665",  # rve            integrase core
    "PF13683",  # rve_3
    "PF03732",  # Retrotrans_gag
    "PF14223",  # Retrotran_gag_2
    "PF13975",  # gag-asp_proteas
    "PF08284",  # RVP_2          retroviral aspartyl protease
    "PF00075",  # RNase_H
}
TE_EXEMPT_FAMILIES = {"FAR1"}


def base_acc(a):
    return a.split(".")[0]


def run_hmmsearch(proteins, pfam, domtbl, cpu, force):
    if domtbl.exists() and not force:
        print(f"[hmmsearch] reuse existing {domtbl}", file=sys.stderr)
        return
    domtbl.parent.mkdir(parents=True, exist_ok=True)
    hmmsearch = resolve_executable("hmmsearch")
    if hmmsearch is None:
        raise FileNotFoundError(
            "HMMER hmmsearch not found; run `boltzscan doctor --fix`"
        )
    cmd = [
        hmmsearch, "--cut_ga", "--cpu", str(cpu),
        "--domtblout", str(domtbl), "-o", "/dev/null",
        str(pfam), str(proteins),
    ]
    print("[hmmsearch] " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def parse_domtbl(domtbl):
    """protein -> {pfam_acc: [count, best_domain_score]}, and acc -> pfam name."""
    per_prot = defaultdict(lambda: defaultdict(lambda: [0, 0.0]))
    names = {}
    with open(domtbl) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            c = ln.split()
            prot = c[0]                 # target name = our protein
            pname = c[3]                # query name = Pfam model name
            pacc = base_acc(c[4])       # query accession = PFxxxxx
            dscore = float(c[13])       # this-domain bit score
            names[pacc] = pname
            rec = per_prot[prot][pacc]
            rec[0] += 1
            rec[1] = max(rec[1], dscore)
    return per_prot, names


def _present(req, counts):
    """Is a single all_of requirement satisfied by this protein's pfam counts?"""
    if isinstance(req, tuple):                 # ("PFxxxxx", min_count)
        acc, n = req
        return counts.get(acc, [0])[0] >= n
    if isinstance(req, list):                  # any-of
        return any(a in counts for a in req)
    return req in counts                       # single acc, >=1 copy


def assign_family(counts):
    """Return (family, deciding_pfam_acc) or (None, None) if not a TF."""
    for fam, all_of, none_of in RULES:
        if any(a in counts for a in none_of):
            continue
        if all(_present(r, counts) for r in all_of):
            # deciding pfam = highest-scoring DBD acc present (for evidence)
            dbd_present = [a for a in counts if a in DBD_ACCS]
            if dbd_present:
                deciding = max(dbd_present, key=lambda a: counts[a][1])
            else:
                deciding = None
            return fam, deciding
    return None, None


@dataclass
class FindTfSummary:
    annotation_tsv: Path
    tf_fasta: Path
    domtbl: Path
    n_total: int
    n_tf: int
    n_te_excluded: int
    family_counts: dict = field(default_factory=dict)


def find_transcription_factors(proteins, out_dir, pfam=DEFAULT_PFAM,
                               cpu=8, force=False, domtbl=None):
    """Identify plant TFs in a protein FASTA. Returns :class:`FindTfSummary`.

    If ``domtbl`` is given, that existing hmmsearch domtblout is reused and no
    new search is run; otherwise hmmsearch is run into ``<out_dir>/pfam.domtbl``.
    """
    proteins = Path(proteins)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    domtbl = Path(domtbl) if domtbl else out / "pfam.domtbl"
    if not (domtbl.exists() and not force):
        run_hmmsearch(proteins, Path(pfam), domtbl, cpu, force)
    elif domtbl.exists():
        print(f"[hmmsearch] reuse existing {domtbl}", file=sys.stderr)

    per_prot, names = parse_domtbl(domtbl)
    protein_ids = {
        rec.id for rec in SeqIO.parse(str(proteins), "fasta")
    }
    per_prot = {
        prot: counts for prot, counts in per_prot.items()
        if prot in protein_ids
    }

    tf_rows = []          # (gene_id, family, dbd_pfam, dbd_name, dbd_score, evidence)
    fam_count = defaultdict(int)
    n_te_excluded = 0
    for prot, counts in per_prot.items():
        fam, deciding = assign_family(counts)
        if fam is None:
            continue
        if fam not in TE_EXEMPT_FAMILIES and any(a in TE_FORBIDDEN for a in counts):
            n_te_excluded += 1
            continue
        fam_count[fam] += 1
        ev = ";".join(
            f"{a}({names.get(a, a)})x{counts[a][0]}"
            for a in sorted(counts, key=lambda a: -counts[a][1])
        )
        dname = names.get(deciding, "") if deciding else ""
        dscore = f"{counts[deciding][1]:.1f}" if deciding else "0.0"
        tf_rows.append((prot, fam, deciding or "", dname, dscore, ev))

    tf_ids = {r[0] for r in tf_rows}

    ann = out / "tf_annotation.tsv"
    with open(ann, "w") as fh:
        fh.write("gene_id\tfamily\tdbd_pfam\tdbd_name\tdbd_score\tall_pfam_evidence\n")
        for r in sorted(tf_rows, key=lambda r: (r[1], r[0])):
            fh.write("\t".join(map(str, r)) + "\n")

    tf_fasta = out / "tf_proteins.fasta"
    n_total = 0
    n_written = 0
    with open(tf_fasta, "w") as fh:
        for rec in SeqIO.parse(str(proteins), "fasta"):
            n_total += 1
            if rec.id in tf_ids:
                SeqIO.write(rec, fh, "fasta")
                n_written += 1

    return FindTfSummary(
        annotation_tsv=ann,
        tf_fasta=tf_fasta,
        domtbl=domtbl,
        n_total=n_total,
        n_tf=len(tf_ids),
        n_te_excluded=n_te_excluded,
        family_counts=dict(fam_count),
    )


def print_report(summary):
    pct = 100 * summary.n_tf / summary.n_total if summary.n_total else 0.0
    print(f"\n[find_tf] {summary.n_tf}/{summary.n_total} proteins classified as TF "
          f"({pct:.1f}%); {summary.n_te_excluded} DBD-bearing proteins dropped "
          f"as transposon-derived", file=sys.stderr)
    print(f"[find_tf] wrote {summary.annotation_tsv}", file=sys.stderr)
    print(f"[find_tf] wrote {summary.tf_fasta}", file=sys.stderr)
    print("\n[find_tf] per-family counts:", file=sys.stderr)
    for fam in sorted(summary.family_counts, key=lambda k: -summary.family_counts[k]):
        print(f"  {fam:14} {summary.family_counts[fam]}", file=sys.stderr)
