"""Reference-level motif clustering (build our own non-redundant clustered set).

Mirrors the cisTarget/SCENIC+ "nr_clust" idea: collapse redundant PWMs into
clusters and keep one representative each, so a TF that maps to dozens of
near-duplicate reference motifs is reduced to its distinct binding preferences.

Per DBD family (from ref_index) we run Tomtom self-comparison and greedily
cluster: motifs are processed most-supported first; a motif joins the most-
supported already-chosen representative it is similar to (Tomtom q < qthresh),
else it becomes a new representative. The result is cached as
``_refs/motif_clusters.tsv`` (motif_id, family, representative_id) and consumed
by Stage B's ``--collapse-clusters``.
"""
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

_MEME_ENV_TOMTOM = "/home/zlab/miniconda3/envs/meme/bin/tomtom"
_MEME_HEADER = (
    "MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
    "Background letter frequencies (from uniform background):\n"
    "A 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n"
)


@dataclass
class ClusterSummary:
    clusters_tsv: Path
    n_motifs: int
    n_clusters: int
    n_families: int


def resolve_tomtom(tomtom=None):
    return tomtom or shutil.which("tomtom") or _MEME_ENV_TOMTOM


def _safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name) or "NA"


def build_motif_meta(ref_index_tsv):
    """Return (motif_family, motif_support) from ref_index.tsv.

    motif_family[m] = most common family among the ref TFs carrying motif m;
    motif_support[m] = number of distinct ref TFs carrying motif m.
    """
    fam = defaultdict(Counter)
    support = defaultdict(set)
    with open(ref_index_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        fi, mi, ri = header.index("family"), header.index("motif_ids"), header.index("ref_id")
        for ln in f:
            c = ln.rstrip("\n").split("\t")
            if len(c) <= mi:
                continue
            for m in c[mi].split(";"):
                if m:
                    fam[m][c[fi]] += 1
                    support[m].add(c[ri])
    motif_family = {m: cnt.most_common(1)[0][0] for m, cnt in fam.items()}
    motif_support = {m: len(s) for m, s in support.items()}
    return motif_family, motif_support


def _combined_meme(motif_ids, meme_dir, out_path):
    """Concatenate per-motif MEME files into one MEME db (single header)."""
    meme_dir = Path(meme_dir)
    with open(out_path, "w") as out:
        out.write(_MEME_HEADER)
        for m in motif_ids:
            p = meme_dir / f"{m}.meme"
            if not p.exists():
                continue
            text = p.read_text()
            idx = text.find("MOTIF")
            if idx >= 0:
                out.write(text[idx:].rstrip() + "\n\n")
    return out_path


def parse_tomtom_edges(tsv_text):
    """Parse `tomtom -text` output into a symmetric adjacency dict (excl. self)."""
    adj = defaultdict(set)
    for ln in tsv_text.splitlines():
        if not ln or ln.startswith("#") or ln.startswith("Query_ID"):
            continue
        c = ln.split("\t")
        if len(c) < 2:
            continue
        q, t = c[0], c[1]
        if q and t and q != t:
            adj[q].add(t)
            adj[t].add(q)
    return adj


def run_tomtom_self(meme_file, tomtom, qthresh):
    proc = subprocess.run(
        [tomtom, "-text", "-no-ssc", "-thresh", str(qthresh), str(meme_file), str(meme_file)],
        capture_output=True, text=True,
    )
    return parse_tomtom_edges(proc.stdout)


def greedy_cluster(motif_ids, support, adj):
    """Greedy single-pass clustering. Returns {motif_id: representative_id}.

    Motifs are visited most-supported first; each joins the highest-support
    representative it is adjacent to, else starts a new cluster as its own rep.
    """
    order = sorted(motif_ids, key=lambda m: (-support.get(m, 0), m))
    reps = []
    rep_of = {}
    for m in order:
        cand = [r for r in reps if r in adj.get(m, ())]
        if cand:
            rep_of[m] = max(cand, key=lambda r: support.get(r, 0))
        else:
            reps.append(m)
            rep_of[m] = m
    return rep_of


def _cluster_one_bucket(fam, motifs, meme_dir, work, tomtom, qthresh, support):
    if len(motifs) == 1:
        return {motifs[0]: motifs[0]}
    combined = _combined_meme(motifs, meme_dir, work / f"{_safe(fam)}.meme")
    adj = run_tomtom_self(combined, tomtom, qthresh)
    return greedy_cluster(motifs, support, adj)


def cluster_reference_motifs(refs_dir, tomtom=None, qthresh=0.05, cpu=8):
    refs_dir = Path(refs_dir)
    ref_index = refs_dir / "ref_index.tsv"
    meme_dir = refs_dir / "motif_store" / "meme"
    tomtom = resolve_tomtom(tomtom)

    motif_family, motif_support = build_motif_meta(ref_index)
    have_meme = [m for m in motif_family if (meme_dir / f"{m}.meme").exists()]

    buckets = defaultdict(list)
    for m in have_meme:
        buckets[motif_family[m]].append(m)

    work = refs_dir / "_cluster_work"
    work.mkdir(parents=True, exist_ok=True)
    # Process families concurrently; each runs its own (single-threaded) tomtom.
    motif2rep = {}
    with ThreadPoolExecutor(max_workers=cpu) as ex:
        futs = [ex.submit(_cluster_one_bucket, fam, motifs, meme_dir, work,
                          tomtom, qthresh, motif_support)
                for fam, motifs in buckets.items()]
        for fut in futs:
            motif2rep.update(fut.result())

    # any motif without a meme (or family) is its own representative
    for m in motif_family:
        motif2rep.setdefault(m, m)

    out = refs_dir / "motif_clusters.tsv"
    with open(out, "w") as fh:
        fh.write("motif_id\tfamily\trepresentative_id\n")
        for m in sorted(motif2rep):
            fh.write(f"{m}\t{motif_family.get(m, '?')}\t{motif2rep[m]}\n")

    n_clusters = len(set(motif2rep.values()))
    return ClusterSummary(out, len(motif2rep), n_clusters, len(buckets))


def load_clusters(refs_dir):
    """Return {motif_id: representative_id}; empty dict if not clustered yet."""
    p = Path(refs_dir) / "motif_clusters.tsv"
    if not p.exists():
        return {}
    out = {}
    with open(p) as f:
        f.readline()
        for ln in f:
            c = ln.rstrip("\n").split("\t")
            if len(c) >= 3:
                out[c[0]] = c[2]
    return out
