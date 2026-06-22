"""CpG-island scanning over a promoter FASTA.

For CXXC-domain / CpG-binding factors (e.g. DDR1) there is no position-specific
PWM motif — their binding feature is the unmethylated CpG island. This module
finds CpG islands per promoter using the Gardiner-Garden & Frommer (1987)
criteria (sliding window >= `window` bp with GC% >= `min_gc` and observed/
expected CpG >= `min_oe`), merges overlapping passing windows into islands, and
emits a per-promoter table + island BED + a candidate target-gene list. This is
the CpG-binding analogue of a cisTarget motif DB.

Observed/expected CpG for a region of length L = n_CpG / (n_C * n_G / L).
"""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio import SeqIO


@dataclass
class CpgSummary:
    out_dir: Path
    n_promoters: int
    n_with_island: int
    n_islands: int


def find_cpg_islands(seq, window=200, min_gc=0.5, min_oe=0.6):
    """Return list of islands as dicts: {start,end,length,gc,oe,n_cpg} (0-based,
    end-exclusive). Window is slid at step 1; passing windows are merged and the
    merged span is re-tested against the same criteria."""
    s = seq.upper()
    L = len(s)
    if L < window:
        return []
    a = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
    is_c = (a == ord("C")).astype(np.int32)
    is_g = (a == ord("G")).astype(np.int32)
    is_cg = ((a[:-1] == ord("C")) & (a[1:] == ord("G"))).astype(np.int32)  # len L-1

    pc = np.concatenate(([0], np.cumsum(is_c)))
    pg = np.concatenate(([0], np.cumsum(is_g)))
    pcg = np.concatenate(([0], np.cumsum(is_cg)))

    def region_stats(start, end):
        n_c = int(pc[end] - pc[start])
        n_g = int(pg[end] - pg[start])
        n_cg = int(pcg[end - 1] - pcg[start])  # CpG dinucs fully inside [start,end)
        ln = end - start
        gc = (n_c + n_g) / ln
        exp = (n_c * n_g) / ln if (n_c and n_g) else 0.0
        oe = (n_cg / exp) if exp > 0 else 0.0
        return gc, oe, n_cg

    # vectorized window pass/fail over all window starts
    starts = np.arange(0, L - window + 1)
    ends = starts + window
    cc = pc[ends] - pc[starts]
    gg = pg[ends] - pg[starts]
    cg = pcg[ends - 1] - pcg[starts]
    gc = (cc + gg) / window
    with np.errstate(divide="ignore", invalid="ignore"):
        exp = (cc * gg) / window
        oe = np.where(exp > 0, cg / exp, 0.0)
    passing = (gc >= min_gc) & (oe >= min_oe)

    # merge runs of consecutive passing window-starts into islands
    islands = []
    i = 0
    n = len(passing)
    while i < n:
        if not passing[i]:
            i += 1
            continue
        j = i
        while j < n and passing[j]:
            j += 1
        start = int(starts[i])
        end = int(ends[j - 1])            # span covered by the run of windows
        g, o, ncg = region_stats(start, end)
        if (end - start) >= window and g >= min_gc and o >= min_oe:
            islands.append({"start": start, "end": end, "length": end - start,
                            "gc": round(g, 4), "oe": round(o, 4), "n_cpg": ncg})
        i = j
    return islands


def scan_promoters_for_cpg(fasta, out_dir, window=200, min_gc=0.5, min_oe=0.6):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    per = out_dir / "cpg_per_promoter.tsv"
    bed = out_dir / "cpg_islands.bed"
    targets = out_dir / "cpg_target_genes.txt"

    n_prom = n_with = n_isl = 0
    with open(per, "w") as ph, open(bed, "w") as bh, open(targets, "w") as th:
        ph.write("gene_id\tlength\tn_islands\tisland_bp\tmax_oe\tmax_gc\thas_cgi\tislands\n")
        for rec in SeqIO.parse(str(fasta), "fasta"):
            n_prom += 1
            islands = find_cpg_islands(str(rec.seq), window, min_gc, min_oe)
            if islands:
                n_with += 1
                n_isl += len(islands)
                th.write(rec.id + "\n")
                for k, isl in enumerate(islands):
                    bh.write(f"{rec.id}\t{isl['start']}\t{isl['end']}\t"
                             f"{rec.id}_cgi{k}\t{isl['oe']}\t+\n")
            isl_bp = sum(i["length"] for i in islands)
            max_oe = max((i["oe"] for i in islands), default=0.0)
            max_gc = max((i["gc"] for i in islands), default=0.0)
            desc = ";".join(f"{i['start']}-{i['end']}({i['oe']},{i['gc']})" for i in islands)
            ph.write(f"{rec.id}\t{len(rec.seq)}\t{len(islands)}\t{isl_bp}\t"
                     f"{max_oe}\t{max_gc}\t{int(bool(islands))}\t{desc}\n")
    return CpgSummary(out_dir, n_prom, n_with, n_isl)
