"""Stage B: map a species' TFs to reference motifs via DBD %ID (run per species)."""
import json
import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from boltzscan.pwmmap import dbd, align
from boltzscan.pwmmap.refs import load_ref_store, load_ref_index
from boltzscan.pwmmap.thresholds import cutoff_for
from boltzscan.pwmmap.cluster import load_clusters


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
                blastp=None, makeblastdb=None, pfam=None, cpu=8,
                collapse_clusters=False):
    pfam = pfam or dbd.DEFAULT_PFAM
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    txt_dir = out_dir/"txt"; meme_dir = out_dir/"meme"
    txt_dir.mkdir(exist_ok=True); meme_dir.mkdir(exist_ok=True)
    store = load_ref_store(refs_dir)
    index = load_ref_index(refs_dir)            # dbd_seq_id -> row(family, motif_ids,...)
    clusters = load_clusters(refs_dir) if collapse_clusters else {}
    bp, mk = align.resolve_blast_bins(blastp, makeblastdb)

    recs = dbd.extract_dbds(species_fasta, domtbl=domtbl, pfam=pfam, cpu=cpu, work_dir=out_dir)
    sp_dbd_fa = out_dir/"species_dbd.fasta"
    dbd.write_dbd_fasta(recs, sp_dbd_fa)
    n_species_tfs = len({r.tf_id for r in recs})

    hits = _blast(sp_dbd_fa, store, bp, mk, cpu, min_cov, out_dir)

    # per TF: motif (or its cluster representative) -> best (pct_id, ref row)
    tf_motif = defaultdict(dict)
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
            key = clusters.get(m, m)            # collapse to cluster rep when enabled
            cur = tf_motif[h.query_tf].get(key)
            if cur is None or h.pct_id > cur[0]:
                tf_motif[h.query_tf][key] = (h.pct_id, row["ref_id"], row["source"],
                                             row["species"], h.query_pfam)

    # write tf2pwms.json + report (best ref per kept motif) + copy motif files.
    # Keep only motifs that have a scannable MEME PWM in the store — cisBP ships
    # some empty/header-only PWMs (no usable matrix); those must not appear in the
    # output (they would be dead entries for FIMO/cisTarget downstream).
    def _scannable(m):
        return (Path(store.motif_meme_dir)/f"{m}.meme").exists()

    report = [("species_tf", "pfam_acc", "ref_id", "source", "species", "pct_id", "motif_id")]
    tf2pwms = {}
    for tf, motifs in tf_motif.items():
        kept = {m: v for m, v in motifs.items() if _scannable(m)}
        if not kept:
            continue
        tf2pwms[tf] = sorted(kept)
        for m, (pid, ref_id, source, species, pfam_acc) in sorted(
                kept.items(), key=lambda kv: -kv[1][0]):
            report.append((tf, pfam_acc, ref_id, source, species, f"{pid:.3f}", m))
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
