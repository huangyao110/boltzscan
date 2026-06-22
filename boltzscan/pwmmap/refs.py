"""Stage A: build the reusable cisBP+JASPAR reference DBD store (run once)."""
import json
from dataclasses import dataclass
from pathlib import Path

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
    return RefStore(root, root / "ref_dbd.fasta", root / "ref_index.tsv",
                    root / "motif_store" / "txt", root / "motif_store" / "meme")


def load_ref_index(refs_dir):
    """Load ref_index.tsv keyed by dbd_seq_id (blast subject id); one TF may have multiple DBD rows."""
    out = {}
    with open(Path(refs_dir) / "ref_index.tsv") as f:
        header = f.readline().rstrip("\n").split("\t")
        for ln in f:
            vals = ln.rstrip("\n").split("\t")
            row = dict(zip(header, vals))
            out[row["dbd_seq_id"]] = row
    return out


def build_reference_db(refs_dir="data/pwms/_refs", pfam=None,
                       cpu=8, refresh=False, include_cisbp=True, include_jaspar=True):
    pfam = pfam or dbd.DEFAULT_PFAM
    root = Path(refs_dir)

    # idempotency — skip rebuild if both artifacts exist and refresh not requested
    if not refresh and (root / "ref_dbd.fasta").exists() and (root / "ref_index.tsv").exists():
        print(f"[refs] reuse existing store at {root}")
        return load_ref_store(root)

    if not Path(pfam).exists():
        raise SystemExit(f"Pfam HMM not found: {pfam} — build/point --pfam to a hmmpress-ed Pfam-A.hmm")

    root.mkdir(parents=True, exist_ok=True)
    txt_dir = root / "motif_store" / "txt"
    meme_dir = root / "motif_store" / "meme"
    txt_dir.mkdir(parents=True, exist_ok=True)
    meme_dir.mkdir(parents=True, exist_ok=True)

    refs = []
    if include_cisbp:
        cis_refs, pwms_dir = _gather_cisbp(root / "cisbp", refresh)
        refs += cis_refs
        # copy needed cisBP motif txt -> motif_store, convert to meme
        if pwms_dir:
            pwms_path = Path(pwms_dir)
            wanted = {m for r in cis_refs for m in r.motif_ids}
            for m in wanted:
                src = pwms_path / f"{m}.txt"
                if src.exists():
                    try:
                        pwmio.copy_cisbp_pwm(src, m, txt_dir, meme_dir)
                    except Exception as e:  # degenerate cisBP txt
                        print(f"[refs] skip cisBP {m}: {e}")
    if include_jaspar:
        refs += _gather_jaspar(root / "jaspar", txt_dir, meme_dir, refresh)

    seqs = _resolve_seqs(refs, root / "uniprot_cache.json")

    # write reference protein fasta (only refs with a sequence)
    prot_fa = root / "ref_proteins.fasta"
    by_id = {}
    with open(prot_fa, "w") as fh:
        for ref in refs:
            s = seqs.get(ref.ref_id)
            if not s:
                continue
            fh.write(f">{ref.ref_id}\n{s}\n")
            by_id[ref.ref_id] = ref

    # extract reference DBDs (hmmsearch over ref proteins), write ref_dbd.fasta + index
    # Fix 3: capture written ids from write_dbd_fasta to keep fasta and index in lockstep
    recs = dbd.extract_dbds(prot_fa, domtbl=None, pfam=pfam, cpu=cpu, work_dir=root)
    written_ids = dbd.write_dbd_fasta(recs, root / "ref_dbd.fasta")
    with open(root / "ref_index.tsv", "w") as fh:
        fh.write("ref_id\tsource\tspecies\tfamily\tpfam_acc\tdbd_seq_id\tmotif_ids\n")
        for r, sid in zip(recs, written_ids):
            ref = by_id.get(r.tf_id)
            if not ref:
                continue
            fh.write(f"{ref.ref_id}\t{ref.source}\t{ref.species}\t{ref.family}\t"
                     f"{r.pfam_acc}\t{sid}\t{';'.join(ref.motif_ids)}\n")

    (root / "build_manifest.json").write_text(json.dumps({
        "n_refs": len(refs), "n_with_seq": len(by_id), "n_dbd": len(recs),
        "include_cisbp": include_cisbp, "include_jaspar": include_jaspar,
    }, indent=2))
    print(f"[refs] {len(refs)} refs, {len(by_id)} with seq, {len(recs)} DBDs -> {root}")
    return load_ref_store(root)
