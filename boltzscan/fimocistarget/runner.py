"""High-level cisTarget DB builder + TF-level aggregation.

Wraps PyMEMESuiteCisTargetBuilder with:
  - per-target ACGT background frequencies computed from the input FASTA
  - optional restriction to motifs belonging to a given TF set (--tf-fasta)
  - motif → TF rollup via a tf2pwms JSON mapping
"""
import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def compute_acgt_background(fasta_path):
    """Count A/C/G/T in `fasta_path` and return frequencies summing to 1.

    N and ambiguous bases are excluded from both numerator and denominator.
    """
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for record in SeqIO.parse(fasta_path, 'fasta'):
        seq = str(record.seq).upper()
        for nt in 'ACGT':
            counts[nt] += seq.count(nt)
    total = sum(counts.values())
    if total == 0:
        raise ValueError(f"No ACGT bases found in {fasta_path}")
    return {nt: counts[nt] / total for nt in 'ACGT'}


def _motif_filter_from_tf_fasta(tf_fasta, tf2pwms):
    """Read TF IDs from `tf_fasta` headers, then collect the union of motif
    IDs assigned to those TFs in `tf2pwms`. Returns a set of motif IDs."""
    tf_ids = [record.id for record in SeqIO.parse(tf_fasta, 'fasta')]
    motifs = set()
    missing = []
    for tf_id in tf_ids:
        hits = tf2pwms.get(tf_id)
        if hits:
            motifs.update(hits)
        else:
            missing.append(tf_id)
    if missing:
        print(
            f"Warning: {len(missing)}/{len(tf_ids)} TF(s) in {tf_fasta} have no "
            f"tf2pwms entry (first few: {missing[:5]})"
        )
    if not motifs:
        raise ValueError(
            f"No motifs matched any TF in {tf_fasta}; check that IDs in the FASTA "
            f"headers match keys in the tf2pwms JSON."
        )
    return motifs


_TF_AGG_FNS = {
    'max': lambda sub: sub.max(axis=1),
    'sum': lambda sub: sub.sum(axis=1),
    'mean': lambda sub: sub.mean(axis=1),
}


def _aggregate_motifs_to_tfs(score_matrix, tf2pwms, tf_agg):
    if tf_agg not in _TF_AGG_FNS:
        raise ValueError(f"Unknown tf_agg: {tf_agg}")
    fn = _TF_AGG_FNS[tf_agg]
    motif_cols = set(score_matrix.columns)
    tf_columns = {}
    skipped = 0
    for tf_id, motif_ids in tf2pwms.items():
        cols = [m for m in motif_ids if m in motif_cols]
        if not cols:
            skipped += 1
            continue
        tf_columns[tf_id] = fn(score_matrix[cols])
    tf_score_matrix = pd.DataFrame(tf_columns, index=score_matrix.index)
    tf_score_matrix.index.name = 'sequence_name'
    return tf_score_matrix, skipped


def run_cistarg_with_tf_aggregation(
    target_fasta,
    motif_dir,
    tf2pwms_path,
    output_dir,
    motif_agg='max',
    tf_agg='max',
    pvalue_thresh=1e-4,
    max_stored_scores=500000,
    motif_pseudo=0.1,
    top_n=3,
    analyze_hits=True,
    tf_fasta=None,
    reuse=False,
):
    """
    Build a cisTarget DB from `target_fasta` + `motif_dir`, then roll the
    motif-level score matrix up to a TF-level score matrix using
    `tf2pwms_path` ({TF_id: [motif_id, ...]}).

    If `tf_fasta` is given, only motifs that map to TFs in that FASTA's
    headers are scanned (subset build, no need to scan the full library).

    The FIMO background frequencies come from `target_fasta` itself.

    If `reuse=True` and `<output_dir>/score_matrix.feather` exists, the FIMO
    rebuild is skipped and only the TF rollup runs.

    Writes into `output_dir`:
      - (from the builder) fimo_results_raw.csv.gz, hit_statistics.csv,
        score_matrix.{feather,csv.gz}, ranking_db.{feather,csv.gz}
      - tf_score_matrix.{feather,csv.gz}, tf_ranking_db.{feather,csv.gz}
      - background_acgt.json (only on rebuild)
    """
    from boltzscan.fimocistarget.fimo_cistarg import PyMEMESuiteCisTargetBuilder

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    score_path = out_dir / 'score_matrix.feather'

    with open(tf2pwms_path) as f:
        tf2pwms = json.load(f)

    motif_filter = None
    if tf_fasta:
        motif_filter = _motif_filter_from_tf_fasta(tf_fasta, tf2pwms)
        print(
            f"--tf-fasta given: restricting scan to {len(motif_filter)} motif(s) "
            f"({len(tf2pwms)} TFs in tf2pwms_dct)"
        )

    if reuse and score_path.exists():
        print(f"Reusing existing score matrix: {score_path}")
        score_matrix = pd.read_feather(score_path).set_index('sequence_name')
    else:
        custom_bg = compute_acgt_background(target_fasta)
        print(
            f"Background ACGT frequencies (from {target_fasta}): "
            f"A={custom_bg['A']:.4f} C={custom_bg['C']:.4f} "
            f"G={custom_bg['G']:.4f} T={custom_bg['T']:.4f}"
        )
        with open(out_dir / 'background_acgt.json', 'w') as f:
            json.dump(custom_bg, f, indent=2)

        builder = PyMEMESuiteCisTargetBuilder(
            fasta_file=target_fasta,
            motif_dir=motif_dir,
            output_dir=str(out_dir),
            custom_bg=custom_bg,
            pvalue_thresh=pvalue_thresh,
            max_stored_scores=max_stored_scores,
            motif_pseudo=motif_pseudo,
            motif_filter=motif_filter,
        )
        score_matrix, _ = builder.build_database(
            aggregation_method=motif_agg,
            analyze_hits=analyze_hits,
            top_n=top_n,
        )

    tf_score_matrix, skipped_tfs = _aggregate_motifs_to_tfs(score_matrix, tf2pwms, tf_agg)
    tf_ranking_db = tf_score_matrix.rank(axis=0, ascending=False, method='min')

    tf_score_matrix.reset_index().to_feather(out_dir / 'tf_score_matrix.feather')
    tf_ranking_db.reset_index().to_feather(out_dir / 'tf_ranking_db.feather')
    tf_score_matrix.to_csv(out_dir / 'tf_score_matrix.csv.gz', compression='gzip')
    tf_ranking_db.to_csv(out_dir / 'tf_ranking_db.csv.gz', compression='gzip')

    print(
        f"TF aggregation done: {tf_score_matrix.shape[0]} seqs × "
        f"{tf_score_matrix.shape[1]} TFs (skipped {skipped_tfs} TFs with no matching motifs)"
    )
    return tf_score_matrix, tf_ranking_db
