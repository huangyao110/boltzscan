"""Post-FIMO pipeline.

Turn the raw FIMO hits produced by `PyMEMESuiteCisTargetBuilder` into the
per-row inputs that Boltz expects:
  - one FASTA per (TF, motif, promoter, hit) with the TF protein chain (+ MSA),
    the cropped DNA motif, and its reverse complement;
  - distributed across N subdirectories for parallel scheduling.

The function mirrors the notebook flow in `boltzscan.ipynb`.
"""
import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def _read_fasta_to_str(path):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, 'fasta')}


def _build_msa_dct(msa_dir, a3m_name='0.a3m'):
    """Walk `msa_dir` and return {subdir_name: <subdir>/<a3m_name>} for every
    subdirectory that actually contains the a3m file."""
    msa_dir = Path(msa_dir)
    dct = {}
    if not msa_dir.exists():
        return dct
    for sub in msa_dir.iterdir():
        if sub.is_dir():
            a3m = sub / a3m_name
            if a3m.exists():
                dct[sub.name] = str(a3m)
    return dct


def _invert_tf2pwms(tf2pwms):
    """{tf: [pwm, ...]} -> {pwm: [tf, ...]}"""
    inv = {}
    for tf, pwms in tf2pwms.items():
        for pwm in pwms:
            inv.setdefault(pwm, []).append(tf)
    return inv


def build_boltz_input_from_fimo(
    fimo_csv,
    tf_fasta,
    tf2pwms_path,
    promoter_fasta,
    msa_dir,
    output_dir,
    pvalue_thresh=1e-5,
    overlap_thresh=0.0,
    extend_range=5,
    layout='auto',
    auto_threshold=20000,
    num_subdirs=100,
    max_workers=8,
):
    """
    1. Load `fimo_csv`, filter by p-value <= `pvalue_thresh`.
    2. Annotate each row with promoter_seq (from `promoter_fasta`) and tf_name
       (from inverted `tf2pwms_path`), exploding when one motif maps to many TFs.
    3. Drop rows whose TF isn't in `tf_fasta` (keep only buildable TFs).
    4. Dedupe via `filter_fimo_by_seq_and_overlap` (sequence-content +
       physical overlap, per (sequence_name, tf_name) group).
    5. Attach tf_seq (from `tf_fasta`) and msa_path (resolved against
       `msa_dir` using the `.`->`_` normalised TF id as subdirectory name).
    6. Write per-row Boltz FASTAs via `to_boltz_fasta`, then split them into
       `num_subdirs` buckets via `distribute_files`.

    Writes into `output_dir`:
        fimo_results_filtered_<p>.csv.gz
        fimo_results_filtered_with_seq_<p>.csv.gz
        boltzscan_input/BOLTZ_INPUT/<boltz_name>.fasta   (then redistributed)
        boltzscan_input/subdir_NNN/*.fasta
    """
    from boltzscan.fimocistarget.cistarg_impl import filter_fimo_by_seq_and_overlap
    from boltzscan.utils.io_utils import to_boltz_fasta, distribute_files

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    save_boltz_dir = output_dir / 'boltzscan_input'
    save_boltz_dir.mkdir(parents=True, exist_ok=True)

    hit = pd.read_csv(fimo_csv)
    print(f'Total hits: {hit.shape[0]}')
    fimo_df = hit[hit['pvalue'] <= pvalue_thresh].copy()
    print(f'Hits after filtering with pvalue <= {pvalue_thresh}: {fimo_df.shape[0]}')

    promoter_seqs = _read_fasta_to_str(promoter_fasta)
    pep_seqs = _read_fasta_to_str(tf_fasta)
    with open(tf2pwms_path) as f:
        tf2pwms = json.load(f)
    pwm2tf = _invert_tf2pwms(tf2pwms)

    fimo_df['promoter_seq'] = fimo_df['sequence_name'].map(promoter_seqs)
    fimo_df = fimo_df.dropna(subset=['promoter_seq'])
    fimo_df['tf_name'] = fimo_df['motif_id'].map(pwm2tf)
    fimo_df = fimo_df.explode('tf_name').dropna(subset=['tf_name'])
    fimo_df = fimo_df[fimo_df['tf_name'].isin(pep_seqs)]
    print(f'Hits after restricting to TFs in {tf_fasta}: {fimo_df.shape[0]}')

    fimo_df = filter_fimo_by_seq_and_overlap(
        fimo_df, overlap_threshold=overlap_thresh, extend_range=extend_range,
    )
    print(f'Hits after removing redundancy: {fimo_df.shape[0]}')

    filtered_path = output_dir / f'fimo_results_filtered_{pvalue_thresh}.csv.gz'
    fimo_df.to_csv(filtered_path, index=False)

    fimo_df['boltz_name'] = (
        fimo_df['tf_name'] + '-' + fimo_df['motif_id'] + '-'
        + fimo_df['sequence_name'] + '-' + fimo_df['start'].astype(str)
        + '-' + fimo_df['stop'].astype(str)
    )
    fimo_df['tf_name_new'] = fimo_df['tf_name'].str.replace('.', '_', regex=False)
    tf_msa_dct = _build_msa_dct(msa_dir)
    fimo_df['msa_path'] = fimo_df['tf_name_new'].map(tf_msa_dct)
    fimo_df['tf_seq'] = fimo_df['tf_name'].map(pep_seqs)

    with_seq_path = output_dir / f'fimo_results_filtered_with_seq_{pvalue_thresh}.csv.gz'
    fimo_df.to_csv(with_seq_path, index=False)

    boltz_input_dir = save_boltz_dir / 'BOLTZ_INPUT'
    print(f'Writing Boltz input FASTAs to {boltz_input_dir}')
    to_boltz_fasta(fimo_df, save_boltz_dir, max_workers=max_workers)

    if layout == 'auto':
        n_files = sum(1 for _ in boltz_input_dir.iterdir())
        chosen = 'sub' if n_files > auto_threshold else 'all'
        print(
            f'Layout=auto: {n_files} files vs threshold {auto_threshold} '
            f'→ keeping {"sub" if chosen == "sub" else "all"} layout'
        )
        layout = chosen

    if layout == 'sub':
        distribute_files(boltz_input_dir, save_boltz_dir, num_subdirs=num_subdirs)
    elif layout == 'all':
        print(f'Layout=all: leaving FASTAs flat in {boltz_input_dir}')
    else:
        raise ValueError(f'Unknown layout: {layout!r} (expected sub|all|auto)')

    return fimo_df
