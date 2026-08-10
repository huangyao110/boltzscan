from dataclasses import dataclass
import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO


@dataclass
class CandidateTfFastaSummary:
    output: Path
    hit_motifs: int
    matched_motifs: int
    candidate_tfs: int
    written_tfs: int
    missing_tfs: list[str]
    unmatched_motifs: list[str]


def build_candidate_tf_fasta(
    hit_statistics,
    tf2pwms_path,
    protein_fasta,
    output,
    motif_col='motif_id',
):
    """Build a candidate TF protein FASTA from motif hits and TF->PWM mapping."""
    hit_statistics = Path(hit_statistics)
    tf2pwms_path = Path(tf2pwms_path)
    protein_fasta = Path(protein_fasta)
    output = Path(output)

    hit_df = pd.read_csv(hit_statistics)
    if motif_col not in hit_df.columns:
        raise ValueError(f"motif column not found in hit statistics: {motif_col}")

    hit_motifs = set(hit_df[motif_col].dropna().astype(str))
    with open(tf2pwms_path) as f:
        tf2pwms = json.load(f)

    motif_to_tfs = _invert_tf2pwms(tf2pwms)
    matched_motifs = hit_motifs.intersection(motif_to_tfs)
    unmatched_motifs = sorted(hit_motifs.difference(motif_to_tfs))

    candidate_tfs = sorted(
        {
            tf
            for motif in matched_motifs
            for tf in motif_to_tfs[motif]
        }
    )

    protein_records = _read_fasta_records(protein_fasta)
    output_records = []
    missing_tfs = []
    for tf in candidate_tfs:
        record = protein_records.get(tf)
        if record is None:
            missing_tfs.append(tf)
            continue
        output_records.append(record)

    if not output_records:
        raise ValueError(
            "No candidate TF protein sequences were found. "
            "Check that TF IDs in tf2pwms match protein FASTA record IDs."
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(output_records, output, 'fasta')

    return CandidateTfFastaSummary(
        output=output,
        hit_motifs=len(hit_motifs),
        matched_motifs=len(matched_motifs),
        candidate_tfs=len(candidate_tfs),
        written_tfs=len(output_records),
        missing_tfs=missing_tfs,
        unmatched_motifs=unmatched_motifs,
    )


def _invert_tf2pwms(tf2pwms):
    motif_to_tfs = {}
    for tf, motifs in tf2pwms.items():
        for motif in motifs:
            motif_to_tfs.setdefault(str(motif), []).append(str(tf))
    return motif_to_tfs


def _read_fasta_records(fasta_path):
    records = {}
    for record in SeqIO.parse(fasta_path, 'fasta'):
        records.setdefault(record.id, record)
    return records
