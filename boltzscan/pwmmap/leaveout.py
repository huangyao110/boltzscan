"""Create leakage-controlled subsets of a PWM reference store.

The reference store contains the DBD sequences used for homology transfer and
the PWM files that may subsequently be copied into a species-specific mapping.
For a leave-species-out benchmark, both pieces must be filtered *before*
BLAST/internal motif clustering; deleting PWM IDs after mapping would leave target-species
information in the reference alignment and in cluster-representative choices.
"""

from __future__ import annotations

import csv
import json
import re
import shutil
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO


@dataclass(frozen=True)
class ReferenceSubsetSummary:
    """Auditable counts for a filtered reference store."""

    out_dir: Path
    excluded_species: tuple[str, ...]
    n_input_dbd_rows: int
    n_retained_dbd_rows: int
    n_excluded_dbd_rows: int
    n_retained_motifs: int
    n_excluded_motifs: int
    n_copied_txt_pwms: int
    n_copied_meme_pwms: int


def _normalise_species(value: str) -> str:
    """Normalise common database spelling differences without fuzzy matching."""

    value = value.replace("_", " ").casefold()
    return re.sub(r"\s+", " ", value).strip()


def _has_excluded_species(value: str, excluded: set[str]) -> bool:
    # JASPAR can store more than one species in one cell.  Exact matching after
    # normalisation intentionally avoids excluding a taxon just because its
    # name happens to contain another taxon's name.
    return any(_normalise_species(part) in excluded for part in value.split(";"))


def _motif_ids(row: dict[str, str]) -> set[str]:
    return {motif for motif in row["motif_ids"].split(";") if motif}


def create_reference_subset(
    source_refs: str | Path,
    out_dir: str | Path,
    *,
    exclude_species: list[str] | tuple[str, ...],
) -> ReferenceSubsetSummary:
    """Write an isolated, species-excluded PWM reference store.

    The source reference store is never modified.  Only the retained DBD FASTA
    records and PWM files are copied to ``out_dir``.  ``motif_clusters.tsv`` is
    deliberately not copied: callers must recluster the subset before using
    ``map-pwm --collapse-clusters``.

    A PWM appearing in both an excluded and retained reference is rejected.
    Treating that situation as an error is conservative: otherwise a target
    species' PWM might silently remain available through a shared motif ID.
    """

    source_refs = Path(source_refs)
    out_dir = Path(out_dir)
    if not exclude_species:
        raise ValueError("At least one --exclude-species value is required")
    if not (source_refs / "ref_index.tsv").exists():
        raise FileNotFoundError(f"Missing reference index: {source_refs / 'ref_index.tsv'}")
    if not (source_refs / "ref_dbd.fasta").exists():
        raise FileNotFoundError(f"Missing reference DBD FASTA: {source_refs / 'ref_dbd.fasta'}")
    if out_dir.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing subset store: {out_dir}. "
            "Choose a new output directory so its provenance remains immutable."
        )

    excluded_normalised = {_normalise_species(s) for s in exclude_species}
    with open(source_refs / "ref_index.tsv", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames or not {"species", "dbd_seq_id", "motif_ids"}.issubset(fieldnames):
            raise ValueError("ref_index.tsv is missing required columns")
        all_rows = list(reader)

    excluded_rows = [
        row for row in all_rows
        if _has_excluded_species(row["species"], excluded_normalised)
    ]
    retained_rows = [row for row in all_rows if row not in excluded_rows]
    if not retained_rows:
        raise ValueError("Species filter removed every reference DBD row")

    excluded_motifs = set().union(*(_motif_ids(row) for row in excluded_rows)) if excluded_rows else set()
    retained_motifs = set().union(*(_motif_ids(row) for row in retained_rows))
    overlap = sorted(excluded_motifs & retained_motifs)
    if overlap:
        preview = ", ".join(overlap[:10])
        raise ValueError(
            "A PWM ID occurs in both excluded and retained references; refusing "
            f"a potentially leaky subset (first IDs: {preview})"
        )

    retained_dbd_ids = {row["dbd_seq_id"] for row in retained_rows}
    out_dir.mkdir(parents=True)
    try:
        with open(out_dir / "ref_index.tsv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(retained_rows)

        n_written_dbd = 0
        with open(out_dir / "ref_dbd.fasta", "w") as handle:
            for record in SeqIO.parse(str(source_refs / "ref_dbd.fasta"), "fasta"):
                if record.id in retained_dbd_ids:
                    SeqIO.write(record, handle, "fasta")
                    n_written_dbd += 1
        if n_written_dbd != len(retained_dbd_ids):
            raise ValueError(
                "Reference index/FASTA mismatch: expected "
                f"{len(retained_dbd_ids)} retained DBD records, wrote {n_written_dbd}"
            )

        copied = {"txt": 0, "meme": 0}
        for kind in ("txt", "meme"):
            source_dir = source_refs / "motif_store" / kind
            target_dir = out_dir / "motif_store" / kind
            target_dir.mkdir(parents=True)
            suffix = f".{kind}"
            for motif_id in sorted(retained_motifs):
                source_file = source_dir / f"{motif_id}{suffix}"
                if source_file.exists():
                    shutil.copy2(source_file, target_dir / source_file.name)
                    copied[kind] += 1

        summary = ReferenceSubsetSummary(
            out_dir=out_dir,
            excluded_species=tuple(exclude_species),
            n_input_dbd_rows=len(all_rows),
            n_retained_dbd_rows=len(retained_rows),
            n_excluded_dbd_rows=len(excluded_rows),
            n_retained_motifs=len(retained_motifs),
            n_excluded_motifs=len(excluded_motifs),
            n_copied_txt_pwms=copied["txt"],
            n_copied_meme_pwms=copied["meme"],
        )
        manifest = {
            **asdict(summary),
            "out_dir": str(out_dir),
            "source_refs": str(source_refs.resolve()),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "cluster_policy": (
                "No motif_clusters.tsv copied. The unified run clusters this "
                "filtered subset internally before representative PWM mapping."
            ),
            "leakage_check": {
                "excluded_retained_motif_id_overlap": 0,
                "status": "passed",
            },
        }
        (out_dir / "subset_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        return summary
    except Exception:
        # A failed build must not look like a valid immutable reference store.
        shutil.rmtree(out_dir, ignore_errors=True)
        raise
