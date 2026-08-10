#!/usr/bin/env python
"""Causal ablation of HMM boundary choice versus sequence cropping context.

The experiment is intentionally separate from production BoltzScan.  It uses
the same seed and Boltz2 inference command for four new conditions:

* HMM Pfam envelope union + 40 or 80 aa
* native protein residues within 5 A of the selected DNA chain + 20 or 80 aa

Existing full-length and HMM+20 predictions provide the two anchors.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon
import yaml

from boltzscan.msa import crop_a3m_file, read_a3m_query
from boltzscan.runlog import CommandLog


HERE = Path(__file__).resolve().parent
BENCHMARK = HERE / "benchmark.py"
NEW_CONDITIONS = ("hmm40", "hmm80", "oracle20", "oracle80")
ALL_CONDITIONS = ("full", "hmm20", *NEW_CONDITIONS)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def write_yaml(path: Path, protein: str, dna: str, msa: Path) -> None:
    payload = {
        "version": 1,
        "sequences": [
            {"protein": {"id": ["A"], "sequence": protein, "msa": str(msa)}},
            {"dna": {"id": ["B"], "sequence": dna}},
            {"dna": {"id": ["C"], "sequence": reverse_complement(dna)}},
        ],
    }
    path.write_text(yaml.safe_dump(payload, sort_keys=False))


def hmm_core(dbd_hits: str) -> tuple[int, int]:
    intervals = []
    for hit in str(dbd_hits).split(";"):
        match = re.fullmatch(r"[^:]+:(\d+)-(\d+)", hit)
        if not match:
            raise ValueError(f"Cannot parse DBD hit: {hit}")
        intervals.append(tuple(map(int, match.groups())))
    return min(start for start, _ in intervals), max(stop for _, stop in intervals)


def heavy_atom_coordinates(residues) -> np.ndarray:
    coords = [
        atom.coord
        for residue in residues
        for atom in residue.get_atoms()
        if str(getattr(atom, "element", "")).upper() != "H"
    ]
    return np.asarray(coords, dtype=float)


def native_contacts(row, structure, cutoff: float) -> list[int]:
    protein = structure.child_dict.get(str(row.native_protein_chain))
    dna = structure.child_dict.get(str(row.native_dna_chain))
    if protein is None or dna is None:
        raise ValueError(
            f"Missing native chain for {row.job_id}: "
            f"protein={row.native_protein_chain}, DNA={row.native_dna_chain}"
        )
    native_sequence = str(getattr(protein, "sequence", ""))
    input_sequence = str(row.protein_sequence_full)
    if native_sequence != input_sequence:
        raise ValueError(
            f"Native/input sequence mismatch for {row.job_id}: "
            f"{len(native_sequence)} != {len(input_sequence)}"
        )
    dna_coords = heavy_atom_coordinates(dna.get_residues())
    if not len(dna_coords):
        return []
    tree = cKDTree(dna_coords)
    contacts = []
    for index, residue in enumerate(protein.get_residues(), start=1):
        coords = heavy_atom_coordinates([residue])
        if len(coords) and float(tree.query(coords, k=1)[0].min()) <= cutoff:
            contacts.append(index)
    return contacts


def interval(core_start: int, core_stop: int, flank: int, length: int):
    return max(1, core_start - flank), min(length, core_stop + flank)


def load_manifests(full_run: Path, hmm20_run: Path):
    full = pd.read_csv(
        full_run / "benchmark_manifest.csv",
        dtype={"native_protein_chain": str, "native_dna_chain": str},
    )
    hmm = pd.read_csv(
        hmm20_run / "benchmark_manifest.csv",
        dtype={"native_protein_chain": str, "native_dna_chain": str},
    )
    shortened = hmm[hmm["protein_length"] < hmm["full_protein_length"]].copy()
    full_columns = [
        "job_id", "protein_sequence", "protein_length", "msa_path",
    ]
    merged = shortened.drop(columns=[
        column for column in ("full_protein_sequence", "full_protein_length")
        if column in shortened
    ]).merge(
        full[full_columns],
        on="job_id",
        how="inner",
        suffixes=("_hmm20", "_full"),
        validate="one_to_one",
    )
    return full, hmm, merged


def prepare(args) -> None:
    from DockQ.DockQ import load_PDB

    full_run = Path(args.full_run)
    hmm20_run = Path(args.hmm20_run)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    with CommandLog(out / "prepare.log", sys.argv, args):
        full_manifest, _, candidates = load_manifests(full_run, hmm20_run)
        cache = {}
        audits, excluded = [], []
        for row in candidates.itertuples(index=False):
            native_path = str(row.native_structure)
            if native_path not in cache:
                cache[native_path] = load_PDB(native_path)
            contacts = native_contacts(row, cache[native_path], args.contact_cutoff)
            if not contacts:
                excluded.append({
                    "job_id": row.job_id,
                    "reason": (
                        f"no protein heavy atom within {args.contact_cutoff:g} A "
                        "of selected native DNA chain"
                    ),
                })
                continue
            hmm_start, hmm_stop = hmm_core(row.dbd_hits)
            audit = row._asdict()
            audit.update({
                "hmm_core_start": hmm_start,
                "hmm_core_stop": hmm_stop,
                "oracle_contact_start": min(contacts),
                "oracle_contact_stop": max(contacts),
                "oracle_contact_residue_count": len(contacts),
                "oracle_contact_positions": ";".join(map(str, contacts)),
                "contact_cutoff_angstrom": args.contact_cutoff,
            })
            audits.append(audit)

        audit = pd.DataFrame(audits)
        excluded = pd.DataFrame(excluded, columns=["job_id", "reason"])
        audit.to_csv(out / "coordinate_audit.csv", index=False)
        excluded.to_csv(out / "excluded_oracle.csv", index=False)
        if audit.empty:
            raise ValueError("No oracle-contact inputs were retained")

        specs = {
            "hmm40": ("hmm_core_start", "hmm_core_stop", 40),
            "hmm80": ("hmm_core_start", "hmm_core_stop", 80),
            "oracle20": ("oracle_contact_start", "oracle_contact_stop", 20),
            "oracle80": ("oracle_contact_start", "oracle_contact_stop", 80),
        }
        condition_rows = []
        for condition, (start_column, stop_column, flank) in specs.items():
            run = out / condition
            input_dir = run / "inputs" / "boltz2"
            input_dir.mkdir(parents=True, exist_ok=True)
            for stale in input_dir.glob("*.yaml"):
                stale.unlink()
            retained = []
            for row in audit.itertuples(index=False):
                protein = str(row.protein_sequence_full)
                crop_start, crop_stop = interval(
                    int(getattr(row, start_column)),
                    int(getattr(row, stop_column)),
                    flank,
                    len(protein),
                )
                crop_sequence = protein[crop_start - 1:crop_stop]
                crop_msa = (run / "msa" / row.job_id / "crop_0.a3m").resolve()
                crop_a3m_file(
                    Path(row.msa_path_full),
                    crop_msa,
                    interval=(crop_start - 1, crop_stop),
                )
                if read_a3m_query(crop_msa) != crop_sequence:
                    raise ValueError(f"Protein/A3M query mismatch for {condition}/{row.job_id}")
                record = {
                    "dataset_row": row.dataset_row,
                    "job_id": row.job_id,
                    "complex_id": row.complex_id,
                    "native_structure": row.native_structure,
                    "native_protein_chain": row.native_protein_chain,
                    "native_dna_chain": row.native_dna_chain,
                    "protein_sequence": crop_sequence,
                    "dna_sequence": row.dna_sequence,
                    "protein_length": len(crop_sequence),
                    "dna_length": row.dna_length,
                    "total_residues": len(crop_sequence) + 2 * int(row.dna_length),
                    "msa_path": str(crop_msa),
                    "full_protein_length": len(protein),
                    "crop_start_1based_inclusive": crop_start,
                    "crop_stop_1based_inclusive": crop_stop,
                    "crop_definition": condition.rstrip("0123456789"),
                    "crop_flank_aa": flank,
                }
                retained.append(record)
                condition_rows.append({
                    "condition": condition,
                    "job_id": row.job_id,
                    "full_length": len(protein),
                    "crop_length": len(crop_sequence),
                    "crop_start": crop_start,
                    "crop_stop": crop_stop,
                    "actually_shortened": len(crop_sequence) < len(protein),
                })
                write_yaml(
                    input_dir / f"{row.job_id}.yaml",
                    crop_sequence,
                    str(row.dna_sequence),
                    crop_msa,
                )

            manifest = pd.DataFrame(retained)
            manifest.to_csv(run / "benchmark_manifest.csv", index=False)
            config = {
                "schema_version": 1,
                "benchmark_mode": "crop_causality",
                "condition": condition,
                "source_full_run": str(full_run.resolve()),
                "source_hmm20_run": str(hmm20_run.resolve()),
                "configs": ["boltz2"],
                "selected_interfaces": len(manifest),
                "seed": args.seed,
                "contact_cutoff_angstrom": args.contact_cutoff,
                "historical_predictions_used": False,
            }
            (run / "benchmark_config.json").write_text(json.dumps(config, indent=2) + "\n")

        condition_audit = pd.DataFrame(condition_rows)
        condition_audit.to_csv(out / "condition_intervals.csv", index=False)
        summary = condition_audit.groupby("condition", sort=False).agg(
            n=("job_id", "size"),
            actually_shortened=("actually_shortened", "sum"),
            mean_full_length=("full_length", "mean"),
            mean_crop_length=("crop_length", "mean"),
            median_crop_length=("crop_length", "median"),
        ).reset_index()
        summary.to_csv(out / "input_summary.csv", index=False)
        metadata = {
            "generated": datetime.now().astimezone().isoformat(timespec="seconds"),
            "candidate_hmm20_shortened": len(candidates),
            "retained_native_contact": len(audit),
            "excluded_no_selected_chain_contact": len(excluded),
            "contact_definition": (
                f"protein heavy atom within {args.contact_cutoff:g} A of any heavy atom "
                "in the selected native DNA chain"
            ),
            "conditions": list(ALL_CONDITIONS),
            "seed": args.seed,
        }
        (out / "experiment.json").write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"Prepared {len(audit)} causal-ablation jobs x {len(specs)} new conditions -> {out}")


def run_benchmark_stage(args, stage: str) -> None:
    out = Path(args.out)
    with CommandLog(out / f"{stage}.log", sys.argv, args):
        for condition in NEW_CONDITIONS:
            run = out / condition
            command = [
                sys.executable,
                str(BENCHMARK),
                stage,
                "--run",
                str(run),
                "--configs",
                "boltz2",
            ]
            if stage == "predict":
                command.extend(["--seed", str(args.seed)])
            else:
                command.extend(["--processes", str(args.processes)])
            print("[causality] " + " ".join(command))
            subprocess.run(command, check=True)


def paired_test(left, right):
    delta = right - left
    nonzero = delta != 0
    pvalue = np.nan
    if nonzero.any():
        pvalue = wilcoxon(right[nonzero], left[nonzero], zero_method="wilcox").pvalue
    return {
        "n": len(delta),
        "left_mean": left.mean(),
        "right_mean": right.mean(),
        "mean_delta_right_minus_left": delta.mean(),
        "median_delta_right_minus_left": delta.median(),
        "right_wins": int((delta > 0).sum()),
        "ties": int((delta == 0).sum()),
        "right_losses": int((delta < 0).sum()),
        "wilcoxon_p": pvalue,
    }


def markdown(frame: pd.DataFrame) -> str:
    def render(value):
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{value:.4g}"
        return str(value)
    rows = [
        "| " + " | ".join(frame.columns) + " |",
        "| " + " | ".join("---" for _ in frame.columns) + " |",
    ]
    rows.extend(
        "| " + " | ".join(render(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    )
    return "\n".join(rows)


def analyze(args) -> None:
    out = Path(args.out)
    full_run = Path(args.full_run)
    hmm20_run = Path(args.hmm20_run)
    with CommandLog(out / "analyze.log", sys.argv, args):
        metrics = {}
        sources = {
            "full": full_run / "metrics" / "metrics_per_case.csv",
            "hmm20": hmm20_run / "metrics" / "metrics_per_case.csv",
            **{
                condition: out / condition / "metrics" / "metrics_per_case.csv"
                for condition in NEW_CONDITIONS
            },
        }
        retained = set(pd.read_csv(out / "coordinate_audit.csv")["job_id"])
        for condition, path in sources.items():
            frame = pd.read_csv(path)
            frame = frame[
                frame["config"].eq("boltz2") & frame["job_id"].isin(retained)
            ][["job_id", "dockq_interface"]].copy()
            metrics[condition] = frame.rename(
                columns={"dockq_interface": f"dockq_{condition}"}
            )

        paired = metrics["full"]
        for condition in ALL_CONDITIONS[1:]:
            paired = paired.merge(metrics[condition], on="job_id", how="inner")
        score_columns = [f"dockq_{condition}" for condition in ALL_CONDITIONS]
        paired = paired.dropna(subset=score_columns).copy()
        paired.to_csv(out / "paired_dockq_per_case.csv", index=False)

        condition_rows = []
        for condition in ALL_CONDITIONS:
            values = paired[f"dockq_{condition}"]
            condition_rows.append({
                "condition": condition,
                "n": len(values),
                "dockq_mean": values.mean(),
                "dockq_median": values.median(),
                "acceptable_or_better": (values >= 0.23).mean(),
                "medium_or_better": (values >= 0.49).mean(),
                "high": (values >= 0.80).mean(),
            })
        condition_summary = pd.DataFrame(condition_rows)
        condition_summary.to_csv(out / "condition_dockq_summary.csv", index=False)

        comparisons = [
            ("full", "hmm20", "total observed HMM+20 crop effect"),
            ("full", "oracle20", "native-contact crop effect at flank 20"),
            ("hmm20", "oracle20", "boundary-definition effect at flank 20"),
            ("hmm20", "hmm40", "HMM flank recovery: 20 to 40"),
            ("hmm40", "hmm80", "HMM flank recovery: 40 to 80"),
            ("hmm20", "hmm80", "HMM flank recovery: 20 to 80"),
            ("oracle20", "oracle80", "oracle flank recovery: 20 to 80"),
            ("oracle80", "full", "residual context recovery: oracle80 to full"),
            ("hmm80", "oracle80", "boundary-definition effect at flank 80"),
        ]
        comparison_rows = []
        for left, right, interpretation in comparisons:
            result = paired_test(paired[f"dockq_{left}"], paired[f"dockq_{right}"])
            result.update({
                "left": left,
                "right": right,
                "interpretation": interpretation,
            })
            comparison_rows.append(result)
        comparison = pd.DataFrame(comparison_rows)[[
            "left", "right", "interpretation", "n", "left_mean", "right_mean",
            "mean_delta_right_minus_left", "median_delta_right_minus_left",
            "right_wins", "ties", "right_losses", "wilcoxon_p",
        ]]
        comparison.to_csv(out / "paired_comparisons.csv", index=False)

        input_summary = pd.read_csv(out / "input_summary.csv")
        indexed = comparison.set_index(["left", "right"])
        total = indexed.loc[("full", "hmm20")]
        oracle_crop20 = indexed.loc[("full", "oracle20")]
        boundary20 = indexed.loc[("hmm20", "oracle20")]
        hmm_flank = indexed.loc[("hmm20", "hmm80")]
        oracle_flank = indexed.loc[("oracle20", "oracle80")]
        residual = indexed.loc[("oracle80", "full")]
        report = f"""# FoldBench crop-causality ablation

- Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}
- Predictor: Boltz2, seed {args.seed}
- Candidate set: proteins actually shortened by the original HMM+20 rule
- Native oracle: protein heavy atoms within 5 A of the selected native DNA chain
- Strict paired scoreable jobs across all six conditions: {len(paired)}
- Full-length fallback inside a crop condition: reported, never hidden

## Input lengths

{markdown(input_summary.round(4))}

## Strict paired DockQ

{markdown(condition_summary.round(4))}

## Pre-specified paired contrasts

{markdown(comparison)}

## Causal interpretation

The total HMM+20 crop change was {total['mean_delta_right_minus_left']:+.4f}
DockQ relative to full length. Replacing HMM boundaries with native-contact
boundaries at the same 20-aa flank changed DockQ by
{boundary20['mean_delta_right_minus_left']:+.4f}. Cropping around the native
contact boundary at flank 20 changed DockQ by
{oracle_crop20['mean_delta_right_minus_left']:+.4f} relative to full length.
Expanding only the HMM flank
from 20 to 80 aa recovered {hmm_flank['mean_delta_right_minus_left']:+.4f},
whereas expanding the native-contact flank recovered
{oracle_flank['mean_delta_right_minus_left']:+.4f}. Moving from oracle+80 back
to full length changed DockQ by {residual['mean_delta_right_minus_left']:+.4f}.

These effects are not assumed to be additive: boundary choice changes both
the retained residues and crop length. The oracle uses native structural
information and is therefore a causal diagnostic, not a deployable predictor.
"""
        (out / "crop_causality_report.md").write_text(report)
        print(f"Analyzed {len(paired)} strict paired jobs -> {out}")


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest="stage", required=True)
    for stage in ("prepare", "predict", "score", "analyze"):
        child = sub.add_parser(stage)
        child.add_argument("--full-run", required=True)
        child.add_argument("--hmm20-run", required=True)
        child.add_argument("--out", required=True)
        child.add_argument("--seed", type=int, default=42)
        child.add_argument("--contact-cutoff", type=float, default=5.0)
        child.add_argument("--processes", type=int, default=12)
    return root


def main() -> None:
    args = parser().parse_args()
    if args.stage == "prepare":
        prepare(args)
    elif args.stage in {"predict", "score"}:
        run_benchmark_stage(args, args.stage)
    else:
        analyze(args)


if __name__ == "__main__":
    main()
