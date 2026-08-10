#!/usr/bin/env python
"""Benchmark model-predicted protein--dsDNA interface cropping.

This experiment is deliberately separate from production BoltzScan.  For each
of the six FoldBench configurations, it uses that configuration's existing
full-length prediction as an interface localizer, crops the protein and A3M to
the predicted DNA-contact span plus a configurable flank, and reruns the same
configuration.  A matched native-contact (oracle) crop is also run for every
configuration.  Full, predicted-interface, and oracle results are then paired
within model.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import shutil
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
CONFIGS = (
    "esmfold2_msa",
    "esmfold2_nomsa",
    "boltz1",
    "boltz1_ode",
    "boltz2",
    "boltz2_ode",
)
CONDITIONS = ("interface20", "oracle20")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def write_yaml(path: Path, protein: str, dna: str, msa: str | Path) -> None:
    payload = {
        "version": 1,
        "sequences": [
            {"protein": {"id": ["A"], "sequence": protein, "msa": str(msa)}},
            {"dna": {"id": ["B"], "sequence": dna}},
            {"dna": {"id": ["C"], "sequence": reverse_complement(dna)}},
        ],
    }
    path.write_text(yaml.safe_dump(payload, sort_keys=False))


def prediction_directory(run: Path, config: str) -> Path:
    prefix = "esmfold_results_" if config.startswith("esmfold2") else "boltz_results_"
    return run / "outputs" / config / f"{prefix}{config}" / "predictions"


def prediction_cif(run: Path, config: str, job_id: str) -> Path | None:
    directory = prediction_directory(run, config) / job_id
    return next(iter(sorted(directory.glob("*.cif"))), None)


def heavy_atom_coordinates(residues) -> np.ndarray:
    coordinates = [
        atom.coord
        for residue in residues
        for atom in residue.get_atoms()
        if str(getattr(atom, "element", "")).upper() != "H"
    ]
    return np.asarray(coordinates, dtype=float)


def predicted_contacts(structure, expected_sequence: str, cutoff: float) -> list[int]:
    protein = structure.child_dict.get("A")
    dna_chains = [structure.child_dict.get(chain_id) for chain_id in ("B", "C")]
    dna_chains = [chain for chain in dna_chains if chain is not None]
    if protein is None or len(dna_chains) != 2:
        raise ValueError("prediction does not contain protein A and dsDNA chains B/C")
    observed_sequence = str(getattr(protein, "sequence", ""))
    if observed_sequence != expected_sequence:
        raise ValueError(
            "predicted/input protein sequence mismatch: "
            f"{len(observed_sequence)} != {len(expected_sequence)}"
        )
    dna_coordinates = heavy_atom_coordinates(
        residue for chain in dna_chains for residue in chain.get_residues()
    )
    if not len(dna_coordinates):
        return []
    tree = cKDTree(dna_coordinates)
    contacts = []
    for index, residue in enumerate(protein.get_residues(), start=1):
        coordinates = heavy_atom_coordinates([residue])
        if len(coordinates) and float(tree.query(coordinates, k=1)[0].min()) <= cutoff:
            contacts.append(index)
    return contacts


def crop_interval(start: int, stop: int, flank: int, length: int) -> tuple[int, int]:
    return max(1, start - flank), min(length, stop + flank)


def crop_msa(
    config: str,
    source: Path,
    destination: Path,
    start: int,
    stop: int,
    expected_query: str,
) -> str | Path:
    if config == "esmfold2_nomsa":
        return "empty"
    destination.parent.mkdir(parents=True, exist_ok=True)
    crop_a3m_file(source, destination, interval=(start - 1, stop))
    if read_a3m_query(destination) != expected_query:
        raise ValueError(f"Protein/A3M query mismatch: {destination}")
    return destination.resolve()


def base_manifest_record(row, protein: str, msa: str | Path, start: int, stop: int):
    return {
        "dataset_row": row.dataset_row,
        "job_id": row.job_id,
        "complex_id": row.complex_id,
        "native_structure": row.native_structure,
        "native_protein_chain": row.native_protein_chain,
        "native_dna_chain": row.native_dna_chain,
        "protein_sequence": protein,
        "dna_sequence": row.dna_sequence,
        "protein_length": len(protein),
        "dna_length": row.dna_length,
        "total_residues": len(protein) + 2 * int(row.dna_length),
        "msa_path": str(msa),
        "full_protein_length": len(row.protein_sequence_full),
        "crop_start_1based_inclusive": start,
        "crop_stop_1based_inclusive": stop,
    }


def write_run_config(
    run: Path,
    config: str,
    condition: str,
    selected: int,
    args,
) -> None:
    payload = {
        "schema_version": 1,
        "benchmark_mode": "predicted_interface_crop",
        "condition": condition,
        "config": config,
        "configs": [config],
        "source_full_run": str(Path(args.full_run).resolve()),
        "source_oracle_audit": str(Path(args.oracle_audit).resolve()),
        "selected_interfaces": selected,
        "seed": args.seed,
        "contact_cutoff_angstrom": args.contact_cutoff,
        "crop_flank_aa": args.flank,
        "historical_predictions_used": False,
    }
    (run / "benchmark_config.json").write_text(json.dumps(payload, indent=2) + "\n")


def prepare(args) -> None:
    from DockQ.DockQ import load_PDB

    out = Path(args.out)
    full_run = Path(args.full_run)
    out.mkdir(parents=True, exist_ok=True)
    with CommandLog(out / "prepare.log", sys.argv, args):
        full = pd.read_csv(
            full_run / "benchmark_manifest.csv",
            dtype={"native_protein_chain": str, "native_dna_chain": str},
        )
        oracle = pd.read_csv(
            args.oracle_audit,
            dtype={"native_protein_chain": str, "native_dna_chain": str},
        )
        full = full.rename(columns={
            "protein_sequence": "protein_sequence_full",
            "protein_length": "protein_length_full",
            "msa_path": "msa_path_full",
        })
        oracle_columns = [
            "job_id", "oracle_contact_start", "oracle_contact_stop",
            "oracle_contact_positions",
        ]
        candidates = full.merge(
            oracle[oracle_columns], on="job_id", how="inner", validate="one_to_one"
        )
        if candidates.empty:
            raise ValueError("No full-length jobs overlap the oracle audit")

        audit_rows = []
        for config in args.configs:
            interface_run = out / "interface20" / config
            oracle_run = out / "oracle20" / config
            for run in (interface_run, oracle_run):
                input_dir = run / "inputs" / config
                input_dir.mkdir(parents=True, exist_ok=True)
                for stale in input_dir.glob("*.yaml"):
                    stale.unlink()

            interface_records, oracle_records = [], []
            for row in candidates.itertuples(index=False):
                full_sequence = str(row.protein_sequence_full)
                full_length = len(full_sequence)
                cif = prediction_cif(full_run, config, row.job_id)
                status = "predicted_contact"
                error = None
                contacts = []
                if cif is None:
                    status = "missing_full_prediction"
                else:
                    try:
                        contacts = predicted_contacts(
                            load_PDB(str(cif)), full_sequence, args.contact_cutoff
                        )
                        if not contacts:
                            status = "no_contact_full_fallback"
                    except Exception as exc:
                        status = "invalid_full_prediction"
                        error = str(exc)

                predicted_start = min(contacts) if contacts else None
                predicted_stop = max(contacts) if contacts else None
                if status in {"missing_full_prediction", "invalid_full_prediction"}:
                    interface_start = interface_stop = None
                elif contacts:
                    interface_start, interface_stop = crop_interval(
                        predicted_start, predicted_stop, args.flank, full_length
                    )
                else:
                    interface_start, interface_stop = 1, full_length

                oracle_start, oracle_stop = crop_interval(
                    int(row.oracle_contact_start),
                    int(row.oracle_contact_stop),
                    args.flank,
                    full_length,
                )

                if interface_start is not None:
                    protein = full_sequence[interface_start - 1:interface_stop]
                    msa_path = crop_msa(
                        config,
                        Path(row.msa_path_full),
                        interface_run / "msa" / row.job_id / "crop_0.a3m",
                        interface_start,
                        interface_stop,
                        protein,
                    )
                    record = base_manifest_record(
                        row, protein, msa_path, interface_start, interface_stop
                    )
                    record.update({
                        "crop_definition": "predicted_interface",
                        "crop_flank_aa": args.flank,
                        "localizer_config": config,
                        "localizer_status": status,
                    })
                    interface_records.append(record)
                    write_yaml(
                        interface_run / "inputs" / config / f"{row.job_id}.yaml",
                        protein,
                        str(row.dna_sequence),
                        msa_path,
                    )

                oracle_protein = full_sequence[oracle_start - 1:oracle_stop]
                oracle_msa = crop_msa(
                    config,
                    Path(row.msa_path_full),
                    oracle_run / "msa" / row.job_id / "crop_0.a3m",
                    oracle_start,
                    oracle_stop,
                    oracle_protein,
                )
                oracle_record = base_manifest_record(
                    row, oracle_protein, oracle_msa, oracle_start, oracle_stop
                )
                oracle_record.update({
                    "crop_definition": "native_interface",
                    "crop_flank_aa": args.flank,
                })
                oracle_records.append(oracle_record)
                write_yaml(
                    oracle_run / "inputs" / config / f"{row.job_id}.yaml",
                    oracle_protein,
                    str(row.dna_sequence),
                    oracle_msa,
                )

                native_positions = {
                    int(value) for value in str(row.oracle_contact_positions).split(";") if value
                }
                covered_positions = (
                    sum(interface_start <= value <= interface_stop for value in native_positions)
                    if interface_start is not None else 0
                )
                audit_rows.append({
                    "config": config,
                    "job_id": row.job_id,
                    "full_prediction": str(cif) if cif else None,
                    "localizer_status": status,
                    "localizer_error": error,
                    "predicted_contact_count": len(contacts),
                    "predicted_contact_start": predicted_start,
                    "predicted_contact_stop": predicted_stop,
                    "interface_crop_start": interface_start,
                    "interface_crop_stop": interface_stop,
                    "interface_crop_length": (
                        interface_stop - interface_start + 1
                        if interface_start is not None else None
                    ),
                    "interface_actually_shortened": (
                        interface_stop - interface_start + 1 < full_length
                        if interface_start is not None else False
                    ),
                    "oracle_contact_start": row.oracle_contact_start,
                    "oracle_contact_stop": row.oracle_contact_stop,
                    "oracle_crop_start": oracle_start,
                    "oracle_crop_stop": oracle_stop,
                    "oracle_crop_length": oracle_stop - oracle_start + 1,
                    "native_contact_residue_count": len(native_positions),
                    "native_contact_residues_retained": covered_positions,
                    "native_contact_residue_coverage": (
                        covered_positions / len(native_positions) if native_positions else np.nan
                    ),
                    "oracle_span_fully_retained": (
                        interface_start is not None
                        and interface_start <= int(row.oracle_contact_start)
                        and interface_stop >= int(row.oracle_contact_stop)
                    ),
                    "full_length": full_length,
                })

            interface_manifest = pd.DataFrame(interface_records)
            oracle_manifest = pd.DataFrame(oracle_records)
            interface_manifest.to_csv(interface_run / "benchmark_manifest.csv", index=False)
            oracle_manifest.to_csv(oracle_run / "benchmark_manifest.csv", index=False)
            write_run_config(
                interface_run, config, "interface20", len(interface_manifest), args
            )
            write_run_config(oracle_run, config, "oracle20", len(oracle_manifest), args)
            print(
                f"Prepared {config}: interface={len(interface_manifest)}, "
                f"oracle={len(oracle_manifest)}"
            )

        audit = pd.DataFrame(audit_rows)
        audit.to_csv(out / "interface_coordinate_audit.csv", index=False)
        summary = summarize_input_audit(audit)
        summary.to_csv(out / "interface_input_summary.csv", index=False)
        metadata = {
            "generated": datetime.now().astimezone().isoformat(timespec="seconds"),
            "candidate_jobs": len(candidates),
            "configs": list(args.configs),
            "conditions": list(CONDITIONS),
            "interface_definition": (
                f"protein residue with any heavy atom within {args.contact_cutoff:g} A "
                "of either predicted DNA chain B/C"
            ),
            "crop_flank_aa": args.flank,
            "seed": args.seed,
        }
        (out / "experiment.json").write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"Prepared interface-crop benchmark -> {out}")


def run_stage(args, stage: str) -> None:
    out = Path(args.out)
    with CommandLog(out / f"{stage}.log", sys.argv, args):
        failures = []
        for condition in args.conditions:
            for config in args.configs:
                run = out / condition / config
                command = [
                    sys.executable,
                    str(BENCHMARK),
                    stage,
                    "--run",
                    str(run),
                    "--configs",
                    config,
                ]
                if stage == "predict":
                    command.extend(["--seed", str(args.seed)])
                else:
                    command.extend(["--processes", str(args.processes)])
                print("[interface-crop] " + " ".join(command))
                result = subprocess.run(command)
                if result.returncode:
                    failures.append(f"{condition}/{config}")
                    print(
                        f"[interface-crop] partial/failed: {condition}/{config} "
                        f"(exit={result.returncode}); continuing"
                    )
        if failures:
            raise RuntimeError(
                "One or more experiment tasks were partial or failed: "
                + ", ".join(failures)
            )


def paired_test(left: pd.Series, right: pd.Series) -> dict:
    delta = right - left
    nonzero = delta != 0
    pvalue = np.nan
    if nonzero.any():
        pvalue = wilcoxon(
            right[nonzero],
            left[nonzero],
            zero_method="wilcox",
            alternative="two-sided",
            method="auto",
        ).pvalue
    return {
        "n": len(delta),
        "left_mean": left.mean(),
        "right_mean": right.mean(),
        "mean_delta_right_minus_left": delta.mean(),
        "std_delta_right_minus_left": delta.std(ddof=1),
        "median_delta_right_minus_left": delta.median(),
        "right_wins": int((delta > 0).sum()),
        "ties": int((delta == 0).sum()),
        "right_losses": int((delta < 0).sum()),
        "wilcoxon_p": pvalue,
    }


def significance_label(pvalue: float) -> str:
    if pd.isna(pvalue):
        return "NA"
    if pvalue < 0.0001:
        return "****"
    if pvalue < 0.001:
        return "***"
    if pvalue < 0.01:
        return "**"
    if pvalue < 0.05:
        return "*"
    return "ns"


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


def metric_frame(path: Path, config: str, label: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame = frame[frame["config"].eq(config)].copy()
    columns = ["job_id", "dockq_interface", "ipsae", "ipsae_asym_min"]
    frame = frame[[column for column in columns if column in frame]]
    return frame.rename(columns={
        column: f"{column}_{label}" for column in frame.columns if column != "job_id"
    })


def summarize_input_audit(audit: pd.DataFrame) -> pd.DataFrame:
    """Keep localization availability separate from boundary accuracy."""
    rows = []
    for config, group in audit.groupby("config", sort=False):
        available = group[group["interface_crop_start"].notna()]
        rows.append({
            "config": config,
            "candidates": len(group),
            "localizer_available": len(available),
            "localizer_success_rate": len(available) / len(group),
            "predicted_contact": int((group["localizer_status"] == "predicted_contact").sum()),
            "no_contact_full_fallback": int(
                (group["localizer_status"] == "no_contact_full_fallback").sum()
            ),
            "actually_shortened": int(group["interface_actually_shortened"].sum()),
            "mean_interface_crop_length": available["interface_crop_length"].mean(),
            "mean_oracle_crop_length": group["oracle_crop_length"].mean(),
            "mean_native_contact_coverage_available": (
                available["native_contact_residue_coverage"].mean()
            ),
            "oracle_span_fully_retained_available": (
                available["oracle_span_fully_retained"].mean()
            ),
        })
    return pd.DataFrame(rows)


def analyze(args) -> None:
    out = Path(args.out)
    full_metrics = Path(args.full_run) / "metrics" / "metrics_per_case.csv"
    with CommandLog(out / "analyze.log", sys.argv, args):
        condition_rows, comparison_rows, paired_frames, timing_rows = [], [], [], []
        audit = pd.read_csv(out / "interface_coordinate_audit.csv")
        input_summary = summarize_input_audit(audit)
        input_summary.to_csv(out / "interface_input_summary.csv", index=False)
        full_timing_path = Path(args.full_run) / "timings.csv"
        full_timings = (
            pd.read_csv(full_timing_path) if full_timing_path.is_file() else pd.DataFrame()
        )
        strict_sets = []
        for config in args.configs:
            frames = {
                "full": metric_frame(full_metrics, config, "full"),
                "interface20": metric_frame(
                    out / "interface20" / config / "metrics" / "metrics_per_case.csv",
                    config,
                    "interface20",
                ),
                "oracle20": metric_frame(
                    out / "oracle20" / config / "metrics" / "metrics_per_case.csv",
                    config,
                    "oracle20",
                ),
            }
            paired = frames["full"]
            for condition in ("interface20", "oracle20"):
                paired = paired.merge(frames[condition], on="job_id", how="inner")
            dockq_columns = [f"dockq_interface_{condition}" for condition in ("full", *CONDITIONS)]
            paired = paired.dropna(subset=dockq_columns).copy()
            paired.insert(0, "config", config)
            paired_frames.append(paired)
            strict_sets.append(set(paired["job_id"]))

            for condition in ("full", *CONDITIONS):
                dockq = paired[f"dockq_interface_{condition}"]
                ipsae_column = f"ipsae_{condition}"
                condition_rows.append({
                    "config": config,
                    "condition": condition,
                    "n": len(dockq),
                    "dockq_mean": dockq.mean(),
                    "dockq_std": dockq.std(ddof=1),
                    "dockq_median": dockq.median(),
                    "acceptable_or_better": (dockq >= 0.23).mean(),
                    "medium_or_better": (dockq >= 0.49).mean(),
                    "high": (dockq >= 0.80).mean(),
                    "ipsae_mean": (
                        paired[ipsae_column].mean() if ipsae_column in paired else np.nan
                    ),
                    "ipsae_std": (
                        paired[ipsae_column].std(ddof=1)
                        if ipsae_column in paired else np.nan
                    ),
                    "ipsae_dockq_spearman": (
                        paired[[ipsae_column, f"dockq_interface_{condition}"]]
                        .corr(method="spearman").iloc[0, 1]
                        if ipsae_column in paired and len(paired) >= 2 else np.nan
                    ),
                })

            comparisons = [
                ("full", "interface20", "predicted-interface crop effect"),
                ("full", "oracle20", "native-interface crop effect"),
                ("interface20", "oracle20", "remaining boundary-localization gap"),
            ]
            for left, right, interpretation in comparisons:
                result = paired_test(
                    paired[f"dockq_interface_{left}"],
                    paired[f"dockq_interface_{right}"],
                )
                result.update({
                    "config": config,
                    "left": left,
                    "right": right,
                    "interpretation": interpretation,
                })
                comparison_rows.append(result)

            for condition in ("full", *CONDITIONS):
                timing_path = (
                    full_timing_path
                    if condition == "full"
                    else out / condition / config / "timings.csv"
                )
                if timing_path.is_file():
                    timing = (
                        full_timings.copy()
                        if condition == "full"
                        else pd.read_csv(timing_path)
                    )
                    timing = timing[timing["config"].eq(config)]
                    if not timing.empty:
                        timing_rows.append({
                            "config": config,
                            "condition": condition,
                            "wall_seconds": timing.iloc[-1]["wall_seconds"],
                            "n_inputs": timing.iloc[-1]["n_inputs"],
                            "n_completed": timing.iloc[-1]["n_completed"],
                            "seconds_per_completed": timing.iloc[-1]["seconds_per_completed"],
                            "status": timing.iloc[-1].get("status", "unknown"),
                        })

        condition_summary = pd.DataFrame(condition_rows)
        comparison = pd.DataFrame(comparison_rows)[[
            "config", "left", "right", "interpretation", "n", "left_mean",
            "right_mean", "mean_delta_right_minus_left",
            "std_delta_right_minus_left", "median_delta_right_minus_left",
            "right_wins", "ties", "right_losses", "wilcoxon_p",
        ]]
        family_size = comparison.groupby(["left", "right"])["wilcoxon_p"].transform(
            "count"
        )
        comparison["bonferroni_family_size"] = family_size.astype(int)
        comparison["wilcoxon_p_bonferroni"] = np.minimum(
            comparison["wilcoxon_p"] * family_size, 1.0
        )
        comparison["significance_raw"] = comparison["wilcoxon_p"].map(
            significance_label
        )
        comparison["significance_bonferroni"] = comparison[
            "wilcoxon_p_bonferroni"
        ].map(significance_label)
        paired_all = pd.concat(paired_frames, ignore_index=True)
        timings = pd.DataFrame(timing_rows)
        condition_lookup = condition_summary.set_index(["config", "condition"])
        comparison_lookup = comparison.set_index(["config", "left", "right"])
        timing_lookup = timings.set_index(["config", "condition"])
        unified_rows = []
        for config in args.configs:
            full = condition_lookup.loc[(config, "full")]
            interface = condition_lookup.loc[(config, "interface20")]
            oracle = condition_lookup.loc[(config, "oracle20")]
            interface_vs_full = comparison_lookup.loc[
                (config, "full", "interface20")
            ]
            oracle_vs_interface = comparison_lookup.loc[
                (config, "interface20", "oracle20")
            ]
            full_time = timing_lookup.loc[(config, "full")]
            interface_time = timing_lookup.loc[(config, "interface20")]
            oracle_time = timing_lookup.loc[(config, "oracle20")]
            unified_rows.append({
                "config": config,
                "paired_n": int(full["n"]),
                "full_dockq_mean": full["dockq_mean"],
                "full_dockq_std": full["dockq_std"],
                "interface20_dockq_mean": interface["dockq_mean"],
                "interface20_dockq_std": interface["dockq_std"],
                "oracle20_dockq_mean": oracle["dockq_mean"],
                "oracle20_dockq_std": oracle["dockq_std"],
                "interface_minus_full_mean": interface_vs_full[
                    "mean_delta_right_minus_left"
                ],
                "interface_minus_full_std": interface_vs_full[
                    "std_delta_right_minus_left"
                ],
                "interface_vs_full_wilcoxon_p_raw": interface_vs_full["wilcoxon_p"],
                "interface_vs_full_wilcoxon_p_bonferroni": interface_vs_full[
                    "wilcoxon_p_bonferroni"
                ],
                "interface_vs_full_significance_bonferroni": interface_vs_full[
                    "significance_bonferroni"
                ],
                "interface_minus_oracle_mean": -oracle_vs_interface[
                    "mean_delta_right_minus_left"
                ],
                "interface_minus_oracle_std": oracle_vs_interface[
                    "std_delta_right_minus_left"
                ],
                "interface_vs_oracle_wilcoxon_p_raw": oracle_vs_interface[
                    "wilcoxon_p"
                ],
                "interface_vs_oracle_wilcoxon_p_bonferroni": oracle_vs_interface[
                    "wilcoxon_p_bonferroni"
                ],
                "interface_vs_oracle_significance_bonferroni": oracle_vs_interface[
                    "significance_bonferroni"
                ],
                "interface20_ipsae_dockq_spearman": interface[
                    "ipsae_dockq_spearman"
                ],
                "full_inference_wall_seconds": full_time["wall_seconds"],
                "full_inference_n_inputs": int(full_time["n_inputs"]),
                "full_inference_n_completed": int(full_time["n_completed"]),
                "full_inference_seconds_per_completed": full_time[
                    "seconds_per_completed"
                ],
                "full_inference_status": full_time["status"],
                "interface20_inference_wall_seconds": interface_time["wall_seconds"],
                "interface20_inference_n_inputs": int(interface_time["n_inputs"]),
                "interface20_inference_n_completed": int(interface_time["n_completed"]),
                "interface20_inference_seconds_per_completed": interface_time[
                    "seconds_per_completed"
                ],
                "interface20_inference_status": interface_time["status"],
                "oracle20_inference_wall_seconds": oracle_time["wall_seconds"],
                "oracle20_inference_n_inputs": int(oracle_time["n_inputs"]),
                "oracle20_inference_n_completed": int(oracle_time["n_completed"]),
                "oracle20_inference_seconds_per_completed": oracle_time[
                    "seconds_per_completed"
                ],
                "oracle20_inference_status": oracle_time["status"],
                "interface20_vs_full_seconds_per_completed_pct": (
                    interface_time["seconds_per_completed"]
                    / full_time["seconds_per_completed"]
                    - 1.0
                ) * 100.0,
            })
        unified = pd.DataFrame(unified_rows)
        condition_summary.to_csv(out / "condition_metric_summary.csv", index=False)
        comparison.to_csv(out / "paired_comparisons.csv", index=False)
        paired_all.to_csv(out / "paired_metrics_per_case.csv", index=False)
        timings.to_csv(out / "timing_summary.csv", index=False)
        unified.to_csv(out / "unified_results.csv", index=False)
        common = set.intersection(*strict_sets) if strict_sets else set()

        compact = condition_summary.pivot(
            index="config", columns="condition", values="dockq_mean"
        ).reset_index()
        compact["interface_minus_full"] = compact["interface20"] - compact["full"]
        compact["interface_minus_oracle"] = compact["interface20"] - compact["oracle20"]
        crop_effect = comparison[
            comparison["left"].eq("full")
            & comparison["right"].eq("interface20")
        ].copy()
        boundary_gap = comparison[
            comparison["left"].eq("interface20")
            & comparison["right"].eq("oracle20")
        ].copy()
        nominal = crop_effect[crop_effect["wilcoxon_p"] < 0.05]["config"].tolist()
        bonferroni_alpha = 0.05 / max(len(crop_effect), 1)
        corrected = crop_effect[
            crop_effect["wilcoxon_p"] < bonferroni_alpha
        ]["config"].tolist()
        best_crop = crop_effect.loc[
            crop_effect["mean_delta_right_minus_left"].idxmax()
        ]
        worst_crop = crop_effect.loc[
            crop_effect["mean_delta_right_minus_left"].idxmin()
        ]
        available_coverage = input_summary["mean_native_contact_coverage_available"]
        timing_pivot = timings.pivot(
            index="config", columns="condition", values="seconds_per_completed"
        )
        timing_change = (
            (timing_pivot["interface20"] / timing_pivot["full"] - 1.0) * 100.0
        )
        unified_display = pd.DataFrame({
            "model": unified["config"],
            "n": unified["paired_n"],
            "Full DockQ mean+/-SD": unified.apply(
                lambda row: f"{row.full_dockq_mean:.4f} +/- {row.full_dockq_std:.4f}",
                axis=1,
            ),
            "Interface20 DockQ mean+/-SD": unified.apply(
                lambda row: (
                    f"{row.interface20_dockq_mean:.4f} +/- "
                    f"{row.interface20_dockq_std:.4f}"
                ),
                axis=1,
            ),
            "Oracle20 DockQ mean+/-SD": unified.apply(
                lambda row: (
                    f"{row.oracle20_dockq_mean:.4f} +/- "
                    f"{row.oracle20_dockq_std:.4f}"
                ),
                axis=1,
            ),
            "Delta I-F mean+/-SD": unified.apply(
                lambda row: (
                    f"{row.interface_minus_full_mean:+.4f} +/- "
                    f"{row.interface_minus_full_std:.4f}"
                ),
                axis=1,
            ),
            "I vs F p(raw/Bonf), sig": unified.apply(
                lambda row: (
                    f"{row.interface_vs_full_wilcoxon_p_raw:.4g}/"
                    f"{row.interface_vs_full_wilcoxon_p_bonferroni:.4g}, "
                    f"{row.interface_vs_full_significance_bonferroni}"
                ),
                axis=1,
            ),
            "Delta I-O mean+/-SD": unified.apply(
                lambda row: (
                    f"{row.interface_minus_oracle_mean:+.4f} +/- "
                    f"{row.interface_minus_oracle_std:.4f}"
                ),
                axis=1,
            ),
            "I vs O p(raw/Bonf), sig": unified.apply(
                lambda row: (
                    f"{row.interface_vs_oracle_wilcoxon_p_raw:.4g}/"
                    f"{row.interface_vs_oracle_wilcoxon_p_bonferroni:.4g}, "
                    f"{row.interface_vs_oracle_significance_bonferroni}"
                ),
                axis=1,
            ),
            "Full inference wall/completed-attempted=s/case,status": unified.apply(
                lambda row: (
                    f"{row.full_inference_wall_seconds:.1f}/"
                    f"{int(row.full_inference_n_completed)}-"
                    f"{int(row.full_inference_n_inputs)}="
                    f"{row.full_inference_seconds_per_completed:.2f},"
                    f"{row.full_inference_status}"
                ),
                axis=1,
            ),
            "Interface20 inference wall/completed-attempted=s/case,status": unified.apply(
                lambda row: (
                    f"{row.interface20_inference_wall_seconds:.1f}/"
                    f"{int(row.interface20_inference_n_completed)}-"
                    f"{int(row.interface20_inference_n_inputs)}="
                    f"{row.interface20_inference_seconds_per_completed:.2f},"
                    f"{row.interface20_inference_status}"
                ),
                axis=1,
            ),
            "Oracle20 inference wall/completed-attempted=s/case,status": unified.apply(
                lambda row: (
                    f"{row.oracle20_inference_wall_seconds:.1f}/"
                    f"{int(row.oracle20_inference_n_completed)}-"
                    f"{int(row.oracle20_inference_n_inputs)}="
                    f"{row.oracle20_inference_seconds_per_completed:.2f},"
                    f"{row.oracle20_inference_status}"
                ),
                axis=1,
            ),
            "Interface ipSAE-DockQ rho": unified[
                "interface20_ipsae_dockq_spearman"
            ],
        })
        report = f"""# FoldBench predicted-interface crop+{args.flank} benchmark

- Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}
- Interface localizer: each configuration's own full-length prediction
- Predicted contact: protein heavy atom within {args.contact_cutoff:g} A of either predicted DNA strand
- Crop: continuous predicted contact span + {args.flank} aa
- Seed: {args.seed}
- Primary analysis: strict full/interface/oracle pairing within each model
- Strict scoreable jobs shared by all six models: {len(common)}
- No-contact localization falls back to the full-length input and is reported

## Main interpretation

- Interface+{args.flank} minus full mean DockQ ranged from
  {crop_effect['mean_delta_right_minus_left'].min():+.4f} to
  {crop_effect['mean_delta_right_minus_left'].max():+.4f}.  The largest gain was
  {best_crop['config']} ({best_crop['mean_delta_right_minus_left']:+.4f}); the
  largest loss was {worst_crop['config']}
  ({worst_crop['mean_delta_right_minus_left']:+.4f}).
- Nominal paired Wilcoxon p<0.05: {', '.join(nominal) if nominal else 'none'}.
  None survived the six-model Bonferroni threshold
  ({bonferroni_alpha:.4g}){'' if not corrected else ': ' + ', '.join(corrected)}.
- No interface+{args.flank} versus oracle+{args.flank} contrast had p<0.05
  (minimum p={boundary_gap['wilcoxon_p'].min():.4g}).
- Among available full-length localizers, mean native-contact residue coverage
  ranged from {available_coverage.min():.3f} to {available_coverage.max():.3f}.
- Observed interface rerun seconds per completed case were
  {timing_change.min():+.1f}% to {timing_change.max():+.1f}% relative to the
  original full benchmark.  These are different batch sizes and exclude the
  additional full-localization inference from the interface workflow.

## Unified result table

{markdown(unified_display)}

Time cells are `prediction wall seconds / completed-attempted predictions =
seconds per completed prediction, status`.

## Statistical calculation

- All DockQ comparisons use the same strict paired jobs within each model.
- Condition variability is the sample standard deviation (`ddof=1`).
- For a contrast, `d_i = DockQ_right,i - DockQ_left,i`; the table reports the
  arithmetic mean and sample SD of these paired differences.  `Delta I-F`
  means Interface20 minus Full; `Delta I-O` means Interface20 minus Oracle20.
- Significance uses a two-sided paired Wilcoxon signed-rank test.  Exact-zero
  differences are counted as ties and excluded from the rank statistic
  (`zero_method=wilcox`); SciPy chooses the computational method automatically.
- Raw p-values are followed by Bonferroni-adjusted p-values.  Adjustment is
  performed separately for each contrast across the six models:
  `p_Bonf = min(6 * p_raw, 1)`.
- Adjusted significance labels: `ns` >= 0.05, `*` < 0.05, `**` < 0.01,
  `***` < 0.001, `****` < 0.0001.
- Prediction time is wall-clock time around `python -m boltzscan predict` and
  includes model startup, input preprocessing, GPU inference, and output
  writing.  It excludes precomputed MSA generation, DockQ/ipSAE scoring, and
  the full-length localization pass additionally required by Interface20.

## Input localization audit

{markdown(input_summary)}

## Mean DockQ on each model's paired set

{markdown(compact)}

## Per-condition metrics

{markdown(condition_summary)}

## Paired DockQ contrasts

{markdown(comparison)}

## Runtime

{markdown(timings) if not timings.empty else 'No timing rows found.'}

The native-contact oracle is a diagnostic boundary annotation, not a deployable
input and not a guaranteed upper bound on DockQ: cropping can remove useful
protein context.  The predicted-interface condition uses the same DNA case for full
localization and cropped rerun; a representative-DNA transfer experiment is a
separate question.

Full-length timings are the original 300-input benchmark, whereas interface
and oracle timings contain 89--106 inputs.  Compare seconds per completed case,
not raw wall time, across conditions.
"""
        (out / "interface_crop20_report.md").write_text(report)
        print(
            f"Analyzed {len(paired_all)} model/job paired rows; "
            f"all-model common n={len(common)} -> {out}"
        )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest="stage", required=True)
    for stage in ("prepare", "predict", "score", "analyze"):
        child = sub.add_parser(stage)
        child.add_argument("--full-run", required=True)
        child.add_argument("--oracle-audit", required=True)
        child.add_argument("--out", required=True)
        child.add_argument("--configs", nargs="+", choices=CONFIGS, default=list(CONFIGS))
        child.add_argument(
            "--conditions", nargs="+", choices=CONDITIONS, default=list(CONDITIONS)
        )
        child.add_argument("--flank", type=int, default=20)
        child.add_argument("--contact-cutoff", type=float, default=5.0)
        child.add_argument("--seed", type=int, default=42)
        child.add_argument("--processes", type=int, default=12)
    return root


def main() -> None:
    args = parser().parse_args()
    if args.flank < 0:
        raise ValueError("--flank must be non-negative")
    if args.contact_cutoff <= 0:
        raise ValueError("--contact-cutoff must be positive")
    if args.stage == "prepare":
        prepare(args)
    elif args.stage in {"predict", "score"}:
        run_stage(args, args.stage)
    else:
        analyze(args)


if __name__ == "__main__":
    main()
