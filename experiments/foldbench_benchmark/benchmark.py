#!/usr/bin/env python
"""Prepare, run, and score the experiment-only FoldBench benchmark."""

from __future__ import annotations

import argparse
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
from difflib import SequenceMatcher
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import yaml


CONFIGS = (
    "esmfold2_msa",
    "esmfold2_nomsa",
    "boltz1",
    "boltz1_ode",
    "boltz2",
    "boltz2_ode",
)
MODELS = {
    "esmfold2_msa": "esmfold2",
    "esmfold2_nomsa": "esmfold2",
    "boltz1": "boltz1",
    "boltz1_ode": "boltz1_ode",
    "boltz2": "boltz2",
    "boltz2_ode": "boltz2_ode",
}
REQUIRED_COLUMNS = {
    "cif_id",
    "protein_chain_id",
    "protein_sequence",
    "dna_chain_id",
    "dna_sequence",
}
DNA_RESIDUES = frozenset({"DA", "DC", "DG", "DT", "DI", "A", "C", "G", "T", "I"})
DEFAULT_ESMFOLD2_WEIGHTS = Path("/media/zlab/lexar4t/esm_weights/ESMFold2")
DEFAULT_ESMC_WEIGHTS = Path("/media/zlab/lexar4t/esm_weights/ESMC-6B")
HERE = Path(__file__).resolve().parent


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
            stream.flush()
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


class StageLog:
    def __init__(self, path, command):
        self.path = Path(path)
        self.command = command

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = self.path.open("w", encoding="utf-8")
        self.started = time.perf_counter()
        self.handle.write("FoldBench experiment".center(72, "=") + "\n")
        self.handle.write(f"Start:   {datetime.now().astimezone().isoformat(timespec='seconds')}\n")
        self.handle.write(f"Command: {self.command}\nOutput:\n")
        self.stdout, self.stderr = sys.stdout, sys.stderr
        self.stdout_tee = Tee(self.stdout, self.handle)
        self.stderr_tee = Tee(self.stderr, self.handle)
        self.stdout_context = redirect_stdout(self.stdout_tee)
        self.stderr_context = redirect_stderr(self.stderr_tee)
        self.stdout_context.__enter__()
        self.stderr_context.__enter__()
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.stderr_context.__exit__(exc_type, exc, traceback)
        self.stdout_context.__exit__(exc_type, exc, traceback)
        self.handle.write(f"Status:  {'FAILED' if exc_type else 'OK'}\n")
        if exc is not None:
            self.handle.write(f"Error:   {exc}\n")
        self.handle.write(f"Elapsed: {time.perf_counter() - self.started:.2f} seconds\n")
        self.handle.close()
        return False


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def safe_chain_id(chain_id):
    return "".join(char if char.isalnum() or char in {"-", "."} else "-" for char in chain_id)


def msa_query_sequence(path):
    query = []
    in_query = False
    with Path(path).open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if in_query:
                    break
                in_query = True
                continue
            if in_query:
                query.append(line)
    if not query:
        raise ValueError(f"MSA contains no query sequence: {path}")
    return "".join(
        char for char in "".join(query) if char != "-" and not char.islower()
    ).upper()


def native_path(native_dir, complex_id):
    for suffix in (".cif", ".pdb"):
        path = native_dir / f"{complex_id}{suffix}"
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError(f"Native structure not found for {complex_id} under {native_dir}")


def raw_manifest(data_dir):
    data_dir = Path(data_dir).resolve()
    table_path = data_dir / "foldbench_pep_dna.csv"
    msa_dir = data_dir / "msa"
    native_dir = data_dir / "ground_truth_1522" / "ground_truth_20250520"
    table = pd.read_csv(table_path, dtype=str)
    missing = sorted(REQUIRED_COLUMNS.difference(table.columns))
    if missing:
        raise ValueError("FoldBench table is missing columns: " + ", ".join(missing))

    records, excluded = [], []
    for dataset_row, row in table.iterrows():
        complex_id = str(row["cif_id"])
        protein_chain = str(row["protein_chain_id"])
        dna_chain = str(row["dna_chain_id"])
        protein = str(row["protein_sequence"]).upper()
        dna = str(row["dna_sequence"]).upper()
        invalid_protein = sorted(set(protein).difference("ARNDCEQGHILKMFPSTWYVXJBZOU"))
        invalid_dna = sorted(set(dna).difference("ACGTN"))
        if invalid_protein or invalid_dna:
            reasons = []
            if invalid_protein:
                reasons.append(f"illegal protein characters: {','.join(invalid_protein)}")
            if invalid_dna:
                reasons.append(f"illegal DNA characters: {','.join(invalid_dna)}")
            excluded.append({
                "dataset_row": dataset_row,
                "complex_id": complex_id,
                "native_protein_chain": protein_chain,
                "native_dna_chain": dna_chain,
                "reason": "; ".join(reasons),
            })
            continue

        base_id = f"{complex_id}_{protein_chain}"
        job_id = f"{base_id}__dna_{safe_chain_id(dna_chain)}"
        msa_path = (msa_dir / base_id / "0.a3m").resolve()
        if not msa_path.is_file():
            raise FileNotFoundError(f"MSA not found for {base_id}: {msa_path}")
        query = msa_query_sequence(msa_path)
        if query != protein:
            raise ValueError(
                f"Protein/MSA query mismatch for {job_id}: "
                f"protein length {len(protein)}, query length {len(query)}"
            )
        records.append({
            "dataset_row": int(dataset_row),
            "job_id": job_id,
            "complex_id": complex_id,
            "native_structure": str(native_path(native_dir, complex_id)),
            "native_protein_chain": protein_chain,
            "native_dna_chain": dna_chain,
            "protein_sequence": protein,
            "dna_sequence": dna,
            "protein_length": len(protein),
            "dna_length": len(dna),
            "total_residues": len(protein) + 2 * len(dna),
            "msa_path": str(msa_path),
        })
    manifest = pd.DataFrame(records)
    if manifest["job_id"].duplicated().any():
        raise ValueError("Generated FoldBench job IDs are not unique")
    excluded = pd.DataFrame(excluded, columns=[
        "dataset_row", "complex_id", "native_protein_chain", "native_dna_chain", "reason",
    ])
    return manifest, excluded


def length_stratified(manifest, limit):
    if limit is None or limit >= len(manifest):
        return manifest.copy()
    if limit < 1:
        raise ValueError("--limit must be positive")
    ordered = manifest.sort_values(["total_residues", "job_id"], kind="mergesort").reset_index(drop=True)
    selected = ordered.iloc[np.linspace(0, len(ordered) - 1, num=limit, dtype=int)]
    return selected.sort_values("dataset_row", kind="mergesort").reset_index(drop=True)


def write_yaml(path, protein, dna, msa):
    payload = {
        "version": 1,
        "sequences": [
            {"protein": {"id": ["A"], "sequence": protein, "msa": str(msa)}},
            {"dna": {"id": ["B"], "sequence": dna}},
            {"dna": {"id": ["C"], "sequence": reverse_complement(dna)}},
        ],
    }
    path.write_text(yaml.safe_dump(payload, sort_keys=False))


def write_input_report(run_dir, data_dir, manifest, excluded, limit):
    protein, dna = manifest["protein_length"], manifest["dna_length"]
    report = f"""# FoldBench input audit

- Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}
- Raw source: `{Path(data_dir).resolve()}`
- Raw interface rows: {len(manifest) + len(excluded)}
- Eligible interfaces: {len(manifest)}
- Excluded noncanonical interfaces: {len(excluded)}
- Unique native complexes: {manifest['complex_id'].nunique()}
- Length-stratified limit: {limit if limit is not None else 'none (full set)'}
- Historical FoldBench predictions used: no

| polymer | min | median | mean | p95 | max |
|---|---:|---:|---:|---:|---:|
| protein | {protein.min()} | {protein.median():.1f} | {protein.mean():.1f} | {protein.quantile(.95):.1f} | {protein.max()} |
| modeled DNA strand | {dna.min()} | {dna.median():.1f} | {dna.mean():.1f} | {dna.quantile(.95):.1f} | {dna.max()} |

Every retained row has legal production characters, a native structure, a
local A3M, an exact protein/MSA query match, and a unique job ID. Excluded rows
are recorded without sequence mutation in `excluded_interfaces.csv`.
"""
    (run_dir / "foldbench_eda_report.md").write_text(report)


def prepare(args):
    run_dir = Path(args.run)
    run_dir.mkdir(parents=True, exist_ok=True)
    with StageLog(run_dir / "prepare.log", " ".join(sys.argv)):
        complete, excluded = raw_manifest(args.data)
        manifest = length_stratified(complete, args.limit)
        manifest.to_csv(run_dir / "benchmark_manifest.csv", index=False)
        excluded.to_csv(run_dir / "excluded_interfaces.csv", index=False)
        for config in CONFIGS:
            input_dir = run_dir / "inputs" / config
            input_dir.mkdir(parents=True, exist_ok=True)
            for stale in input_dir.glob("*.yaml"):
                stale.unlink()
            for row in manifest.itertuples(index=False):
                msa = "empty" if config == "esmfold2_nomsa" else row.msa_path
                write_yaml(
                    input_dir / f"{row.job_id}.yaml",
                    row.protein_sequence,
                    row.dna_sequence,
                    msa,
                )
        config = {
            "schema_version": 1,
            "data_dir": str(Path(args.data).resolve()),
            "configs": list(CONFIGS),
            "raw_interfaces": len(complete) + len(excluded),
            "eligible_interfaces": len(complete),
            "selected_interfaces": len(manifest),
            "excluded_interfaces": len(excluded),
            "limit": args.limit,
            "historical_predictions_used": False,
        }
        (run_dir / "benchmark_config.json").write_text(json.dumps(config, indent=2) + "\n")
        write_input_report(run_dir, args.data, manifest, excluded, args.limit)
        print(f"Prepared {len(manifest)} fresh interfaces x {len(CONFIGS)} configurations -> {run_dir}")


def _write_fasta(path, manifest):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in manifest.itertuples(index=False):
            handle.write(f">{row.job_id}\n")
            sequence = str(row.protein_sequence)
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")


def _parse_whitelist_domtbl(domtbl, whitelist):
    allowed = set(
        pd.read_csv(whitelist, sep="\t", dtype=str)["pfam_acc"].astype(str)
    )
    records = []
    with Path(domtbl).open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            accession = fields[4].split(".")[0]
            if accession in allowed:
                records.append((fields[0], accession, int(fields[19]), int(fields[20])))
    return records


def write_crop_report(run_dir, source_run, source_manifest, manifest, excluded, flank):
    full_lengths = pd.to_numeric(manifest["full_protein_length"])
    crop_lengths = pd.to_numeric(manifest["protein_length"])
    report = f"""# FoldBench DBD crop input audit

- Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}
- Source full-length run: `{Path(source_run).resolve()}`
- Source interfaces: {len(source_manifest)}
- Interfaces with a curated DBD Pfam hit: {len(manifest)}
- Interfaces without a curated DBD Pfam hit: {len(excluded)}
- Crop rule: union of all DBD Pfam envelope coordinates plus {flank} aa per side
- Full-length fallback used: no
- MSA generated during prediction: no

| protein form | min | median | mean | p95 | max |
|---|---:|---:|---:|---:|---:|
| full length | {full_lengths.min()} | {full_lengths.median():.1f} | {full_lengths.mean():.1f} | {full_lengths.quantile(.95):.1f} | {full_lengths.max()} |
| DBD+{flank} | {crop_lengths.min()} | {crop_lengths.median():.1f} | {crop_lengths.mean():.1f} | {crop_lengths.quantile(.95):.1f} | {crop_lengths.max()} |

Every retained A3M was cropped with the same 1-based inclusive protein
coordinates, and its ungapped query was verified against the cropped protein.
"""
    (Path(run_dir) / "foldbench_crop_audit_report.md").write_text(report)


def prepare_crop(args):
    from boltzscan.msa import crop_a3m_file, read_a3m_query
    from boltzscan.pwmmap.dbd import dbd_crop_interval, parse_domtbl_dbds

    source_run = Path(args.source_run)
    _, source_manifest = load_run(source_run)
    run_dir = Path(args.run)
    run_dir.mkdir(parents=True, exist_ok=True)
    with StageLog(run_dir / "prepare-crop.log", " ".join(sys.argv)):
        dbd_dir = run_dir / "dbd"
        protein_fasta = dbd_dir / "proteins.fasta"
        _write_fasta(protein_fasta, source_manifest)
        find_tf_command = [
            sys.executable, "-m", "boltzscan", "find-tf",
            "--proteins", str(protein_fasta),
            "--output", str(dbd_dir),
            "--pfam", str(Path(args.pfam).resolve()),
            "--cpu", str(args.cpu),
            "--force",
        ]
        print("[experiment] " + " ".join(find_tf_command))
        subprocess.run(find_tf_command, check=True)

        domtbl = dbd_dir / "pfam.domtbl"
        dbd_records = (
            _parse_whitelist_domtbl(domtbl, args.dbd_whitelist)
            if args.dbd_whitelist
            else parse_domtbl_dbds(domtbl)
        )
        by_job = {}
        for job_id, pfam_acc, start, stop in dbd_records:
            by_job.setdefault(job_id, []).append((pfam_acc, int(start), int(stop)))

        retained, excluded = [], []
        for row in source_manifest.itertuples(index=False):
            hits = by_job.get(row.job_id, [])
            if not hits:
                excluded.append({
                    "dataset_row": row.dataset_row,
                    "job_id": row.job_id,
                    "complex_id": row.complex_id,
                    "reason": "no curated DBD Pfam hit",
                })
                continue
            crop_lo, crop_hi = dbd_crop_interval(
                dbd_records,
                row.job_id,
                protein_length=len(row.protein_sequence),
                flank=args.flank,
            )
            crop_sequence = row.protein_sequence[crop_lo - 1:crop_hi]
            crop_msa = (run_dir / "msa" / row.job_id / "crop_0.a3m").resolve()
            crop_a3m_file(
                Path(row.msa_path),
                crop_msa,
                interval=(crop_lo - 1, crop_hi),
            )
            if read_a3m_query(crop_msa) != crop_sequence:
                raise ValueError(f"Cropped protein/MSA query mismatch for {row.job_id}")
            record = row._asdict()
            record.update({
                "full_protein_sequence": row.protein_sequence,
                "full_protein_length": len(row.protein_sequence),
                "protein_sequence": crop_sequence,
                "protein_length": len(crop_sequence),
                "total_residues": len(crop_sequence) + 2 * len(row.dna_sequence),
                "msa_path": str(crop_msa),
                "crop_start_1based_inclusive": crop_lo,
                "crop_stop_1based_inclusive": crop_hi,
                "dbd_flank_aa": args.flank,
                "dbd_hits": ";".join(
                    f"{pfam_acc}:{start}-{stop}" for pfam_acc, start, stop in hits
                ),
            })
            retained.append(record)

        manifest = pd.DataFrame(retained)
        if manifest.empty:
            raise ValueError("No FoldBench interface has a curated DBD Pfam hit")
        excluded = pd.DataFrame(excluded)
        manifest.to_csv(run_dir / "benchmark_manifest.csv", index=False)
        excluded.to_csv(run_dir / "excluded_no_dbd.csv", index=False)
        for config in CONFIGS:
            input_dir = run_dir / "inputs" / config
            input_dir.mkdir(parents=True, exist_ok=True)
            for stale in input_dir.glob("*.yaml"):
                stale.unlink()
            for row in manifest.itertuples(index=False):
                msa = "empty" if config == "esmfold2_nomsa" else row.msa_path
                write_yaml(
                    input_dir / f"{row.job_id}.yaml",
                    row.protein_sequence,
                    row.dna_sequence,
                    msa,
                )
        config = {
            "schema_version": 1,
            "benchmark_mode": "dbd_crop",
            "source_run": str(source_run.resolve()),
            "data_dir": str(Path(load_run(source_run)[0]["data_dir"]).resolve()),
            "configs": list(CONFIGS),
            "source_interfaces": len(source_manifest),
            "selected_interfaces": len(manifest),
            "excluded_no_dbd": len(excluded),
            "crop_flank_aa": args.flank,
            "pfam": str(Path(args.pfam).resolve()),
            "dbd_whitelist": (
                str(Path(args.dbd_whitelist).resolve()) if args.dbd_whitelist else None
            ),
            "historical_predictions_used": False,
        }
        (run_dir / "benchmark_config.json").write_text(json.dumps(config, indent=2) + "\n")
        write_crop_report(
            run_dir, source_run, source_manifest, manifest, excluded, args.flank
        )
        print(
            f"Prepared {len(manifest)}/{len(source_manifest)} DBD+{args.flank} interfaces "
            f"x {len(CONFIGS)} configurations -> {run_dir}"
        )


def load_run(run_dir):
    run_dir = Path(run_dir)
    config_path = run_dir / "benchmark_config.json"
    manifest_path = run_dir / "benchmark_manifest.csv"
    if not config_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError(f"Experiment run is not prepared: {run_dir}")
    config = json.loads(config_path.read_text())
    manifest = pd.read_csv(
        manifest_path,
        dtype={"native_protein_chain": str, "native_dna_chain": str},
    )
    return config, manifest


def predictions_dir(run_dir, config):
    prefix = "esmfold_results_" if config.startswith("esmfold2") else "boltz_results_"
    return Path(run_dir) / "outputs" / config / f"{prefix}{config}" / "predictions"


def completed_predictions(path):
    path = Path(path)
    if not path.is_dir():
        return 0
    return sum(1 for child in path.iterdir() if child.is_dir() and any(child.glob("*.cif")))


def save_timing(run_dir, row):
    path = Path(run_dir) / "timings.csv"
    frame = pd.read_csv(path) if path.is_file() else pd.DataFrame()
    if not frame.empty:
        frame = frame[frame["config"].astype(str) != row["config"]]
    pd.concat([frame, pd.DataFrame([row])], ignore_index=True).to_csv(path, index=False)


def prediction_environment(args, config):
    env = dict(os.environ)
    if config.startswith("esmfold2"):
        weights = Path(args.esmfold2_weights).resolve()
        esmc = Path(args.esmc_weights).resolve()
        if not (weights / "model.safetensors").is_file():
            raise FileNotFoundError(f"ESMFold2 weights not found: {weights}")
        if not esmc.is_dir():
            raise FileNotFoundError(f"ESMC weights not found: {esmc}")
        env["BOLTZSCAN_ESMFOLD2_WEIGHTS"] = str(weights)
        env["BOLTZSCAN_ESMC_WEIGHTS"] = str(esmc)
        quiet = str(HERE / "quiet_runtime")
        env["PYTHONPATH"] = quiet + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    return env


def predict(args):
    run_dir = Path(args.run)
    _, manifest = load_run(run_dir)
    configs = tuple(args.configs or CONFIGS)
    failures = []
    for config in configs:
        input_dir = run_dir / "inputs" / config
        output_dir = run_dir / "outputs" / config
        output_dir.mkdir(parents=True, exist_ok=True)
        model = MODELS[config]
        command = [
            sys.executable, "-m", "boltzscan", "predict",
            "-i", str(input_dir), "-o", str(output_dir), "-m", model,
            "--seed", str(args.seed),
        ]
        started_at = datetime.now().astimezone()
        started = time.perf_counter()
        print(f"[experiment] start {config}: {' '.join(command)}")
        result = subprocess.run(command, env=prediction_environment(args, config))
        elapsed = time.perf_counter() - started
        generic_log = output_dir / "predict.log"
        if generic_log.is_file():
            shutil.copyfile(generic_log, run_dir / f"predict-{config}.log")
        completed = completed_predictions(predictions_dir(run_dir, config))
        status = "ok" if result.returncode == 0 else "failed"
        save_timing(run_dir, {
            "config": config,
            "model": model,
            "seed": args.seed,
            "start": started_at.isoformat(timespec="seconds"),
            "end": datetime.now().astimezone().isoformat(timespec="seconds"),
            "wall_seconds": round(elapsed, 6),
            "n_inputs": len(manifest),
            "n_completed": completed,
            "seconds_per_completed": round(elapsed / completed, 6) if completed else np.nan,
            "status": status,
        })
        print(f"[experiment] {config}: {completed}/{len(manifest)} in {elapsed:.2f}s ({status})")
        if result.returncode:
            failures.append(config)
    if failures:
        raise RuntimeError("Prediction failed or was incomplete for: " + ", ".join(failures))


def is_dna_chain(chain):
    residues = list(chain.get_residues())
    return bool(residues) and (
        sum(res.resname.strip().upper() in DNA_RESIDUES for res in residues) / len(residues)
    ) >= 0.8


def sequence_similarity(left, right):
    if not left or not right:
        return 0.0
    if left in right or right in left:
        return min(len(left), len(right)) / max(len(left), len(right))
    return SequenceMatcher(None, left, right).ratio()


def complementary_native_chain(native, selected_chain, dna_sequence):
    target = reverse_complement(dna_sequence)
    candidates = [
        (sequence_similarity(target, chain.sequence), chain_id)
        for chain_id, chain in native.child_dict.items()
        if chain_id != selected_chain and is_dna_chain(chain)
    ]
    if not candidates:
        return None, 0.0
    similarity, chain_id = max(candidates, key=lambda item: (item[0], item[1]))
    return (chain_id if similarity >= 0.75 else None), float(similarity)


def dockq_pair(model, native, native_left, native_right, model_left, model_right):
    from DockQ.DockQ import run_on_all_native_interfaces

    result, _ = run_on_all_native_interfaces(
        model,
        native,
        chain_map={native_left: model_left, native_right: model_right},
    )
    return float(next(iter(result.values()))["DockQ"]) if result else None


def artifacts(root, job_id):
    directory = Path(root) / job_id
    cif = next(iter(sorted(directory.glob("*.cif"))), None)
    confidence = next(iter(sorted(directory.glob("confidence_*.json"))), None)
    return cif, confidence


def confidence_metrics(path):
    if path is None:
        return {"global_iptm": None, "global_ptm": None}
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {"global_iptm": None, "global_ptm": None}
    return {"global_iptm": payload.get("iptm"), "global_ptm": payload.get("ptm")}


def ipsae_metrics(root, output, processes):
    from boltzscan.utils.ipsae_score import score_ipsae_table

    summary = score_ipsae_table(root, output=output, processes=processes)
    return pd.read_csv(summary.output).set_index("model_id").to_dict(orient="index")


def correlations(metrics):
    confidence_columns = [
        "global_iptm", "boltz_iptm", "ipsae", "ipsae_asym_min",
        "ipsae_iptm", "ipsae_iptm_asym_min", "pDockQ",
    ]
    rows = []
    for config, group in metrics.groupby("config", sort=False):
        for confidence in confidence_columns:
            pair = group[["dockq_interface", confidence]].apply(pd.to_numeric, errors="coerce").dropna()
            rows.append({
                "config": config,
                "confidence_metric": confidence,
                "n": len(pair),
                "pearson_r": pair["dockq_interface"].corr(pair[confidence]) if len(pair) >= 2 else np.nan,
                "spearman_rho": pair["dockq_interface"].rank().corr(pair[confidence].rank()) if len(pair) >= 2 else np.nan,
            })
    return pd.DataFrame(rows)


def model_summary(metrics, timings):
    scored_sets = {
        config: set(group.loc[pd.to_numeric(group["dockq_interface"], errors="coerce").notna(), "job_id"])
        for config, group in metrics.groupby("config", sort=False)
    }
    common = set.intersection(*scored_sets.values()) if scored_sets else set()
    rows = []
    for config, group in metrics.groupby("config", sort=False):
        dockq = pd.to_numeric(group["dockq_interface"], errors="coerce")
        common_dockq = pd.to_numeric(
            group.loc[group["job_id"].isin(common), "dockq_interface"], errors="coerce"
        )
        timing = timings[timings["config"].astype(str) == config]
        timing = timing.iloc[-1] if not timing.empty else {}
        rows.append({
            "config": config,
            "n_expected": len(group),
            "n_predicted": int(group["prediction_status"].eq("ok").sum()),
            "n_dockq": int(dockq.notna().sum()),
            "n_common": len(common_dockq),
            "success_rate": group["prediction_status"].eq("ok").mean(),
            "dockq_mean_all": dockq.mean(),
            "dockq_median_all": dockq.median(),
            "dockq_mean_common": common_dockq.mean(),
            "dockq_median_common": common_dockq.median(),
            "acceptable_or_better_common": (common_dockq >= 0.23).mean(),
            "medium_or_better_common": (common_dockq >= 0.49).mean(),
            "high_common": (common_dockq >= 0.80).mean(),
            "wall_seconds": timing.get("wall_seconds", np.nan),
            "seconds_per_completed": timing.get("seconds_per_completed", np.nan),
        })
    return pd.DataFrame(rows)


def markdown_table(frame):
    """Render a small DataFrame without adding pandas' optional tabulate dependency."""
    display = frame.copy()
    display = display.where(pd.notna(display), "")
    headers = [str(column) for column in display.columns]
    rows = [[str(value).replace("|", "\\|") for value in row] for row in display.itertuples(index=False, name=None)]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def write_benchmark_report(run_dir, manifest, summary, corr):
    candidates = summary.dropna(subset=["dockq_mean_common"])
    if candidates.empty:
        recommendation = "No paired common-case result is available yet."
    else:
        best = candidates.sort_values(
            ["dockq_mean_common", "success_rate", "seconds_per_completed"],
            ascending=[False, False, True],
            na_position="last",
        ).iloc[0]
        fastest = candidates.sort_values(
            ["seconds_per_completed", "dockq_mean_common"],
            ascending=[True, False],
            na_position="last",
        ).iloc[0]
        recommendation = (
            f"Quality-first candidate: `{best['config']}` with paired mean interface "
            f"DockQ {best['dockq_mean_common']:.3f}. Throughput-first candidate: "
            f"`{fastest['config']}` at {fastest['seconds_per_completed']:.3f} seconds per "
            "completed input. Mean DockQ and CASP quality classes are preferred over the "
            "median alone because this benchmark has a strongly zero-inflated DockQ distribution."
        )
    top_corr = corr.sort_values("spearman_rho", ascending=False).head(14)
    summary_table = markdown_table(summary.round(4))
    correlation_table = markdown_table(top_corr.round(4))
    report = f"""# BoltzScan FoldBench benchmark

- Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}
- Selected raw interfaces: {len(manifest)}
- Historical FoldBench predictions used: no
- Primary endpoint: native protein–DNA interface DockQ

## Model comparison

{summary_table}

## Confidence relationship with DockQ

{correlation_table}

## Interpretation

{recommendation}

Correlations are calculated within each model. The paired summary uses only
jobs scored successfully by every available configuration, preventing failed
long inputs from making one model look artificially better.
"""
    (Path(run_dir) / "benchmark_report.md").write_text(report)


def score(args):
    from DockQ.DockQ import load_PDB

    run_dir = Path(args.run)
    stored, manifest = load_run(run_dir)
    configs = tuple(args.configs or stored["configs"])
    metrics_dir = run_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    with StageLog(run_dir / "score.log", " ".join(sys.argv)):
        rows, native_cache = [], {}
        for config in configs:
            root = predictions_dir(run_dir, config)
            if not root.is_dir():
                print(f"Skip {config}: no predictions at {root}")
                continue
            ipsae = ipsae_metrics(root, metrics_dir / f"{config}_ipsae.csv", args.processes)
            for row in manifest.itertuples(index=False):
                cif, confidence = artifacts(root, row.job_id)
                result = {
                    "config": config,
                    "model": MODELS[config],
                    "job_id": row.job_id,
                    "dataset_row": row.dataset_row,
                    "complex_id": row.complex_id,
                    "native_protein_chain": row.native_protein_chain,
                    "native_dna_chain": row.native_dna_chain,
                    "protein_length": row.protein_length,
                    "dna_length": row.dna_length,
                    "prediction_status": "ok" if cif else "missing",
                    "prediction_structure": str(cif) if cif else None,
                    "dockq_interface": None,
                    "dockq_complement": None,
                    "dockq_tf_dna_mean": None,
                    "dockq_dna_duplex": None,
                    "native_complement_chain": None,
                    "complement_similarity": None,
                    "score_error": None,
                }
                result.update(confidence_metrics(confidence))
                result.update(ipsae.get(row.job_id, {}))
                if cif:
                    try:
                        native_file = str(row.native_structure)
                        if native_file not in native_cache:
                            native_cache[native_file] = load_PDB(native_file)
                        native = native_cache[native_file]
                        model = load_PDB(str(cif))
                        selected = dockq_pair(
                            model, native, row.native_protein_chain, row.native_dna_chain, "A", "B"
                        )
                        complement, similarity = complementary_native_chain(
                            native, row.native_dna_chain, row.dna_sequence
                        )
                        result["dockq_interface"] = selected
                        result["native_complement_chain"] = complement
                        result["complement_similarity"] = similarity
                        if complement:
                            complement_score = dockq_pair(
                                model, native, row.native_protein_chain, complement, "A", "C"
                            )
                            result["dockq_complement"] = complement_score
                            result["dockq_dna_duplex"] = dockq_pair(
                                model, native, row.native_dna_chain, complement, "B", "C"
                            )
                            values = [value for value in (selected, complement_score) if value is not None]
                            result["dockq_tf_dna_mean"] = np.mean(values) if values else None
                        else:
                            result["dockq_tf_dna_mean"] = selected
                    except Exception as exc:
                        result["score_error"] = str(exc)
                rows.append(result)

        metrics = pd.DataFrame(rows)
        if metrics.empty:
            raise ValueError("No fresh predictions were available to score")
        metrics.to_csv(metrics_dir / "metrics_per_case.csv", index=False)
        timing_path = run_dir / "timings.csv"
        timings = pd.read_csv(timing_path) if timing_path.is_file() else pd.DataFrame(columns=["config"])
        summary = model_summary(metrics, timings)
        corr = correlations(metrics)
        summary.to_csv(metrics_dir / "model_summary.csv", index=False)
        corr.to_csv(metrics_dir / "metric_correlations.csv", index=False)
        write_benchmark_report(run_dir, manifest, summary, corr)
        print(f"Scored {len(metrics)} model/case rows -> {metrics_dir}")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="stage", required=True)
    prepare_p = sub.add_parser("prepare", help="build fresh raw-data YAML inputs")
    prepare_p.add_argument("--data", default="data/foldbench")
    prepare_p.add_argument("--run", required=True)
    prepare_p.add_argument("--limit", type=int, default=None)

    crop_p = sub.add_parser(
        "prepare-crop",
        help="run hmmsearch and prepare DBD-union+flank FoldBench inputs",
    )
    crop_p.add_argument("--source-run", required=True)
    crop_p.add_argument("--run", required=True)
    crop_p.add_argument("--pfam", default="data/pfam/Pfam-A.hmm")
    crop_p.add_argument("--cpu", type=int, default=max(1, os.cpu_count() or 1))
    crop_p.add_argument("--flank", type=int, default=20)
    crop_p.add_argument(
        "--dbd-whitelist",
        help="TSV with a pfam_acc column; default uses production plant DBD_ACCS",
    )

    predict_p = sub.add_parser("predict", help="run current BoltzScan predict commands")
    predict_p.add_argument("--run", required=True)
    predict_p.add_argument("--configs", nargs="+", choices=CONFIGS, default=None)
    predict_p.add_argument("--seed", type=int, default=42)
    predict_p.add_argument("--esmfold2-weights", default=str(DEFAULT_ESMFOLD2_WEIGHTS))
    predict_p.add_argument("--esmc-weights", default=str(DEFAULT_ESMC_WEIGHTS))

    score_p = sub.add_parser("score", help="calculate DockQ/ipTM/ipSAE and correlations")
    score_p.add_argument("--run", required=True)
    score_p.add_argument("--configs", nargs="+", choices=CONFIGS, default=None)
    score_p.add_argument("--processes", type=int, default=max(1, os.cpu_count() or 1))
    return p


def main():
    args = parser().parse_args()
    try:
        {
            "prepare": prepare,
            "prepare-crop": prepare_crop,
            "predict": predict,
            "score": score,
        }[args.stage](args)
    except Exception as exc:
        raise SystemExit(f"FoldBench experiment {args.stage} failed: {exc}") from None


if __name__ == "__main__":
    main()
