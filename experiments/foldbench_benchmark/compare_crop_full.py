#!/usr/bin/env python
"""Paired FoldBench comparison of full-length and DBD-cropped predictions."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


CONFIG_ORDER = [
    "esmfold2_msa",
    "esmfold2_nomsa",
    "boltz1",
    "boltz1_ode",
    "boltz2",
    "boltz2_ode",
]


def paired_summary(frame: pd.DataFrame, stratum: str) -> dict:
    delta = frame["dockq_crop"] - frame["dockq_full"]
    nonzero = delta[delta != 0]
    pvalue = np.nan
    if len(nonzero):
        pvalue = wilcoxon(
            frame.loc[nonzero.index, "dockq_crop"],
            frame.loc[nonzero.index, "dockq_full"],
            zero_method="wilcox",
            alternative="two-sided",
        ).pvalue
    return {
        "config": frame["config"].iloc[0],
        "stratum": stratum,
        "n_pairs": len(frame),
        "dockq_full_mean": frame["dockq_full"].mean(),
        "dockq_crop_mean": frame["dockq_crop"].mean(),
        "dockq_mean_delta": delta.mean(),
        "dockq_median_delta": delta.median(),
        "crop_wins": int((delta > 0).sum()),
        "ties": int((delta == 0).sum()),
        "crop_losses": int((delta < 0).sum()),
        "wilcoxon_p": pvalue,
    }


def markdown_table(frame: pd.DataFrame) -> str:
    def render(value):
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{value:.4f}"
        return str(value)

    header = "| " + " | ".join(frame.columns) + " |"
    divider = "| " + " | ".join("---" for _ in frame.columns) + " |"
    rows = [
        "| " + " | ".join(render(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, divider, *rows])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-run", required=True)
    parser.add_argument("--crop-run", required=True)
    args = parser.parse_args()

    full_run = Path(args.full_run)
    crop_run = Path(args.crop_run)
    out_dir = crop_run / "metrics"
    out_dir.mkdir(parents=True, exist_ok=True)

    full = pd.read_csv(full_run / "metrics" / "metrics_per_case.csv")
    crop = pd.read_csv(crop_run / "metrics" / "metrics_per_case.csv")
    manifest = pd.read_csv(crop_run / "benchmark_manifest.csv")
    crop_state = manifest[[
        "job_id",
        "protein_length",
        "full_protein_length",
    ]].copy()
    crop_state["actually_shortened"] = (
        crop_state["protein_length"] < crop_state["full_protein_length"]
    )

    columns = ["config", "job_id", "dockq_interface"]
    paired = crop[columns].merge(
        full[columns],
        on=["config", "job_id"],
        how="inner",
        suffixes=("_crop", "_full"),
    ).rename(columns={
        "dockq_interface_crop": "dockq_crop",
        "dockq_interface_full": "dockq_full",
    }).merge(crop_state, on="job_id", how="left")
    paired = paired.dropna(subset=["dockq_crop", "dockq_full"]).copy()
    paired["dockq_delta"] = paired["dockq_crop"] - paired["dockq_full"]
    paired.to_csv(out_dir / "crop_vs_full_per_case.csv", index=False)

    summaries = []
    for config in CONFIG_ORDER:
        config_frame = paired[paired["config"] == config]
        if config_frame.empty:
            continue
        summaries.append(paired_summary(config_frame, "all"))
        for shortened, label in ((True, "shortened"), (False, "unchanged")):
            part = config_frame[config_frame["actually_shortened"] == shortened]
            if not part.empty:
                summaries.append(paired_summary(part, label))
    summary = pd.DataFrame(summaries)
    summary["config"] = pd.Categorical(
        summary["config"], categories=CONFIG_ORDER, ordered=True
    )
    summary = summary.sort_values(["config", "stratum"]).reset_index(drop=True)
    summary["config"] = summary["config"].astype(str)
    summary.to_csv(out_dir / "crop_vs_full_summary.csv", index=False)

    # Strict paired set: every configuration has a crop and full DockQ score.
    strict_jobs = None
    for config in CONFIG_ORDER:
        jobs = set(paired.loc[paired["config"] == config, "job_id"])
        strict_jobs = jobs if strict_jobs is None else strict_jobs & jobs
    strict = paired[paired["job_id"].isin(strict_jobs or set())]
    strict_summary = []
    for config in CONFIG_ORDER:
        part = strict[strict["config"] == config]
        if not part.empty:
            strict_summary.append(paired_summary(part, "six-model-strict"))
    strict_summary = pd.DataFrame(strict_summary)
    strict_summary.to_csv(out_dir / "crop_vs_full_six_model_strict.csv", index=False)

    full_time = pd.read_csv(full_run / "timings.csv").add_suffix("_full")
    crop_time = pd.read_csv(crop_run / "timings.csv").add_suffix("_crop")
    timing = crop_time.merge(
        full_time,
        left_on="config_crop",
        right_on="config_full",
        how="left",
    )
    timing["wall_time_delta_seconds"] = (
        timing["wall_seconds_crop"] - timing["wall_seconds_full"]
    )
    timing["wall_time_ratio_crop_over_full"] = (
        timing["wall_seconds_crop"] / timing["wall_seconds_full"]
    )
    timing.to_csv(out_dir / "crop_vs_full_timings.csv", index=False)

    correlations = pd.read_csv(out_dir / "metric_correlations.csv")
    ipsae = correlations[correlations["confidence_metric"] == "ipsae"][[
        "config", "n", "pearson_r", "spearman_rho"
    ]]
    all_rows = summary[summary["stratum"] == "all"].copy()
    shortened_rows = summary[summary["stratum"] == "shortened"].copy()
    unchanged_rows = summary[summary["stratum"] == "unchanged"].copy()
    timing_report = timing[[
        "config_crop",
        "n_inputs_crop",
        "n_completed_crop",
        "wall_seconds_crop",
        "seconds_per_completed_crop",
        "n_inputs_full",
        "n_completed_full",
        "wall_seconds_full",
        "seconds_per_completed_full",
    ]].rename(columns={"config_crop": "config"})
    report = f"""# FoldBench DBD+20 versus full-length comparison

- General-whitelist interfaces: {len(manifest)}
- Unique complexes: {manifest['complex_id'].nunique()}
- Unique native protein chains: {manifest[['complex_id', 'native_protein_chain']].drop_duplicates().shape[0]}
- Actually shortened: {int(crop_state['actually_shortened'].sum())}
- Unchanged because DBD+20 covered the full chain: {int((~crop_state['actually_shortened']).sum())}
- Six-model strict paired jobs: {len(strict_jobs or set())}
- Seed: 42

## Paired DockQ, all scoreable jobs

{markdown_table(all_rows)}

## Paired DockQ, six-model strict intersection

{markdown_table(strict_summary)}

## Actually shortened proteins only

{markdown_table(shortened_rows)}

## Proteins unchanged by DBD+20

{markdown_table(unchanged_rows)}

## Recorded wall times

{markdown_table(timing_report)}

## Exact three-arm ipSAE relationship with crop DockQ

The `ipsae` column is
`min(max(A->B, B->A), max(A->C, C->A))`, where A is protein and B/C are
the two DNA strands.

{markdown_table(ipsae)}

## Interpretation guardrails

The Wilcoxon test is paired by the exact FoldBench protein-chain/DNA-interface
job.  `shortened` and `unchanged` are reported separately because DBD+20 does
not alter every short protein.  ESMFold2-MSA wall time includes a cold read of
the 25.5 GB local weights; ESMFold2-noMSA immediately followed it and used the
operating-system page cache, so their total wall times are not a clean measure
of the MSA ablation alone.
"""
    (crop_run / "crop_vs_full_report.md").write_text(report)
    print(f"Wrote paired comparison -> {crop_run / 'crop_vs_full_report.md'}")


if __name__ == "__main__":
    main()
