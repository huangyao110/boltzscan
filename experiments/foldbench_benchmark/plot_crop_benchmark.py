#!/usr/bin/env python
"""Create the publication figure for the FoldBench DBD+20 benchmark.

Figure contract
---------------
Core conclusion: a general Pfam whitelist improves DBD coverage, but DBD+20
cropping reduces Boltz2 interface accuracy, favoring full-length Boltz2.
Archetype: asymmetric quantitative grid with a paired DockQ hero panel.
Backend: Python/matplotlib-seaborn only.
Final size: 183 x 169 mm; editable SVG/PDF plus 600-dpi TIFF.
Statistics: exact paired jobs, Wilcoxon signed-rank tests, no multiplicity
correction because the six pre-specified model contrasts are shown explicitly.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Mandatory editable-text settings must precede figure creation.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"

import matplotlib as mpl
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns


ORDER = [
    "esmfold2_msa",
    "esmfold2_nomsa",
    "boltz1",
    "boltz1_ode",
    "boltz2",
    "boltz2_ode",
]
LABELS = {
    "esmfold2_msa": "ESMFold2–MSA",
    "esmfold2_nomsa": "ESMFold2–noMSA",
    "boltz1": "Boltz1",
    "boltz1_ode": "Boltz1–ODE",
    "boltz2": "Boltz2",
    "boltz2_ode": "Boltz2–ODE",
}
SHORT_LABELS = {
    "esmfold2_msa": "ESM–MSA",
    "esmfold2_nomsa": "ESM–noMSA",
    "boltz1": "B1",
    "boltz1_ode": "B1–ODE",
    "boltz2": "B2",
    "boltz2_ode": "B2–ODE",
}
COLORS = {
    "esmfold2_msa": "#367F87",
    "esmfold2_nomsa": "#9BC8CC",
    "boltz1": "#7F4B78",
    "boltz1_ode": "#CFAFC8",
    "boltz2": "#234F88",
    "boltz2_ode": "#91ACD0",
}
DOWN = "#C84E49"
UP = "#378C56"
NEUTRAL = "#7A7A7A"
LIGHT = "#D7D7D7"
INK = "#262626"


def apply_style() -> None:
    mpl.rcParams.update({
        "font.size": 6.4,
        "axes.titlesize": 7.2,
        "axes.titleweight": "bold",
        "axes.labelsize": 6.5,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 2.7,
        "ytick.major.size": 2.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })
    sns.set_style("ticks")


def panel_label(ax, label: str, x: float = -0.14, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        fontweight="bold",
        color=INK,
    )


def p_label(value: float) -> str:
    if value < 0.001:
        return "P < 0.001"
    return f"P = {value:.3f}"


def prepare_source_data(full_run: Path, crop_run: Path, source_dir: Path):
    source_dir.mkdir(parents=True, exist_ok=True)
    full_manifest = pd.read_csv(full_run / "benchmark_manifest.csv")
    crop_manifest = pd.read_csv(crop_run / "benchmark_manifest.csv")

    plant_manifest_path = crop_run.parent / "crop20" / "benchmark_manifest.csv"
    if not plant_manifest_path.is_file():
        raise FileNotFoundError(
            "Plant-whitelist comparison manifest is required for panel a: "
            f"{plant_manifest_path}"
        )
    plant_manifest = pd.read_csv(plant_manifest_path)
    coverage = pd.DataFrame({
        "whitelist": ["Plant TF", "General protein–DNA"],
        "interfaces": [len(plant_manifest), len(crop_manifest)],
        "total_interfaces": [len(full_manifest), len(full_manifest)],
    })
    coverage["coverage_fraction"] = (
        coverage["interfaces"] / coverage["total_interfaces"]
    )

    strict = pd.read_csv(
        crop_run / "metrics" / "crop_vs_full_six_model_strict.csv"
    )
    strict["config"] = pd.Categorical(strict["config"], ORDER, ordered=True)
    strict = strict.sort_values("config").reset_index(drop=True)

    per_case = pd.read_csv(
        crop_run / "metrics" / "crop_vs_full_per_case.csv"
    )
    shortened = per_case[per_case["actually_shortened"]].copy()
    shortened["config"] = pd.Categorical(shortened["config"], ORDER, ordered=True)

    model_summary = pd.read_csv(crop_run / "metrics" / "model_summary.csv")
    speed = model_summary[[
        "config",
        "n_expected",
        "n_predicted",
        "dockq_mean_common",
        "seconds_per_completed",
    ]].copy()
    speed["config"] = pd.Categorical(speed["config"], ORDER, ordered=True)
    speed = speed.sort_values("config").reset_index(drop=True)

    correlations = pd.read_csv(
        crop_run / "metrics" / "metric_correlations.csv"
    )
    correlations = correlations[correlations["confidence_metric"] == "ipsae"].copy()
    correlations["config"] = pd.Categorical(
        correlations["config"], ORDER, ordered=True
    )
    correlations = correlations.sort_values("config").reset_index(drop=True)

    coverage.to_csv(source_dir / "panel_a_coverage.csv", index=False)
    strict.to_csv(source_dir / "panel_b_paired_dockq.csv", index=False)
    shortened.to_csv(source_dir / "panel_c_shortened_per_case.csv", index=False)
    speed.to_csv(source_dir / "panel_d_speed_quality.csv", index=False)
    correlations.to_csv(source_dir / "panel_e_ipsae_correlations.csv", index=False)
    return coverage, strict, shortened, speed, correlations


def draw_coverage(ax, coverage: pd.DataFrame) -> None:
    y = np.array([1, 0])
    values = coverage["interfaces"].to_numpy()
    totals = coverage["total_interfaces"].to_numpy()
    colors = ["#A8A8A8", COLORS["boltz2"]]
    ax.barh(y, totals, color="#EEEEEE", height=0.52, zorder=1)
    ax.barh(y, values, color=colors, height=0.52, zorder=2)
    for yi, val, total in zip(y, values, totals):
        ax.text(
            val + 7,
            yi,
            f"{val}/{total}  ({100 * val / total:.1f}%)",
            ha="left",
            va="center",
            fontsize=6.1,
            color=INK,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(["Plant TF", "General\nprotein–DNA"])
    ax.set_xlim(0, 320)
    ax.set_xticks([0, 100, 200, 300])
    ax.set_xlabel("FoldBench interfaces detected")
    ax.set_title("Whitelist coverage", loc="left", pad=3)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    panel_label(ax, "a")


def draw_paired_dockq(ax, strict: pd.DataFrame) -> None:
    y = np.arange(len(strict))[::-1]
    for yi, row in zip(y, strict.itertuples()):
        full = row.dockq_full_mean
        crop = row.dockq_crop_mean
        direction = UP if crop > full else DOWN
        ax.plot([full, crop], [yi, yi], color=direction, lw=1.35, alpha=0.78, zorder=1)
        ax.scatter(
            full,
            yi,
            s=24,
            marker="o",
            color=COLORS[str(row.config)],
            edgecolor="white",
            linewidth=0.45,
            zorder=3,
        )
        ax.scatter(
            crop,
            yi,
            s=28,
            marker="D",
            facecolor="white",
            edgecolor=COLORS[str(row.config)],
            linewidth=1.0,
            zorder=3,
        )
        delta = crop - full
        right = max(full, crop) + 0.012
        ax.text(
            right,
            yi,
            f"Δ {delta:+.3f}; {p_label(row.wilcoxon_p)}",
            ha="left",
            va="center",
            fontsize=5.45,
            color=direction,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS[str(c)] for c in strict["config"]])
    ax.set_xlim(0.20, 0.45)
    ax.set_xticks([0.20, 0.25, 0.30, 0.35, 0.40, 0.45])
    ax.set_xlabel("Mean interface DockQ")
    ax.set_title("Full-length preserves Boltz2 interface accuracy", loc="left", pad=3)
    ax.grid(axis="x", color="#E7E7E7", lw=0.5, zorder=0)
    ax.tick_params(axis="y", length=0)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=NEUTRAL,
               markeredgecolor="white", markersize=5, label="Full length"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white",
               markeredgecolor=NEUTRAL, markersize=4.5, label="DBD+20"),
    ]
    ax.legend(
        handles=handles,
        loc="lower right",
        bbox_to_anchor=(1.0, 1.005),
        ncol=2,
        fontsize=5.8,
        handletextpad=0.4,
        columnspacing=1.0,
        borderaxespad=0,
    )
    ax.text(
        0.995,
        -0.30,
        "n = 186 paired interfaces; two-sided Wilcoxon test",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.2,
        color=NEUTRAL,
    )
    panel_label(ax, "b", x=-0.09)


def draw_delta_distribution(ax, shortened: pd.DataFrame) -> None:
    plot = shortened.copy()
    plot["method"] = plot["config"].map(LABELS)
    labels = [LABELS[c] for c in ORDER]
    palette = {LABELS[c]: COLORS[c] for c in ORDER}
    sns.boxplot(
        data=plot,
        y="method",
        x="dockq_delta",
        hue="method",
        order=labels,
        hue_order=labels,
        palette=palette,
        legend=False,
        width=0.58,
        saturation=0.85,
        showfliers=False,
        linewidth=0.7,
        medianprops={"color": "white", "linewidth": 1.1},
        whiskerprops={"color": NEUTRAL, "linewidth": 0.65},
        capprops={"color": NEUTRAL, "linewidth": 0.65},
        ax=ax,
    )
    rng = np.random.default_rng(42)
    for yi, config in enumerate(ORDER):
        vals = plot.loc[plot["config"] == config, "dockq_delta"].to_numpy()
        jitter = rng.uniform(-0.19, 0.19, len(vals))
        ax.scatter(
            vals,
            yi + jitter,
            s=3.0,
            color=COLORS[config],
            alpha=0.22,
            linewidth=0,
            rasterized=True,
            zorder=0,
        )
        ax.text(
            0.405,
            yi,
            f"n={len(vals)}",
            ha="right",
            va="center",
            fontsize=5.2,
            color=NEUTRAL,
        )
    ax.axvline(0, color=INK, lw=0.8, ls=(0, (2.5, 2)), zorder=1)
    ax.set_xlim(-0.65, 0.43)
    ax.set_xticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4])
    ax.set_xlabel("ΔDockQ (DBD+20 − full length)")
    ax.set_ylabel("")
    ax.set_yticks(np.arange(len(ORDER)))
    ax.set_yticklabels([SHORT_LABELS[c] for c in ORDER])
    ax.set_title("Effect among actually shortened proteins", loc="left", pad=3)
    ax.grid(axis="x", color="#ECECEC", lw=0.45, zorder=-2)
    ax.tick_params(axis="y", length=0)
    ax.text(
        0.01,
        0.03,
        "← cropping reduces accuracy",
        transform=ax.transAxes,
        fontsize=5.3,
        color=DOWN,
        ha="left",
    )
    panel_label(ax, "c", x=-0.09)


def draw_speed_quality(ax, speed: pd.DataFrame) -> None:
    offsets = {
        "esmfold2_msa": (3, 5),
        "esmfold2_nomsa": (3, -9),
        "boltz1": (3, 5),
        "boltz1_ode": (3, -8),
        "boltz2": (4, 5),
        "boltz2_ode": (4, -9),
    }
    for row in speed.itertuples():
        config = str(row.config)
        ax.scatter(
            row.seconds_per_completed,
            row.dockq_mean_common,
            s=34 if config.startswith("boltz2") else 25,
            color=COLORS[config],
            edgecolor="white",
            linewidth=0.55,
            zorder=3,
        )
        dx, dy = offsets[config]
        ax.annotate(
            SHORT_LABELS[config],
            (row.seconds_per_completed, row.dockq_mean_common),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=5.0,
            color=COLORS[config],
            ha="left",
            va="center",
        )
    ax.set_xlim(0.8, 10.3)
    ax.set_xticks([2, 4, 6, 8, 10])
    ax.set_ylim(0.215, 0.292)
    ax.set_yticks([0.22, 0.24, 0.26, 0.28])
    ax.set_xlabel("Seconds per completed input")
    ax.set_ylabel("Mean DockQ")
    ax.set_title("Crop speed–quality trade-off", loc="left", pad=3)
    ax.grid(color="#ECECEC", lw=0.45, zorder=0)
    panel_label(ax, "d", x=-0.24)


def draw_correlations(ax, corr: pd.DataFrame) -> None:
    y = np.arange(len(corr))[::-1]
    for yi, row in zip(y, corr.itertuples()):
        config = str(row.config)
        lo = min(row.pearson_r, row.spearman_rho)
        hi = max(row.pearson_r, row.spearman_rho)
        ax.plot([lo, hi], [yi, yi], color=COLORS[config], lw=1.1, alpha=0.72)
        ax.scatter(
            row.pearson_r,
            yi,
            marker="o",
            s=22,
            color=COLORS[config],
            edgecolor="white",
            linewidth=0.45,
            zorder=3,
        )
        ax.scatter(
            row.spearman_rho,
            yi,
            marker="s",
            s=20,
            facecolor="white",
            edgecolor=COLORS[config],
            linewidth=0.9,
            zorder=3,
        )
        ax.text(
            0.718,
            yi,
            f"n={int(row.n)}",
            ha="right",
            va="center",
            fontsize=5.2,
            color=NEUTRAL,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS[str(c)] for c in corr["config"]])
    ax.set_xlim(0.48, 0.72)
    ax.set_xticks([0.50, 0.55, 0.60, 0.65, 0.70])
    ax.set_xlabel("Correlation with interface DockQ")
    ax.set_title("Three-arm ipSAE tracks model accuracy", loc="left", pad=3)
    ax.grid(axis="x", color="#ECECEC", lw=0.45, zorder=0)
    ax.tick_params(axis="y", length=0)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=NEUTRAL,
               markeredgecolor="white", markersize=5, label="Pearson r"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="white",
               markeredgecolor=NEUTRAL, markersize=4.5, label="Spearman ρ"),
    ]
    ax.legend(
        handles=handles,
        loc="lower right",
        bbox_to_anchor=(1, 1.005),
        ncol=2,
        fontsize=5.8,
        handletextpad=0.4,
        columnspacing=1.0,
        borderaxespad=0,
    )
    ax.text(
        0.995,
        -0.30,
        "ipSAE = min[max(A→B, B→A), max(A→C, C→A)]",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.2,
        color=NEUTRAL,
    )
    panel_label(ax, "e", x=-0.09)


def write_legend(output_dir: Path) -> None:
    text = """# Figure legend

**General DBD annotation improves coverage, whereas DBD+20 cropping reduces
Boltz2 interface accuracy.** **a,** FoldBench protein–DNA interfaces detected
by the original plant-TF Pfam whitelist and the general protein–DNA Pfam
whitelist. Grey backgrounds denote all 300 eligible interfaces. **b,** Mean
interface DockQ for full-length and DBD+20 predictions in the strict set of
186 jobs scoreable for all six models in both forms. P values are from
two-sided paired Wilcoxon signed-rank tests; no multiple-testing correction was
applied to the six pre-specified model comparisons. **c,** Per-interface DockQ
change for proteins actually shortened by DBD+20. Boxes show the median and
interquartile range, whiskers extend to 1.5× the interquartile range, and
points are individual protein-chain–DNA-interface jobs. Sample sizes vary
because prediction or native-interface scoring failed for some configurations.
**d,** Crop-run mean DockQ on the six-model common set versus observed wall
time per completed input. ESMFold2–MSA includes a cold weight read, whereas the
following no-MSA run used the operating-system page cache; those two timings
are not a controlled MSA speed comparison. **e,** Pearson and Spearman
correlations between the exact three-arm ipSAE score and interface DockQ.
All predictions used seed 42. MSA-enabled inputs consumed precomputed A3M
files; no MSA was generated during inference.
"""
    (output_dir / "figure_legend.md").write_text(text)


def write_qa_notes(output_dir: Path) -> None:
    text = """# Figure QA notes

- Backend: Python/matplotlib-seaborn exclusively.
- Final dimensions: 183 × 169 mm (double-column).
- Primary SVG retains editable text; PDF uses TrueType fonts.
- TIFF is 600 dpi; PNG is a 300-dpi review preview.
- Method-family colors are consistent across panels; red/green encode only
  directional loss/gain and are not the sole condition encoding.
- Panel b uses an exact paired six-model intersection (n=186).
- Panel c includes only proteins whose sequence was actually shortened.
- Center/spread in panel c: median, IQR and 1.5×IQR whiskers.
- Statistical test: two-sided paired Wilcoxon signed-rank test; six planned
  contrasts; no multiple-comparison correction.
- There are no image manipulations or representative-image selections.
- Each quantitative panel has a dedicated source-data CSV.
"""
    (output_dir / "figure_qa.md").write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-run", required=True)
    parser.add_argument("--crop-run", required=True)
    parser.add_argument("--output-dir")
    args = parser.parse_args()

    full_run = Path(args.full_run)
    crop_run = Path(args.crop_run)
    output_dir = Path(args.output_dir) if args.output_dir else crop_run / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    source_dir = output_dir / "source_data"
    coverage, strict, shortened, speed, correlations = prepare_source_data(
        full_run, crop_run, source_dir
    )

    apply_style()
    width_in = 183 / 25.4
    height_in = 169 / 25.4
    fig = plt.figure(figsize=(width_in, height_in))
    grid = fig.add_gridspec(
        3,
        6,
        height_ratios=[1.02, 1.42, 1.12],
        hspace=0.72,
        wspace=1.10,
        left=0.09,
        right=0.985,
        top=0.965,
        bottom=0.09,
    )
    ax_a = fig.add_subplot(grid[0, :2])
    ax_b = fig.add_subplot(grid[0, 2:])
    ax_c = fig.add_subplot(grid[1, :4])
    ax_d = fig.add_subplot(grid[1, 4:])
    ax_e = fig.add_subplot(grid[2, 1:5])

    draw_coverage(ax_a, coverage)
    draw_paired_dockq(ax_b, strict)
    draw_delta_distribution(ax_c, shortened)
    draw_speed_quality(ax_d, speed)
    draw_correlations(ax_e, correlations)

    base = output_dir / "foldbench_dbd_crop_benchmark"
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    fig.savefig(
        base.with_suffix(".png"),
        dpi=300,
    )
    plt.close(fig)
    write_legend(output_dir)
    write_qa_notes(output_dir)
    print(f"Wrote figure bundle -> {output_dir}")


if __name__ == "__main__":
    main()
