# FoldBench predictor benchmark

This directory contains experiment-only code for comparing the six BoltzScan
prediction configurations on protein–dsDNA interfaces:

- `esmfold2_msa`
- `esmfold2_nomsa`
- `boltz1`
- `boltz1_ode`
- `boltz2`
- `boltz2_ode`

It does not add a BoltzScan command or modify production defaults. Historical
`data/foldbench/boltz_input/` and `boltz_results_*` data are not used. Every
prediction is launched through the current `python -m boltzscan predict`
command with the same seed. MSA search is never run during this experiment:
all MSA-enabled configurations consume the same precomputed FoldBench
`msa/<sample>/0.a3m`, while `esmfold2_nomsa` explicitly uses `msa: empty`.

Run from the repository root in the `boltzscan` Conda environment:

```bash
python experiments/foldbench_benchmark/benchmark.py prepare \
  --data data/foldbench \
  --run experiments/foldbench_benchmark/runs/full

python experiments/foldbench_benchmark/benchmark.py predict \
  --run experiments/foldbench_benchmark/runs/full \
  --seed 42

python experiments/foldbench_benchmark/benchmark.py score \
  --run experiments/foldbench_benchmark/runs/full
```

For a length-stratified pilot, add `--limit 12` to `prepare`. The prediction
stage writes `predict-<configuration>.log` at the run root. The scoring stage
writes:

- `metrics/metrics_per_case.csv`
- `metrics/model_summary.csv`
- `metrics/metric_correlations.csv`
- `benchmark_report.md`
- `timings.csv`

The primary accuracy endpoint is DockQ for predicted chains A/B against the
native protein chain and DNA chain named by the raw FoldBench row. The report
also calculates a paired common-case summary, so a model is not rewarded for
failing on difficult inputs.

Each predictor receives an isolated `outputs/<configuration>/` root. This
prevents one Boltz configuration from inheriting another configuration's
parsed-feature cache. Reported wall time begins with a precomputed A3M and
includes BoltzScan input validation, local feature conversion, model loading,
inference, and artifact writing; it excludes remote MSA generation.

Two current raw rows contain inosine (`I`). Because production BoltzScan only
accepts A/C/G/T/N DNA, those rows are recorded in `excluded_interfaces.csv`
and are not silently converted.

## General protein--DNA Pfam whitelist and DBD+20 run

The production `find-tf` accession set is plant-TF-oriented, so it is too
narrow for a cross-kingdom protein--DNA benchmark. Build the experiment's
auditable general whitelist from the official GO ontology and Pfam2GO mapping,
then use exact Pfam accessions to select the union of DBD envelope coordinates:

```bash
python experiments/foldbench_benchmark/build_general_dbd_whitelist.py \
  --output experiments/foldbench_benchmark/runs/crop20_general/references/general_dbd_pfam.tsv \
  --cache-dir experiments/foldbench_benchmark/runs/crop20_general/references/cache

python experiments/foldbench_benchmark/benchmark.py prepare-crop \
  --source-run experiments/foldbench_benchmark/runs/full \
  --run experiments/foldbench_benchmark/runs/crop20_general \
  --pfam data/pfam/Pfam-A.hmm \
  --dbd-whitelist experiments/foldbench_benchmark/runs/crop20_general/references/general_dbd_pfam.tsv \
  --flank 20 --cpu 32

python experiments/foldbench_benchmark/benchmark.py predict \
  --run experiments/foldbench_benchmark/runs/crop20_general --seed 42

python experiments/foldbench_benchmark/benchmark.py score \
  --run experiments/foldbench_benchmark/runs/crop20_general --processes 12

python experiments/foldbench_benchmark/compare_crop_full.py \
  --full-run experiments/foldbench_benchmark/runs/full \
  --crop-run experiments/foldbench_benchmark/runs/crop20_general
```

The primary whitelist rule is a Pfam2GO mapping to `GO:0003677` (DNA binding)
or one of its `is_a` descendants. The builder also retains the production
plant-TF accessions and a small reviewed list of exact canonical DBD
accessions. It writes provenance, source checksums, and optional observed-hit
counts beside the TSV. A sequence without an allowed Pfam hit is excluded;
there is no hidden full-length fallback. Every retained A3M is cropped with
the same 1-based inclusive protein interval and its ungapped query is checked
against the cropped protein before prediction.
