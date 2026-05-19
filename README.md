# BoltzScan

BoltzScan is a workflow for constructing gene regulatory networks based on Boltz1x_ODE (ODE, sampling_step=2) and motif scan.

## Overview

BoltzScan provides a comprehensive pipeline for gene regulatory network (GRN) inference, combining ODE-based modeling with motif scanning capabilities. The workflow is designed to help researchers understand transcriptional regulation mechanisms through systematic analysis of gene expression data and regulatory motifs.

![BoltzScan Workflow](png/boltzscan_pipe.jpg)

## Features

BoltzScan consists of three main modules:

1. **Species-specific Motif Dataset Construction**
   - Build comprehensive motif databases tailored to specific species
   - Integrate known transcription factor binding motifs
   - Support custom motif library creation

2. **Motif Scan**
   - Scan genomic sequences for transcription factor binding sites
   - Identify potential regulatory elements
   - Score and rank motif matches based on statistical significance

3. **Boltz1x_ODE Inference**
   - Infer gene regulatory networks using ordinary differential equations
   - Implement Boltz1x_ODE model with configurable sampling steps (default: 2)(which refer to proteinx-mini)
   - Generate network topology and interaction strengths

## Installation

### Requirements

- Python 3.11+
- Required Python packages (see `pyproject.toml`)
- Sufficient computational resources for large-scale motif scanning

### Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/boltzscan.git
cd boltzscan

# Install required packages
pip install -e .
pip install -e .[cistarg]  # optional FIMO/cisTarget dependencies
```

## IPSAE Ranking

Use the `ipsae` command to generate a structure-derived ranking table from
Boltz prediction results:

```bash
python bscan.py ipsae \
  -r <boltz_results_dir_or_predictions_dir> \
  -o <scored.csv> \
  -p 8
```

To merge the scores into an existing candidate table, pass that table with
`-s/--score-file`. The table must contain a column matching Boltz prediction
directory names; use `-i/--id-col` if the column is not named `boltz_name`,
`name`, or `id`.

For TF-DNA predictions, chain `A` is treated as the transcription factor and
chains `B`/`C` as the DNA duplex. The final ranking score first symmetrizes
each TF-DNA interface with the IPSAE `Type == "max"` row, then requires both
DNA strands to score well:

```text
ipsae = min(max(A->B, B->A), max(A->C, C->A))
```

The DNA-DNA `B-C` interface is excluded. Sort candidates by the output `ipsae`
column in descending order; the command writes the table in this order by
default. The stricter raw directional minimum is also written as
`ipsae_asym_min` for diagnostics. The output includes Boltz `pair_chains_iptm`
as `boltz_iptm`, IPSAE-script `ipTM_af` as `ipsae_iptm`, both computed with the
same TF-DNA pair rule, plus `iptm_diff`. The global Boltz confidence `iptm` is
kept separately as `boltz_iptm_global`.

To compare IPSAE-derived `ipTM_af` with Boltz JSON `pair_chains_iptm`, run:

```bash
python bscan.py valid -r <boltz_results_dir_or_predictions_dir> -o valid.csv -p 8
```
