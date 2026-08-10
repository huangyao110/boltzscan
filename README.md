# BoltzScan

BoltzScan prioritizes candidate TF–DNA interactions by combining sequence readout
(PWM transfer and FIMO) with shape readout (protein–DNA complex prediction and
ipSAE). Its primary command is one reproducible, run-directory based workflow:

```bash
boltzscan run \
  --tf TF.fa \
  --promoters promoters.fa \
  --run RUN \
  --refs data/pwms/_refs \
  --exclude-species 'Rosa chinensis' \
  --crop 20 \
  --model esmfold2
```

`--run` is the only output path supplied by the user. BoltzScan owns every
directory and filename below it. The command performs PWM transfer, promoter
scanning, input preparation, structure prediction, and ipSAE ranking. With
`--crop 20`, it first predicts one full-length TF–dsDNA complex per TF to locate
the DNA-contacting interface, then predicts all candidates using that interface
plus 20 amino acids of flank. Omit `--crop` to predict all candidates with the
full-length TF. Use `boltzscan run --help` for the complete interface.

The reusable PWM reference store defaults to `data/pwms/_refs`. Deployments
that keep shared data elsewhere can pass `--refs /path/to/pwm_refs`; this path
is read as a reference input and is not placed under or copied into `RUN`.

## Installation

BoltzScan supports Python 3.11.

```bash
pip install boltzscan
boltzscan doctor
boltzscan doctor --fix
```

`doctor` is read-only apart from its overwritten setup log. `doctor --fix`
creates a pinned BLAST+/HMMER/MEME Suite toolchain under the user's BoltzScan data
directory; it does not use `sudo`, change the active Conda environment, or edit
shell startup files. FIMO is supplied by that runtime toolchain; no Python
MEME wrapper is installed. Reference maintainers and strict LOSO runs can also
validate Tomtom with:

```bash
boltzscan doctor --fix --profile refs-builder
```

Override the managed location with `BOLTZSCAN_TOOL_DIR` or `doctor --tool-dir`.
BoltzScan resolves an explicitly supplied executable first, then its managed
toolchain, then `PATH`.

Model environments remain separate because their CUDA stacks are substantially
larger than the sequence-analysis tools. Most users should also install a
published, checksum-verified PWM reference release:

```bash
boltzscan install-pwm-refs
```

This downloads the built-in Google Drive release, verifies its fixed SHA256,
and installs it at `data/pwms/_refs`. If that directory already contains a
reference store, add `--replace`; the old store is moved to a timestamped
backup. Developers can override both `--url` and `--sha256` to test another
release.

Representative PWMs are the default everywhere. If `motif_clusters.tsv` is
missing, mapping/scanning stops with a command that fixes the problem. Use
`--no-pwm-cluster` only when an unclustered experiment is intentional.

Structure inference loads the vendored source trees at `boltzscan/boltz/src`
and `boltzscan/esm`; both are included in the wheel. The sibling `boltz` and
`esmfold2` Conda environments provide only their heavy dependencies, model
weights, and CUDA runtime. Standard Boltz1/2 enter the pinned upstream CLI with
native defaults. `boltz1_ode` and `boltz2_ode` are BoltzScan-owned adapters that
inject ODE diffusion parameters into the corresponding base family because the
upstream CLI does not provide these modes. Override the environment interpreters
with `BOLTZSCAN_BOLTZ_PYTHON` or `BOLTZSCAN_ESMFOLD_PYTHON` when they are stored
elsewhere.

ESMFold2 and ESMC-6B weights may also live outside the Hugging Face cache. Set
their directories before prediction; when the ESMFold2 directory contains
`ccd.pkl`, BoltzScan selects it automatically and records all three resolved
locations in `inference_parameters.json`:

```bash
export BOLTZSCAN_ESMFOLD2_WEIGHTS=/path/to/esm_weights/ESMFold2
export BOLTZSCAN_ESMC_WEIGHTS=/path/to/esm_weights/ESMC-6B
boltzscan predict -i RUN/inputs/full -o RUN -m esmfold2 --seed 42
```

Reference maintainers rebuild from the official CIS-BP, JASPAR, and UniProt
web sources. One command builds the DBD store, creates representative PWM
clusters, and packages the runtime-only release plus its SHA256 file:

```bash
boltzscan build-pwm-refs \
  --refs data/pwms/_refs \
  --archive dist/boltzscan_pwm_refs.tar.gz \
  --refresh
```

The runtime archive contains `ref_dbd.fasta`, `ref_index.tsv`,
`ref_proteins.fasta`, `motif_clusters.tsv`, `motif_store/{txt,meme}`, and build
and release manifests. It excludes raw download caches, logs, Pfam domtblout,
and Tomtom work files. Upload the archive and checksum separately, then update
the built-in release URL and SHA256 constants in
`boltzscan/pwmmap/archive.py`.

Before publishing a combined CIS-BP/JASPAR artifact publicly, confirm the
applicable redistribution terms and retain source attribution. JASPAR states
that its data are CC BY 4.0; the release manifest intentionally records all
source endpoints and this redistribution notice.

## Primary TF–DNA workflow

### Required inputs

`--tf` is a protein FASTA containing one or more TFs. `--promoters` is a
promoter FASTA whose record IDs identify the target genes. `--refs` points to
the reusable reference PWM store.

For a cross-species evaluation, pass the target species through
`--exclude-species`. BoltzScan removes that species from the reference store
*before* DBD-based PWM transfer and FIMO scanning. Thus a tomato or rose TF is
not scanned using a PWM directly annotated to that same species; the sequence
readout is transferred from the remaining reference species.

Use `--stage prepare` with the same arguments when only the PWM, FIMO,
candidate tables, and model inputs are needed. This is useful before committing
GPU time or ordering a Y1H panel. Prepare records `model`, `crop`, and `seed` in
`RUN/run_config.json`; later stages restore them automatically:

```bash
boltzscan run --run RUN --stage predict
boltzscan run --run RUN --stage score
```

### Full-length and interface-cropped modes

The two production modes are mutually exclusive:

- without `--crop`, every model-level candidate is predicted with the full TF;
- with `--crop 20`, BoltzScan selects one representative dsDNA per TF, performs
  one auxiliary full-length prediction per TF, finds protein residues within
  5 Å of either DNA strand, and predicts every candidate using the contact span
  plus 20 amino acids on each side.

The auxiliary prediction is only an interface localizer; it is not projected to
promoters or included in the final score table. HMMER DBD hits remain part of
PWM transfer but no longer determine the structural crop. Interface and crop
coordinates are 1-based and inclusive, and are recorded in
`inputs/interface_boundaries.csv` and `inputs/crop20/crop_manifest.csv`.
Boltz inputs receive a synchronously cropped A3M owned by the crop input
directory. If an interface cannot be found, the command fails explicitly
instead of silently falling back to the HMMER interval or the full protein.

Select the predictor directly with
`--model {esmfold2,boltz1,boltz2,boltz1_ode,boltz2_ode}`; there is no separate engine
argument. ESMFold2 uses `msa: empty`. ESMFold2, Boltz1, and Boltz2 retain their
native scientific inference defaults. `boltz1_ode` and `boltz2_ode` apply the
same BoltzScan ODE protocol—two sampling steps and step scale 1.0—to the
Boltz1 and Boltz2 model families respectively. The only public inference
override is `--seed`, which is accepted by every model; omitting it preserves
the model's native seed behavior. With an explicit seed, ESMFold2 also enables
deterministic PyTorch/CUDA execution so repeated runs through either the
vendored source API or BoltzScan produce identical coordinates and confidence
arrays on the same hardware/software stack. Cross-version or cross-GPU
bitwise identity is not promised.

The ODE adapter is part of BoltzScan rather than a hand-edited external Boltz
checkout. It maps `boltz1_ode` to the Boltz1 family and `boltz2_ode` to the
Boltz2 family, then replaces only that family's diffusion parameters. The base
checkpoint, feature pipeline, Pairformer, DataModule, writer, and precision
therefore remain model-family correct.

Before loading a model, `predict` validates every YAML as one protein plus a
reverse-complementary DNA duplex, rejects illegal polymer characters and
missing files, and checks that the first MSA sequence equals the YAML protein
sequence after A3M gap/insertion removal. It prints `MSA`, `no-MSA`, or `mixed`
for the input directory. Every Boltz model requires a matching local `.a3m` or
`.csv` MSA for each protein; only ESMFold2 permits explicit no-MSA input.

For any Boltz model, the primary `run` command generates or reuses these A3M
files automatically under `RUN/msa/`. The lower-level `msa` developer command
can also call the remote service directly; no Protenix installation is needed:

```bash
boltzscan msa --fasta tf_proteins.fasta --output msa
```

For each FASTA record this writes `msa/<protein_id>/0.a3m`. Existing matching
A3M files are reused. Use `--jobs` only
when concurrent server requests are appropriate.

The lower-level `fimo2yaml` command converts the de-duplicated FIMO model table
to full-length inputs:

```bash
boltzscan fimo2yaml \
  --fimo RUN/scan/filtered_model_scan_level_res.csv \
  --msa-dir RUN/msa \
  --output RUN/inputs
```

Interface cropping is a run-aware GPU stage because it must perform one
full-length inference per TF. The primary `run --crop` workflow invokes it
automatically. To inspect or resume it separately:

```bash
boltzscan run --run RUN --tf TF.fa --promoters promoters.fa \
  --model boltz2 --crop 20 --stage prepare
boltzscan find-interface --run RUN
boltzscan run --run RUN --stage predict
```

The last command detects the completed interface stage, reuses it, and predicts
only the generated `crop20` candidate inputs.

Boltz execution resources are internal. Preprocessing uses the CPU cores
available to the current process (Linux CPU affinity is respected). The
DataLoader uses two CPU worker processes to feed inference; this is not a GPU
count and does not increase GPU compute capacity. Keeping it at two avoids the
CUDA stalls and excess host-memory use often caused by large worker counts.

### Candidate rules and mapping tables

BoltzScan intentionally distinguishes the biological candidate table from the
structural-computation table.

At the **promoter level**, FIMO hits are de-duplicated only within the same
promoter and TF (including overlapping representations of the same site). If
identical dsDNA occurs in two different promoters, both promoter–TF edges are
retained. They remain distinct candidates for downstream Y1H validation.

At the **model-input level**, an exact dsDNA duplex and its reverse complement
are canonicalized only to reuse a structural calculation for the same TF. One
model input is predicted, then its result is projected back to every retained
promoter-level edge. This reduces redundant GPU work without silently removing
gene-level regulatory hypotheses.

The run directory contains the audit trail needed to recover either view:

```text
RUN/
├── pwm/                         representative PWM mapping and provenance
├── refs_loso/                   species-excluded PWM reference snapshot
├── scan/
│   ├── raw_fimo_scan_res.csv
│   ├── filtered_promoter_scan_level_res.csv
│   └── filtered_model_scan_level_res.csv
├── inputs/
│   ├── full/                    present only in full-length mode
│   ├── interface/               crop mode: one localizer YAML per TF
│   ├── interface_boundaries.csv crop mode: predicted contact intervals
│   └── crop20/                  crop mode: candidate YAMLs, A3Ms, manifest
├── run.log                      latest top-level invocation (overwritten)
├── map-pwm.log                  latest PWM mapping step
├── fimo-scan.log                latest promoter scan step
├── msa.log                      latest MSA step (Boltz models)
├── fimo2yaml.log                latest YAML preparation step
├── find-interface.log           crop mode: localizer inference and contacts
├── predict-full.log             full-length mode only
├── predict-crop20.log           crop mode only
├── wash.log                     latest output publication step
├── ipsae.log                    latest scoring step
├── <engine-native-results>/     untouched model output, including processed/
├── <method>_prediction/         wash-published prediction view
│   ├── full/                    full-length mode only
│   ├── interface/               crop mode auxiliary predictions
│   ├── crop20/                  crop mode candidate predictions
│   └── wash_manifest.json
└── results/                     model scores and promoter-level rankings
```

Logs are deliberately stored directly in `RUN/`, PLINK-style. Re-running a
command replaces its previous log instead of appending another session. Output
data are not silently deleted; use `--resume` when resuming a prepared run.

Inference always retains the model-native tree: Boltz writes
`boltz_results_<arm>/` (including `processed/`), while ESMFold2 writes
`esmfold_results_<arm>/predictions/`. A low-level prediction needs only the
prepared YAML directory, run directory, and model; a seed is optional:

```bash
boltzscan predict -i RUN/inputs/full -o RUN -m esmfold2
boltzscan predict -i RUN/inputs/full -o RUN -m boltz2 --seed 42
```

Publish the shared layout afterwards:

```bash
boltzscan wash --run RUN --model boltz2
boltzscan wash --run RUN --model boltz2 --mode hard
```

`soft` is the default and creates live relative directory symlinks. `hard`
creates real directories with hard-linked files, falling back to copies only
when the filesystem cannot hard-link; rerun it to refresh the snapshot. Both
modes preserve the native source tree. `wash_manifest.json` records every
source, destination, mode, file count, and available inference parameters.
The concrete published root is `esmfold2_prediction`, `boltz1_prediction`,
`boltz2_prediction`, `boltz1_ode_prediction`, or `boltz2_ode_prediction`.

`model_id` identifies a unique structural job. Do not use
`filtered_model_scan_level_res.csv` as a gene-regulatory edge list. Use
`filtered_promoter_scan_level_res.csv` or a final promoter-level ranking for
biological interpretation and Y1H selection; the promoter table already carries
the `model_id` needed to project one structural calculation back to every edge.

### Interpreting scores

For each predicted complex, the TF is chain `A` and the two DNA strands are
chains `B` and `C`. BoltzScan ranks the TF–DNA interface using the weaker of
the two symmetrized strand scores:

```text
ipsae = min(max(A→B, B→A), max(A→C, C→A))
```

The DNA–DNA interface is excluded. Scores are first calculated once per model
input and then joined through the `model_id` column in
`filtered_promoter_scan_level_res.csv`, preserving all promoter-level candidates
in every available structural-arm ranking.

The standard workflow writes one structural result pair for its active mode:
`full_model_input_ipsae.csv.gz` and `full_promoter_level_ipsae.csv.gz`, or the
corresponding `crop20_*` files. The interface-localizer predictions are never
treated as promoter candidates. Exact dsDNA on different promoters remains as
separate rows in the promoter-level table.

## Lower-level commands

`boltzscan run` is the recommended path for a new TF–DNA experiment. The
lower-level commands remain available when a stage must be inspected or reused:

```bash
boltzscan promoter --help          # extract promoter FASTA from genome + GFF
boltzscan find-tf --help           # identify plant TFs and DBDs
boltzscan build-pwm-refs --help    # build the reusable reference store
boltzscan install-pwm-refs --help  # install a verified prebuilt reference release
boltzscan map-pwm --help           # transfer reference PWMs to a TF FASTA
boltzscan fimo-scan --help         # write the three FIMO scan tables
boltzscan fimo2yaml --help         # model-level FIMO table to full-length YAMLs
boltzscan find-interface --help    # infer one interface per TF and crop all candidates
boltzscan predict --help           # predict a directory of prepared YAMLs
boltzscan wash --help              # publish native predictions as soft/hard views
boltzscan ipsae --help             # score existing prediction outputs
```

The standalone scorer takes a complete run directory. It fixes both ipSAE
cutoffs at 10, reuses existing calculations, detects the prediction model and
whichever `full` or `crop<N>` candidate arm exists, and selects CPU concurrency
automatically:

```bash
boltzscan ipsae --run RUN
boltzscan ipsae --run RUN --model boltz2  # only if RUN has multiple models
```

It first writes one score per structural `model_id`, then joins those scores to
every matching row in `filtered_promoter_scan_level_res.csv`. Final tables for
the available structural arms and model-specific rankings are written under
`RUN/results/`; model-level tables remain intermediate audit files. A production
run contains either `full` or `crop<N>` candidate predictions; the auxiliary
`interface` arm is ignored by scoring.

`hit2fasta` is an optional subset-export step. It reads the retained
promoter-level scan table, maps its `motif_id` values back to TF IDs through
the matching `tf2pwms.json`, and writes only those TF protein records:

```bash
boltzscan hit2fasta \
  --scan-table RUN/scan/filtered_promoter_scan_level_res.csv \
  --tf2pwms RUN/pwm/tf2pwms.json \
  --protein-fasta data/genome/sly/tf_proteins.fasta \
  --output RUN/candidate_tf_proteins.fasta
```

This is not required by `boltzscan run`, which already carries the TF FASTA
through the complete workflow. Do not use the model-level scan table here.

These commands are useful for inspecting individual stages, but they do not
themselves replace the run-level promoter-to-model mapping described above.
