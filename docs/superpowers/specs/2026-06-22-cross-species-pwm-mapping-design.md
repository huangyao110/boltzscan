# Cross-species PWM mapping via DBD homology

**Date:** 2026-06-22
**Status:** Design (approved in brainstorming; pending spec review)

## Goal

Build a species-specific motif/PWM database (the `data/pwms/OB_rose_pwms` layout:
`tf2pwms.json` + `txt/` + `meme/`) for an arbitrary plant species by transferring
known PWMs from cisBP and JASPAR onto that species' transcription factors via
**DNA-binding-domain (DBD) sequence similarity**.

First target: 野菊 / *Chrysanthemum* (gene IDs `CM01G…`), 2174 TFs already called
by `boltzscan find-tf` → `tasks/zyf/野菊基因组信息/tf_out/`. Output: `data/pwms/cm_pwms/`.

**Hard requirement: reusable.** Reference data (cisBP + JASPAR) is downloaded and
processed **once**, cached, and then any new species' TF set can be mapped against
it **without re-downloading**. This is achieved by splitting the work into two
stages with a stable on-disk reference store between them.

## Why this is needed (current state)

- `tair_tf2pwm.ipynb` only works for species already in cisBP (Arabidopsis): it
  reads `TF_Information_all_motifs_plus.txt` and groups `DBID → Motif_ID`. cisBP's
  own DBD inference is already baked in there.
- A novel genome (野菊) is **not** in cisBP, so cisBP's per-gene inference cannot
  be reused. The DBD-homology transfer must be done ourselves.
- The rose build (`data/pwms/OB_rose_pwms`, 969 TFs) shows the target structure but
  its build code was never committed. We build a proper, reusable pipeline instead.

## Method rationale

A TF's DNA-binding specificity is determined by its DBD; DBDs are conserved across
species (and across kingdoms for shared families). cisBP's published method infers
a TF's motif from a reference TF in the **same DBD family** whose DBD amino-acid
identity clears a **family-specific %ID threshold** (Weirauch et al. 2014). Species
is irrelevant — the threshold is the correctness gate. Therefore:

- Matching is keyed by **Pfam DBD accession**, not by family name. Plant-specific
  families (AP2/ERF PF00847, NAC PF02365, WRKY PF03106, GRAS, Dof, SBP, TCP, B3 …)
  simply have no animal/fungal references under the same Pfam, so cross-kingdom
  references are filtered out automatically and cannot mis-map.
- Cross-kingdom references (animals, fungi) are **included** and only contribute to
  shared families (Homeobox, bHLH, bZIP, MYB, C2H2, GATA, HSF, NF-Y, E2F/DP, MADS),
  and only when they clear the family threshold.

## Architecture — two stages

```
STAGE A  build-pwm-refs   (run once; network; cached)
  download cisBP (all species) + JASPAR CORE (all taxa)
    → reference protein sequences + motif PWMs + TF→motif metadata
  hmmscan reference proteins with our Pfam-A → extract reference DBD sequences
    → data/pwms/_refs/   (the reusable store)

STAGE B  map-pwm          (run per species; NO network; repeatable)
  input : species TF FASTA + its pfam.domtbl (from find-tf), or run hmmsearch
  extract species DBD sequences (same Pfam-A → comparable boundaries)
  per Pfam family: hmmalign species+ref DBDs to the Pfam HMM → pairwise %ID
  keep ref TFs with %ID ≥ family threshold → union their motif_ids
    → data/pwms/<species>_pwms/   ({tf2pwms.json, txt/, meme/, map_report.tsv})
```

Swapping species = re-run Stage B only.

## Components

New subpackage `boltzscan/pwmmap/`:

- `dbd.py` — **shared DBD extraction.** Given a protein FASTA + a Pfam domtblout
  (or run `hmmsearch --cut_ga` if none), return per-protein DBD subsequences tagged
  with `(pfam_acc, family, env_start, env_end, seq)`. Used by both stages so the two
  sides use identical domain boundaries. Reuses the `DBD_ACCS` / family map from
  `boltzscan/utils/find_tf.py`.
- `sources/cisbp.py` — download + parse cisBP **entire dataset** from static URLs
  (verified): `…/data/3_10/DataFiles/Bulk_downloads/EntireDataset/TF_Information.zip`
  (all species: `Motif_ID, Family_Name, DBDs, DBID, TF_Status, TF_Species`) and
  `…/PWMs.zip` (all M*.txt). Keep reference TFs with a real `Motif_ID` and
  `TF_Status == 'D'` (experimentally determined). cisBP carries **no protein
  sequences**, so sequences are resolved via UniProt (below). Yields
  `(ref_tf_id, source='cisbp', species, family, pfam_acc?, dbid, [motif_ids])`.
- `sources/jaspar.py` — download JASPAR CORE via the REST API
  (`/api/v1/matrix/?collection=CORE`, all `tax_group`s; verified: gives `pfm` +
  `uniprot_ids` + `family`). Yields the same tuple shape with `source='jaspar'` and
  `uniprot_ids`; converts each `pfm` to a PWM `.txt` + `.meme`.
- `sources/uniprot.py` — uniform protein-sequence resolver. JASPAR refs use their
  `uniprot_ids` directly (`https://rest.uniprot.org/uniprotkb/<acc>.fasta`); cisBP
  refs map `DBID → UniProt` via the UniProt ID-mapping API (batched), then fetch.
  Partial coverage is expected and recorded; refs with no sequence are skipped.
- `refs.py` — **Stage A driver.** Orchestrates downloads, runs one `hmmscan` over
  all reference proteins with `data/pfam/Pfam-A.hmm`, extracts reference DBDs, writes
  the reference store. Idempotent/resumable (skip present artifacts; `--refresh`).
- `align.py` — **blastp engine (standard + fast).** Build one BLAST DB from the
  reference DBD FASTA (`makeblastdb`), then `blastp` all species DBDs against it in a
  single multi-threaded run with tabular output
  (`-outfmt "6 qseqid sseqid pident length nident qlen slen qstart qend sstart send"`).
  This is fast (C implementation, all-vs-all in seconds) and gives the cisBP-standard
  DBD identity. %ID is computed **at the DBD level, not raw HSP pident**: keep the
  best HSP per (species DBD, reference DBD) pair, require query AND subject coverage
  ≥ `--min-cov` (default 0.8) of their DBD lengths, then `pct_id = nident / aln_len`.
  Hits are restricted to pairs sharing the same `pfam_acc` (encoded in the DBD IDs).
  Coverage gating prevents a short high-identity HSP from passing — this is what makes
  the number a real domain %ID rather than a local-patch artifact.
  `Bio.Align.PairwiseAligner` is **not** used (pure-Python all-vs-all is far too slow
  at reference scale). blastp comes from the `boltzscan` conda env (installed via
  bioconda) or an explicit `--blastp` path (e.g. the existing `…/envs/swt/bin/blastp`).
- `thresholds.py` — `family → DBD %ID cutoff` table (Weirauch 2014 values; default
  0.70 when a family is absent). `pfam_acc → family` via `find_tf.DBD_ACCS`.
- `mapper.py` — **Stage B driver.** Loads the reference store, extracts species DBDs,
  runs `align.py` per Pfam family, applies `thresholds.py`, writes `<species>_pwms/`.

CLI (wired into `boltzscan/cli.py`, thin handlers per repo convention):

- `boltzscan build-pwm-refs [--refs data/pwms/_refs] [--species all|plants|<list>]
  [--jaspar-collection CORE] [--pfam data/pfam/Pfam-A.hmm] [--refresh] [-c N]`
- `boltzscan map-pwm -f <species_tfs.fasta> -o data/pwms/<sp>_pwms
  [--domtbl <tf_out/pfam.domtbl>] [--refs data/pwms/_refs]
  [--threshold-mode family|global] [--threshold 0.70] [--min-cov 0.8]
  [--blastp <path>] [-c N]`

## On-disk layout

```
data/pwms/_refs/                      # Stage A output (reusable, built once)
  cisbp/<species>/…                   # raw cisBP archives (TF_Information, pwms, prot_seq)
  jaspar/…                            # raw JASPAR PFMs + metadata + fetched UniProt seqs
  ref_dbd.fasta                       # all reference DBD seqs; id = <ref_tf_id>__<pfam_acc>__<i>
  ref_index.tsv                       # ref_tf_id, source, species, family, pfam_acc, dbd_seq_id, motif_ids
  motif_store/txt/<motif_id>.txt      # every candidate PWM (cisBP M*.txt copied; JASPAR converted)
  motif_store/meme/<motif_id>.meme
  build_manifest.json                 # sources, versions, counts, timestamp

data/pwms/cm_pwms/                    # Stage B output (per species)
  tf2pwms.json                        # CM_gene -> [motif_id, …]   (OB_rose_pwms shape)
  txt/<motif_id>.txt                  # only matched motifs (copied from motif_store)
  meme/<motif_id>.meme
  map_report.tsv                      # CM_tf, family, pfam_acc, ref_tf, source, species, pct_id, motif_id
```

## Mapping algorithm (Stage B)

1. Extract species DBDs (`dbd.py`), grouped by `pfam_acc`.
2. For each `pfam_acc` present in the species set:
   a. Gather reference DBDs with the same `pfam_acc` from `ref_dbd.fasta`.
   b. `hmmalign` species+reference DBDs to that Pfam HMM (one alignment per family).
   c. For each species DBD, compute %ID to every reference DBD over shared match cols.
   d. Keep reference TFs with `%ID ≥ cutoff(family)`; union their `motif_ids`.
3. Write `tf2pwms.json` (species TFs with ≥1 hit), copy matched motifs from
   `motif_store` into `txt/` + `meme/`, write `map_report.tsv` audit trail.

A species TF with no same-Pfam reference above threshold gets no motif (recorded in
the report, omitted from `tf2pwms.json`) — matching OB_rose behaviour.

## Error handling

- Downloads: retry with backoff; resume by skipping present files; fail loudly only
  if a whole source is unreachable. `build_manifest.json` records what succeeded.
- Reference proteins lacking an extractable DBD are skipped (counted).
- JASPAR matrices whose UniProt sequence cannot be fetched are skipped (counted);
  their motif is dropped only if no sequence is obtainable.
- Stage B never touches the network; if `--refs` store is missing it errors with a
  pointer to run `build-pwm-refs` first.

## Testing

- Unit: blastp-tabular → coverage-gated DBD %ID (hand-built rows incl. a short
  high-identity HSP that must be rejected by `--min-cov`); threshold lookup;
  Pfam→family map; PFM→txt/meme conversion round-trip.
- Integration (small, offline): a mini reference store of ~20 DBDs (a few families)
  + ~10 species DBDs including a known Arabidopsis bHLH/Homeodomain → assert it maps
  to the expected cisBP motif at the family threshold; assert a plant-specific family
  (e.g. WRKY) gets no cross-kingdom hit.
- End-to-end smoke: map a 50-TF subset of `cm` against the real store; sanity-check
  per-family hit rates against rose's (`OB_rose_pwms`, ~14 motifs/TF average).

## Risks / open items

1. **Reference protein sequences (main coverage risk).** cisBP ships no sequences,
   so they come from UniProt: JASPAR via `uniprot_ids` (high coverage), cisBP via
   `DBID → UniProt` ID-mapping (partial — some cisBP source IDs won't map). The build
   records resolved/total per source; unresolved refs are simply dropped (their motif
   is still available to any other ref that maps). Restricting cisBP to
   `TF_Status == 'D'` keeps the fetch tractable and the references high-quality.
2. **All-species cisBP download size/time.** Large one-time cost; mitigated by being
   cached and resumable. `--species` lets us stage it (plants first, then the rest).
3. **DBD %ID definition.** Computed from blastp as `nident / aln_len` gated by
   ≥80% query+subject DBD coverage, so it reflects identity over the whole domain
   rather than a local HSP. Family-specific thresholds are calibrated on such
   domain-level identity, so this matches cisBP's intent.
4. **JASPAR ↔ UniProt** mapping gaps; cisBP already integrates much of JASPAR for
   plants, so JASPAR's marginal contribution is small (expected, accepted).

## Out of scope

- Building the FIMO/cisTarget DB from the resulting PWMs — that is the existing
  `boltzscan cistarg` step, which consumes `tf2pwms.json` + a meme dir downstream.
- De-novo motif discovery; SR-model scoring. We only transfer existing PWMs.
