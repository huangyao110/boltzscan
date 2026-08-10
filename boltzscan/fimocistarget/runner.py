"""The small, public FIMO scan workflow.

The scan deliberately publishes only three CSV tables: raw motif hits, retained
promoter-level TF candidates, and de-duplicated model-level inputs.  Auxiliary
score/ranking matrices are not required by the structural workflow.
"""
from dataclasses import dataclass
import json
import logging
from pathlib import Path
import shlex
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO

from boltzscan.toolchain import executable_version, resolve_executable
from boltzscan.utils.io_utils import calculate_fasta_background


RAW_FIMO_FILENAME = 'raw_fimo_scan_res.csv'
PROMOTER_SCAN_FILENAME = 'filtered_promoter_scan_level_res.csv'
MODEL_SCAN_FILENAME = 'filtered_model_scan_level_res.csv'
PWM_MAPPING_MANIFEST = 'pwm_mapping.json'
_MAX_STORED_SCORES = 500_000
_MOTIF_PSEUDOCOUNT = 0.1
_RAW_COLUMNS = (
    'sequence_name', 'start', 'stop', 'strand',
    'score', 'pvalue', 'qvalue', 'motif_id',
)
_MEME_HEADER = (
    'MEME version 5\n\n'
    'ALPHABET= ACGT\n\n'
    'strands: + -\n\n'
    'Background letter frequencies (from uniform background):\n'
    'A 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n'
)
logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class FimoScanSummary:
    output_dir: Path
    raw_fimo: Path
    promoter_level: Path
    model_level: Path
    n_raw_hits: int
    n_promoter_candidates: int
    n_model_inputs: int
    pwm_clustered: bool


class _FimoScanner:
    """Official MEME Suite FIMO adapter for the three-table scan workflow."""

    def __init__(self, fasta_file, motif_dir, output_dir, custom_bg,
                 pvalue_thresh, motif_filter):
        self.fimo = resolve_executable('fimo')
        if self.fimo is None:
            raise FileNotFoundError(
                'MEME Suite FIMO not found; run `boltzscan doctor --fix`'
            )
        self.fasta_file = Path(fasta_file)
        self.motif_dir = Path(motif_dir)
        self.output_dir = Path(output_dir)
        self.pvalue_thresh = pvalue_thresh
        self.motif_filter = set(motif_filter)
        raw_background = {base: float(custom_bg[base]) for base in 'ACGT'}
        at_frequency = (raw_background['A'] + raw_background['T']) / 2
        cg_frequency = (raw_background['C'] + raw_background['G']) / 2
        # FIMO applies this transformation internally for double-stranded DNA.
        # Store the effective model explicitly so the logged input is reproducible.
        self.background = {
            'A': at_frequency,
            'C': cg_frequency,
            'G': cg_frequency,
            'T': at_frequency,
        }
        if not self.fasta_file.is_file():
            raise FileNotFoundError(f"Promoter FASTA not found: {self.fasta_file}")
        if not self.motif_dir.is_dir():
            raise FileNotFoundError(f"Motif directory not found: {self.motif_dir}")

    def _motif_files(self):
        files = sorted(
            path for path in self.motif_dir.glob('*.meme')
            if path.stem in self.motif_filter
        )
        if not files:
            raise ValueError(
                f"No mapped .meme motifs found in {self.motif_dir}"
            )
        return files

    @staticmethod
    def _write_motif_database(motif_files, output):
        """Combine one-motif MEME files and normalize IDs to file stems."""
        written = []
        with Path(output).open('w') as handle:
            handle.write(_MEME_HEADER)
            for motif_path in motif_files:
                lines = motif_path.read_text().splitlines()
                motif_lines = [
                    index for index, line in enumerate(lines)
                    if line.startswith('MOTIF ')
                ]
                matrix_lines = [
                    index for index, line in enumerate(lines)
                    if line.startswith('letter-probability matrix:')
                ]
                if len(motif_lines) != 1 or not matrix_lines:
                    raise ValueError(
                        f'Invalid one-motif MEME file: {motif_path}; expected '
                        'exactly one MOTIF line and a letter-probability matrix'
                    )
                start = motif_lines[0]
                block = lines[start:]
                block[0] = f'MOTIF {motif_path.stem}'
                handle.write('\n'.join(block).rstrip() + '\n\n')
                written.append(motif_path.stem)
        if not written:
            raise ValueError('No valid MEME motifs remained after validation')
        return written

    @staticmethod
    def _read_fimo_tsv(path):
        try:
            hits = pd.read_csv(path, sep='\t', comment='#')
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=_RAW_COLUMNS)
        hits = hits.rename(columns={
            'p-value': 'pvalue',
            'q-value': 'qvalue',
        })
        required = {'sequence_name', 'start', 'stop', 'strand', 'score', 'pvalue', 'motif_id'}
        missing = sorted(required.difference(hits.columns))
        if missing:
            raise ValueError(
                f"Invalid FIMO TSV {path}: missing columns {', '.join(missing)}"
            )
        if 'qvalue' not in hits:
            hits['qvalue'] = pd.NA
        for column in ('start', 'stop'):
            hits[column] = pd.to_numeric(hits[column], errors='raise').astype(int)
        for column in ('score', 'pvalue', 'qvalue'):
            hits[column] = pd.to_numeric(hits[column], errors='coerce')
        return hits.loc[:, _RAW_COLUMNS]

    def run(self):
        if not any(SeqIO.parse(self.fasta_file, 'fasta')):
            raise ValueError(f"Promoter FASTA contains no sequences: {self.fasta_file}")

        motif_files = self._motif_files()
        with tempfile.TemporaryDirectory(prefix='boltzscan-fimo-') as temporary:
            temporary = Path(temporary)
            motif_database = temporary / 'mapped_motifs.meme'
            motif_ids = self._write_motif_database(motif_files, motif_database)
            background = temporary / 'promoter_background.txt'
            background.write_text(
                '\n'.join(
                    f'{base} {self.background[base]:.12g}'
                    for base in 'ACGT'
                ) + '\n'
            )
            fimo_output = temporary / 'fimo_output'
            command = [
                self.fimo,
                '--verbosity', '1',
                '--oc', str(fimo_output),
                '--bfile', str(background),
                '--motif-pseudo', str(_MOTIF_PSEUDOCOUNT),
                '--thresh', str(self.pvalue_thresh),
                '--max-stored-scores', str(_MAX_STORED_SCORES),
                '--no-pgc',
                str(motif_database),
                str(self.fasta_file),
            ]
            print(
                f'FIMO {executable_version("fimo", self.fimo)}: {self.fimo}'
            )
            print(
                f'FIMO configuration: motifs={len(motif_ids)}, '
                f'pvalue<{self.pvalue_thresh:g}, '
                f'motif_pseudocount={_MOTIF_PSEUDOCOUNT:g}, '
                f'max_stored_scores={_MAX_STORED_SCORES}, both_strands=yes'
            )
            print(
                'FIMO effective background: '
                + ' '.join(
                    f'{base}={self.background[base]:.4f}' for base in 'ACGT'
                )
            )
            logger.info(
                'FIMO scanning %d mapped motifs: %s',
                len(motif_ids), shlex.join(command),
            )
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode != 0:
                detail = (result.stderr or result.stdout or 'no diagnostic output').strip()
                raise RuntimeError(
                    f'FIMO failed with exit code {result.returncode}: {detail}'
                )
            fimo_tsv = fimo_output / 'fimo.tsv'
            if not fimo_tsv.is_file():
                raise RuntimeError(f'FIMO did not produce its expected table: {fimo_tsv}')
            hits = self._read_fimo_tsv(fimo_tsv)

        if hits.empty:
            raise ValueError(
                "FIMO found no hits; check the p-value threshold, motif files, "
                "and promoter sequences"
            )
        hits.to_csv(self.output_dir / RAW_FIMO_FILENAME, index=False)
        return hits


def _motif_filter_from_tf_fasta(tf_fasta, tf2pwms):
    """Read TF IDs from `tf_fasta` headers, then collect the union of motif
    IDs assigned to those TFs in `tf2pwms`. Returns a set of motif IDs."""
    tf_ids = [record.id for record in SeqIO.parse(tf_fasta, 'fasta')]
    motifs = set()
    missing = []
    for tf_id in tf_ids:
        hits = tf2pwms.get(tf_id)
        if hits:
            motifs.update(hits)
        else:
            missing.append(tf_id)
    if missing:
        print(
            f"Warning: {len(missing)}/{len(tf_ids)} TF(s) in {tf_fasta} have no "
            f"tf2pwms entry (first few: {missing[:5]})"
        )
    if not motifs:
        raise ValueError(
            f"No motifs matched any TF in {tf_fasta}; check that IDs in the FASTA "
            f"headers match keys in the tf2pwms JSON."
        )
    return motifs


def _require_clustered_pwm_mapping(tf2pwms_path, no_pwm_cluster=False):
    """Return whether the mapping is clustered, failing closed by default."""
    tf2pwms_path = Path(tf2pwms_path)
    manifest_path = tf2pwms_path.parent / PWM_MAPPING_MANIFEST
    if no_pwm_cluster:
        if not manifest_path.is_file():
            return False
        try:
            return json.loads(manifest_path.read_text()).get('pwm_clustered') is True
        except json.JSONDecodeError:
            return False

    hint = (
        "Rerun `boltzscan map-pwm` without `--no-pwm-cluster` using a reference "
        "store from `boltzscan install-pwm-refs` or `boltzscan build-pwm-refs`; "
        "otherwise pass `--no-pwm-cluster` to explicitly scan unclustered PWMs."
    )
    if not manifest_path.is_file():
        raise ValueError(
            f"Cannot verify representative PWMs: missing {manifest_path}. {hint}"
        )
    try:
        manifest = json.loads(manifest_path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid PWM mapping manifest: {manifest_path}") from exc
    if manifest.get('pwm_clustered') is not True:
        refs = manifest.get('reference_dir', '<reference_dir>')
        raise ValueError(
            "FIMO scan requires representative PWMs by default, but this mapping "
            f"is unclustered. Rerun `boltzscan map-pwm --refs {refs}` without "
            "`--no-pwm-cluster`, or explicitly pass `--no-pwm-cluster` to FIMO."
        )
    return True


def run_fimo_scan(
    target_fasta,
    motif_dir,
    tf2pwms_path,
    tf_fasta,
    output_dir,
    pvalue_thresh=1e-4,
    overlap_thresh=0.0,
    dna_flank=5,
    no_pwm_cluster=False,
    reuse=False,
):
    """Scan promoters and write exactly the three public FIMO CSV levels."""
    from boltzscan.utils.boltz_input import (
        build_promoter_candidates_from_fimo,
        collapse_model_inputs,
    )
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = out_dir / RAW_FIMO_FILENAME
    promoter_path = out_dir / PROMOTER_SCAN_FILENAME
    model_path = out_dir / MODEL_SCAN_FILENAME
    pwm_clustered = _require_clustered_pwm_mapping(
        tf2pwms_path,
        no_pwm_cluster=no_pwm_cluster,
    )

    with open(tf2pwms_path) as f:
        tf2pwms = json.load(f)

    motif_filter = _motif_filter_from_tf_fasta(tf_fasta, tf2pwms)
    print(
        f"FIMO motifs: {len(motif_filter)} for {len(tf2pwms)} mapped TF(s); "
        f"PWM mode: {'cluster representatives' if pwm_clustered else 'unclustered (explicit)'}"
    )

    if reuse and raw_path.is_file():
        print(f"Reusing raw FIMO results: {raw_path}")
        raw_hits = pd.read_csv(raw_path)
    else:
        freqs = calculate_fasta_background(target_fasta)
        custom_bg = dict(zip('ACGT', freqs))
        print(
            f"Background ACGT frequencies (from {target_fasta}): "
            f"A={custom_bg['A']:.4f} C={custom_bg['C']:.4f} "
            f"G={custom_bg['G']:.4f} T={custom_bg['T']:.4f}"
        )
        scanner = _FimoScanner(
            fasta_file=target_fasta,
            motif_dir=motif_dir,
            output_dir=out_dir,
            custom_bg=custom_bg,
            pvalue_thresh=pvalue_thresh,
            motif_filter=motif_filter,
        )
        raw_hits = scanner.run()

    promoter_candidates = build_promoter_candidates_from_fimo(
        fimo_csv=raw_path,
        tf_fasta=tf_fasta,
        tf2pwms_path=tf2pwms_path,
        promoter_fasta=target_fasta,
        pvalue_thresh=pvalue_thresh,
        overlap_thresh=overlap_thresh,
        dna_flank=dna_flank,
    )
    promoter_candidates, model_inputs, _ = collapse_model_inputs(promoter_candidates)

    promoter_export = promoter_candidates.drop(
        columns=['promoter_seq', 'tf_seq'], errors='ignore'
    )
    promoter_export.to_csv(promoter_path, index=False)
    model_inputs.to_csv(model_path, index=False)

    print(
        f"FIMO scan: {len(raw_hits)} raw hits -> {len(promoter_export)} promoter "
        f"candidates -> {len(model_inputs)} model inputs"
    )
    print(f"Wrote {raw_path}")
    print(f"Wrote {promoter_path}")
    print(f"Wrote {model_path}")
    return FimoScanSummary(
        output_dir=out_dir,
        raw_fimo=raw_path,
        promoter_level=promoter_path,
        model_level=model_path,
        n_raw_hits=len(raw_hits),
        n_promoter_candidates=len(promoter_export),
        n_model_inputs=len(model_inputs),
        pwm_clustered=pwm_clustered,
    )
