"""CLI entry point for the boltzscan package.

Exposed as the `boltzscan` console script via pyproject.toml and as
`python -m boltzscan` via :mod:`boltzscan.__main__`.
"""
import argparse
import os
import sys
from pathlib import Path

from boltzscan import __version__
from boltzscan.fimocistarget.expro import (
    extract_promoter_regions,
)
from boltzscan.pwmmap.archive import (
    DEFAULT_REFERENCE_RELEASE,
    DEFAULT_REFERENCE_RELEASE_SHA256,
    DEFAULT_REFERENCE_RELEASE_URL,
)
from boltzscan.predict.runners import PREDICTION_MODELS, available_cpu_count
from boltzscan.utils.boltz_input import write_fimo_model_yamls

__author__ = 'HuangYao'
__date__ = '2026-08-08'

BOLTZSCAN_BANNER = r"""
####   ###  #     ##### #####  ####  ####  ###  #   #
#   # #   # #       #       # #     #     #   # ##  #
#   # #   # #       #      #  #     #     #   # # # #
####  #   # #       #     #    ###  #     ##### #  ##
#   # #   # #       #    #        # #     #   # #   #
#   # #   # #       #   #         # #     #   # #   #
####   ###  #####   #   ##### ####   #### #   # #   #
""".strip('\n')
_BANNER_ANSI_256 = (45, 39, 33, 69, 99, 135, 171)


def _terminal_banner(stream=None):
    """Return a gradient banner only when help is printed to a real terminal."""
    stream = stream or sys.stdout
    color_enabled = (
        hasattr(stream, 'isatty')
        and stream.isatty()
        and 'NO_COLOR' not in os.environ
        and os.environ.get('TERM', '') != 'dumb'
    )
    if not color_enabled:
        return BOLTZSCAN_BANNER
    return '\n'.join(
        f'\033[1;38;5;{color}m{line}\033[0m'
        for color, line in zip(_BANNER_ANSI_256, BOLTZSCAN_BANNER.splitlines())
    )


def _build_parser():
    parser = argparse.ArgumentParser(
        prog='boltzscan',
        description=(
            f'{_terminal_banner()}\n\n'
            'BoltzScan: TF-DsDNA candidate prioritization by PWM/FIMO and structural modeling.\n'
            'Start a new experiment with `boltzscan run --help`; '
            'lower-level commands remain available for inspection.'
        ),
        epilog=(
            f'Author : {__author__}\n'
            f'Date   : {__date__}\n'
            f'Version: {__version__}'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-v', '--version', action='version',
        version=f'boltzscan {__version__} (by {__author__}, {__date__})')
    sub_parsers = parser.add_subparsers(dest='command')
    sub_parsers.required = True

    p = sub_parsers.add_parser(
        'doctor',
        help='Check BoltzScan dependencies and optionally install managed external tools',
        description=(
            'Inspect Python packages, external bioinformatics programs, PWM references, '
            'and Pfam data. By default this command is read-only. Pass --fix to install '
            'a pinned BLAST+/HMMER/MEME Suite toolchain in a private BoltzScan directory without '
            'changing the active Python or Conda environment.'
        ),
        epilog=(
            'Examples:\n'
            '  boltzscan doctor\n'
            '  boltzscan doctor --fix\n\n'
            'The managed toolchain includes BLAST+, HMMER, FIMO, and Tomtom.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--fix', action='store_true',
        help='Install/update the external toolchain (never uses sudo)')
    p.add_argument('--tool-dir', default=None,
        help='Managed toolchain directory (default: user data dir; or BOLTZSCAN_TOOL_DIR)')
    p.add_argument('--refs', default='data/pwms/_refs',
        help='PWM reference directory to validate (default: data/pwms/_refs)')
    p.add_argument('--pfam', default=None,
        help='Override Pfam HMM; default: compact library from --refs or package')

    # Primary TF-DNA workflow. Lower-level commands remain available for
    # inspection, but a new experiment should start here.
    p = sub_parsers.add_parser(
        'run',
        help='Run the auditable TF-DNA workflow: PWM transfer -> FIMO -> model inputs -> prediction -> ipSAE',
        description=(
            'Production entry point. Provide one named --run directory; BoltzScan '
            'owns all output directory and file names below it. Initial prepare '
            'also needs --tf and --promoters; predict/score resumes need only --run.'
        ),
        epilog=(
            'Example:\n'
            '  boltzscan run --run tomato_run --tf TF.fa '
            '--promoters promoters.fa --model boltz2\n\n'
            'All step logs are written directly under tomato_run/ and are '
            'overwritten when the same step is run again.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--tf', default=None,
        help='TF protein FASTA (required only when prepare is selected)')
    p.add_argument('--promoters', default=None,
        help='Promoter FASTA; record IDs are target-gene IDs (required only for prepare)')
    p.add_argument('-r', '--run', required=True,
        help='The only output path: BoltzScan owns every directory and filename below it')
    p.add_argument('--refs', default='data/pwms/_refs',
        help='PWM reference directory (default: data/pwms/_refs)')
    p.add_argument('--pfam', default=None,
        help='Override the compact plant-TF Pfam HMM supplied by --refs')
    p.add_argument('--exclude-species', action='append', default=[],
        help='Species to remove from the PWM reference store before transfer; repeat or comma-separate values')
    p.add_argument('--no-pwm-cluster', action='store_true',
        help='Explicitly use all mapped PWMs; default requires cluster representatives')
    p.add_argument('--pvalue', type=float, default=1e-5,
        help='FIMO threshold used for both scanning and candidate preparation (default: 1e-5)')
    p.add_argument('--overlap-thresh', type=float, default=0.0,
        help='Promoter-local FIMO core-overlap redundancy threshold (default: 0.0)')
    p.add_argument('--dna-flank', type=int, default=5,
        help='Bases on each side of the FIMO core in the modeled dsDNA (default: 5)')
    p.add_argument('--crop', type=int, default=None,
        help='Enable interface localization and crop with N AA flank (disabled by default)')
    p.add_argument('-m', '--model', default=None, choices=PREDICTION_MODELS,
        help='Structure model (default: esmfold2 for a new run; restored from RUN later)')
    p.add_argument('--seed', type=int, default=None,
        help='Optional random seed for any prediction model (default: model native behavior)')
    p.add_argument('--stage', default='all',
        choices=['all', 'prepare', 'predict', 'score'],
        help='Workflow stage to run (default: all)')
    p.add_argument('--resume', dest='reuse', action='store_true',
        help='Resume prepare from compatible PWM and raw FIMO outputs already in RUN')
    p = sub_parsers.add_parser(
        'wash',
        help='Publish native Boltz/ESMFold2 outputs in the shared method-named layout',
    )
    p.add_argument('-r', '--run', required=True, help='Run directory containing native predictions')
    p.add_argument('-m', '--model', default='esmfold2', choices=PREDICTION_MODELS,
        help='Prediction model whose native outputs will be published (default: esmfold2)')
    p.add_argument('--mode', default='soft', choices=['soft', 'hard'],
        help='soft: live relative symlinks; hard: real dirs with hard links/copy fallback')
    p.add_argument('--arms', nargs='+', default=None,
        help='Optional native arms to publish, e.g. full crop20 (default: auto-discover)')

    # msa
    p = sub_parsers.add_parser(
        'msa',
        help='Generate protein MSAs with the remote Protenix-compatible server',
        description=(
            'Submit every protein in a FASTA file to the remote MSA server. '
            'Only <output>/<protein_id>/0.a3m is retained; Protenix is not required.'
        ),
        epilog=(
            'Example:\n'
            '  boltzscan msa -f tf_proteins.fasta -o msa\n\n'
            'Output:\n'
            '  msa/<protein_id>/0.a3m\n'
            '  msa.log'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('-f', '--fasta', dest='fasta_file', required=True,
        help='Input protein FASTA; multiple records are supported')
    p.add_argument('-o', '--output', dest='output_dir', required=True,
        help='Output directory containing one <protein_id>/0.a3m per record')
    p.add_argument('-j', '--jobs', dest='jobs', type=int, default=1,
        help='Concurrent server requests (default: 1)')
    p.add_argument('--server-url', default=None,
        help='MSA service base URL (default: Protenix server or BOLTZSCAN_MSA_SERVER_URL)')

    # promoter
    p = sub_parsers.add_parser(
        'promoter',
        help='Extract promoter regions from genome+GFF as BED, FASTA, or both')
    p.add_argument('-gff', '--gff', required=False,
        default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.60.gff3',
        help='Input GFF3 annotation file (can be .gz). Assumes ### gene separator or groups implicitly.')
    p.add_argument('-g', '--genome', required=False,
        default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.dna_sm.toplevel.fa',
        help='Input Genome FASTA file (can be .gz).')
    p.add_argument('-o', '--output', required=False,
        default='./Rosa_chinensis_promoter',
        help='Output prefix; suffixes .bed and/or .fasta are added automatically')
    p.add_argument('-f', '--format', default='both', choices=['bed', 'fasta', 'both'],
        help='Output promoter BED, FASTA, or both (default: both)')
    p.add_argument('-u', '--upstream', type=int, default=2000, help='Bases UPSTREAM of CDS anchor (will be uppercase).')
    p.add_argument('-d', '--downstream', type=int, default=200, help='Bases DOWNSTREAM of CDS anchor, including anchor position (will be lowercase).')

    # fimo-scan
    p = sub_parsers.add_parser(
        'fimo-scan',
        help='Scan promoters and write raw, promoter-level, and model-level CSVs',
    )
    p.add_argument('-f', '--fasta', '--promoters', dest='fasta', required=True,
        help='Target promoter FASTA')
    p.add_argument('-m', '--motif-dir', required=True,
        help='Directory of mapped per-motif .meme files')
    p.add_argument('-t', '--tf2pwms', required=True,
        help='JSON mapping TF_id -> [motif_id, ...] from map-pwm')
    p.add_argument('--tf-fasta', '--tf', dest='tf_fasta', required=True,
        help='TF protein FASTA; record IDs must match tf2pwms keys')
    p.add_argument('-o', '--output', required=True,
        help='Output directory for the three FIMO CSV levels')
    p.add_argument('--pvalue-thresh', type=float, default=1e-4, help='FIMO p-value threshold (default: 1e-4)')
    p.add_argument('--overlap-thresh', type=float, default=0.0,
        help='Promoter-local FIMO core-overlap threshold (default: 0.0)')
    p.add_argument('--dna-flank', type=int, default=5,
        help='Bases retained on each side of a FIMO core (default: 5)')
    p.add_argument('--no-pwm-cluster', action='store_true',
        help='Explicitly allow an unclustered PWM mapping')
    p.add_argument('--reuse', action='store_true',
        help='Reuse raw_fimo_scan_res.csv and rebuild the two filtered levels')

    # fimo2yaml
    p = sub_parsers.add_parser(
        'fimo2yaml',
        help='Convert model-level FIMO results into MSA-backed prediction YAMLs',
        description=(
            'Read filtered_model_scan_level_res.csv from fimo-scan and write '
            'full-length YAMLs with matching A3M files.'
        ),
        epilog=(
            'Example:\n'
            '  boltzscan fimo2yaml \\\n'
            '    --fimo RUN/scan/filtered_model_scan_level_res.csv \\\n'
            '    --msa-dir RUN/msa --output RUN/inputs'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('-f', '--fimo', required=True,
        help='filtered_model_scan_level_res.csv produced by fimo-scan')
    p.add_argument('-m', '--msa-dir', required=True,
        help='Directory containing one <TF_id>/0.a3m per TF')
    p.add_argument('-o', '--output', required=True,
        help='Input root containing full/')

    p = sub_parsers.add_parser(
        'find-interface',
        help='Find one reusable TF-DNA interface and write cropped prediction inputs',
        description=(
            'Run or reuse one full-length TF-dsDNA prediction per TF to find '
            'protein residues within 5 A of DNA, then write interface+flank YAMLs '
            'and synchronously cropped A3M files for every DNA candidate in RUN.'
        ),
        epilog=(
            'Example:\n'
            '  boltzscan find-interface --run tomato_run\n\n'
            'The model and flank are restored from run_config.json. Use --flank '
            'only for a manually assembled run.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('-r', '--run', required=True,
        help='Prepared BoltzScan run containing the scan table and localizer inputs')
    p.add_argument('--flank', type=int, default=None,
        help='AA flank around the predicted contact span (default: RUN configuration)')
    p.add_argument('--cutoff', type=float, default=5.0,
        help='Heavy-atom protein-DNA contact cutoff in angstrom (default: 5.0)')
    p.add_argument('-m', '--model', default=None, choices=PREDICTION_MODELS,
        help='Prediction model (default: RUN configuration)')
    p.add_argument('--seed', type=int, default=None,
        help='Prediction seed (default: RUN configuration)')

    # hit2fasta
    p = sub_parsers.add_parser(
        'hit2fasta',
        help='Extract TF proteins whose motifs have retained promoter hits',
        description=(
            'Extract a smaller TF protein FASTA from a FIMO scan table. '
            'Use filtered_promoter_scan_level_res.csv for the final retained '
            'promoter candidates; the table must contain motif_id.'
        ),
        epilog=(
            'Example:\n'
            '  boltzscan hit2fasta --scan-table '
            'RUN/scan/filtered_promoter_scan_level_res.csv \\\n'
            '    --tf2pwms RUN/pwm/tf2pwms.json \\\n'
            '    --protein-fasta data/genome/sly/tf_proteins.fasta \\\n'
            '    --output RUN/candidate_tf_proteins.fasta\n\n'
            'Output: one FASTA record per matched TF, plus a sibling .log file.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('-s', '--scan-table', dest='hit_statistics',
        metavar='SCAN_TABLE', required=True,
        help='Recommended: filtered_promoter_scan_level_res.csv from fimo-scan')
    p.add_argument('-j', '--tf2pwms', required=True,
        help='tf2pwms.json produced by the matching map-pwm run')
    p.add_argument('-p', '--protein-fasta', required=True,
        help='Source TF protein FASTA; record IDs must equal tf2pwms TF IDs')
    p.add_argument('-o', '--output', required=True,
        help='Output FASTA containing only TFs with retained promoter hits')
    p.add_argument('--motif-col', default='motif_id',
        help='Motif ID column in hit statistics (default: motif_id)')

    # find-tf
    p = sub_parsers.add_parser(
        'find-tf',
        help='Identify plant TFs in a protein FASTA (hmmsearch vs Pfam-A + PlantTFDB-style rules)',
    )
    p.add_argument('-f', '--proteins', required=True,
        help='Protein FASTA to scan (record IDs become gene_id)')
    p.add_argument('-o', '--output', required=True,
        help='Output directory (tf_annotation.tsv, tf_proteins.fasta, pfam.domtbl)')
    p.add_argument('-m', '--pfam', default=None,
        help='Override Pfam HMM (default: packaged compact plant-TF library)')
    p.add_argument('-c', '--cpu', type=int, default=8,
        help='Threads for hmmsearch (default: 8)')
    p.add_argument('--domtbl', default=None,
        help='Reuse an existing hmmsearch domtblout instead of running hmmsearch')
    p.add_argument('--force', action='store_true',
        help='Re-run hmmsearch even if <output>/pfam.domtbl already exists')

    # build-pwm-refs
    p = sub_parsers.add_parser('build-pwm-refs',
        help='Maintainer: download sources, build+cluster PWM refs, and package a release')
    p.add_argument('--refs', default='data/pwms/_refs', help='Reference store dir (default: data/pwms/_refs)')
    p.add_argument('-m', '--pfam', default=None,
        help='Full Pfam-A.hmm maintainer input; a compact runtime subset is generated')
    p.add_argument('-c', '--cpu', type=int, default=8)
    p.add_argument('--no-cisbp', action='store_true')
    p.add_argument('--no-jaspar', action='store_true')
    p.add_argument('--refresh', action='store_true', help='Re-download even if present')
    p.add_argument('--qthresh', type=float, default=0.05,
        help='Tomtom q-value cutoff for representative PWM clustering (default: 0.05)')
    p.add_argument('--tomtom', default=None,
        help='Path to tomtom (default: PATH or the meme conda env)')
    p.add_argument('--archive', default=None,
        help='Release tar.gz path (default: <refs-parent>/boltzscan_pwm_refs.tar.gz)')
    p.add_argument('--no-archive', action='store_true',
        help='Build and cluster the store without creating a release archive')

    # install-pwm-refs
    p = sub_parsers.add_parser(
        'install-pwm-refs',
        help='User: download, SHA256-verify, and install a prebuilt PWM reference release',
    )
    p.add_argument('--url', default=DEFAULT_REFERENCE_RELEASE_URL,
        help=(
            'Release URL or local tar.gz path '
            f'(default: built-in {DEFAULT_REFERENCE_RELEASE} GitHub Release)'
        ))
    p.add_argument('--sha256', default=DEFAULT_REFERENCE_RELEASE_SHA256,
        help=(
            'Expected archive SHA256 '
            f'(default: built-in {DEFAULT_REFERENCE_RELEASE} checksum)'
        ))
    p.add_argument('--refs', default='data/pwms/_refs',
        help='Installation directory (default: data/pwms/_refs)')
    p.add_argument('--replace', action='store_true',
        help='Move an existing store to _refs.backup_<timestamp> before installing')

    # map-pwm
    p = sub_parsers.add_parser(
        'map-pwm',
        help='Stage B: map a species TF FASTA to reference motifs via DBD %%ID (no network)',
        description=(
            'Transfer reference PWMs to TFs by Pfam DBD identity. Processing order: '
            'optional LOSO filtering -> DBD BLAST -> identity threshold -> optional '
            'representative-PWM collapse. Defaults use family-specific thresholds and '
            'cluster representatives.'
        ),
        epilog='''
Examples:
  # Default: all references, family thresholds, representative PWMs
  boltzscan map-pwm -f TF.fa -o OUT

  # Reuse pfam.domtbl from find-tf instead of running hmmsearch again
  boltzscan map-pwm -f TF.fa -o OUT --domtbl pfam.domtbl

  # Keep every transferred original PWM
  boltzscan map-pwm -f TF.fa -o OUT --no-pwm-cluster

  # LOSO before BLAST; recluster retained motifs and output representatives
  boltzscan map-pwm -f TF.fa -o OUT --exclude-species "Solanum lycopersicum"

  # LOSO before BLAST; output all retained original PWMs
  boltzscan map-pwm -f TF.fa -o OUT --exclude-species "Solanum lycopersicum" --no-pwm-cluster

  # Require at least 80% DBD identity for every TF family
  boltzscan map-pwm -f TF.fa -o OUT --threshold-mode global --threshold 0.80

Main outputs:
  OUT/tf2pwms.json       TF -> scan-ready motif IDs
  OUT/map_report.tsv     best reference evidence per transferred motif
  OUT/pwm_mapping.json   mapping provenance and cluster mode
  OUT/{txt,meme}/        copied PWM files
  OUT/refs_loso/         filtered reference subset, only with --exclude-species
''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('-f', '--proteins', required=True, help='Species TF protein FASTA')
    p.add_argument('-o', '--output', required=True, help='Output dir, e.g. data/pwms/cm_pwms')
    p.add_argument('--refs', default='data/pwms/_refs')
    p.add_argument('--exclude-species', action='append', default=[],
        help='LOSO species removed before BLAST; repeat for multiple species')
    p.add_argument('--domtbl', default=None, help='Reuse a find-tf pfam.domtbl (else run hmmsearch)')
    p.add_argument('--threshold-mode', default='family', choices=['family', 'global'])
    p.add_argument('--threshold', type=float, default=0.70, help='Global %%ID cutoff (global mode)')
    p.add_argument('--min-cov', type=float, default=0.8, help='Min DBD coverage for a blast hit')
    p.add_argument('--blastp', default=None)
    p.add_argument('--makeblastdb', default=None)
    p.add_argument('-m', '--pfam', default=None,
        help='Override the compact plant-TF Pfam HMM supplied by --refs')
    p.add_argument('-c', '--cpu', type=int, default=8)
    p.add_argument('--no-pwm-cluster', action='store_true',
        help='Explicitly map all unclustered PWMs; default uses representatives')

    # predict
    p = sub_parsers.add_parser('predict',
        help='Predict protein-DNA complexes from a YAML input directory',
        description=(
            'Predict every *.yaml job in one prepared input directory. Choose the '
            'model directly; ESMFold2/Boltz1/Boltz2 keep native scientific defaults, '
            'while the Boltz1/2 ODE variants use the built-in BoltzScan protocol. '
            'A preflight validates TF-dsDNA sequences and MSA query consistency.'
        ),
        epilog=(
            'Examples:\n'
            '  boltzscan predict -i RUN/inputs/full -o RUN -m esmfold2\n'
            '  boltzscan predict -i RUN/inputs/crop20 -o RUN -m boltz2 --seed 42\n\n'
            'All Boltz models require a matching local MSA for every protein; '
            'ESMFold2 also permits explicit no-MSA input.\n'
            'Boltz resources are automatic: preprocessing uses available CPU cores; '
            'the inference DataLoader uses 2 CPU workers.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('-i', '--input-dir', required=True, help='Directory of Boltz-style input YAMLs')
    p.add_argument('-o', '--output', required=True,
        help='Run root; inference keeps the engine-native output layout')
    p.add_argument('-m', '--model', default='esmfold2', choices=PREDICTION_MODELS,
        help='Prediction model (default: esmfold2)')
    p.add_argument('--seed', type=int, default=None,
        help='Optional random seed (default: model native behavior)')

    # ipsae
    p = sub_parsers.add_parser(
        'ipsae',
        help='Score a completed run and restore promoter-level results',
        description=(
            'Calculate model-level TF-DNA ipSAE scores, then join them through '
            'model_id to every row in filtered_promoter_scan_level_res.csv. '
            'Chain A is the TF and chains B/C are the DNA duplex. Cutoffs, '
            'available full/crop arms, cached results, and CPU concurrency are automatic.'
        ),
        epilog=(
            'Examples:\n'
            '  boltzscan ipsae --run RUN\n'
            '  boltzscan ipsae --run RUN --model boltz2\n\n'
            'Usually --run is sufficient. Use --model only when RUN contains '
            'prediction outputs from more than one model. Final files are written '
            'under RUN/results/ at both model and promoter levels.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('-r', '--run', required=True,
        help='Run directory containing scan/, <model>_prediction/, and results/')
    p.add_argument('-m', '--model', choices=PREDICTION_MODELS, default=None,
        help='Prediction model (default: auto-detect the only model in RUN)')

    # txt2meme
    p = sub_parsers.add_parser(
        'txt2meme',
        help='Convert PWM txt files to MEME motif files',
    )
    p.add_argument('-i', '--input', required=True,
        help='Input PWM txt file or directory containing .txt files')
    p.add_argument('-o', '--output', required=True,
        help='Output .meme file for single input, or output directory for batch input')
    p.add_argument('-n', '--nsites', type=int, default=100,
        help='MEME nsites value to write (default: 100)')
    p.add_argument('-b', '--background', default=None,
        help='Comma-separated A,C,G,T background frequencies (default: 0.25,0.25,0.25,0.25)')
    p.add_argument('--bg-fasta', default=None,
        help='FASTA file used to calculate A,C,G,T background frequencies')
    p.add_argument('-f', '--force', action='store_true',
        help='Overwrite existing .meme files')

    # valid
    p = sub_parsers.add_parser(
        'valid',
        help='Validate IPSAE ipTM_af values against Boltz pair_chains_iptm values',
    )
    p.add_argument('-r', '--res-dir', required=True,
        help='Boltz result directory or its predictions/ subdirectory')
    p.add_argument('-o', '--output', default=None,
        help='Optional output CSV path for the validation report')
    p.add_argument('-t', '--tolerance', type=float, default=0.05,
        help='Allowed absolute difference between IPSAE ipTM_af and Boltz pair iptm (default: 0.05)')
    p.add_argument('-f', '--force', action='store_true',
        help='Recalculate IPSAE even when an existing *_10_10.txt file is present')
    p.add_argument('-p', '--processes', type=int, default=1,
        help='Number of concurrent IPSAE calculations (default: 1)')

    return parser


def _cmd_msa(args):
    from boltzscan.msa import DEFAULT_MSA_SERVER_URL, run_msa

    server_url = args.server_url or DEFAULT_MSA_SERVER_URL
    try:
        summary = run_msa(
            args.fasta_file,
            args.output_dir,
            jobs=args.jobs,
            server_url=server_url,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(f"boltzscan msa: {exc}") from None
    print(
        f"MSA complete: {summary.requested} proteins; "
        f"generated={summary.completed}, reused={summary.skipped}; "
        f"output={args.output_dir}"
    )


def _cmd_doctor(args):
    from boltzscan.doctor import run_doctor

    try:
        summary = run_doctor(
            fix=args.fix,
            refs=args.refs,
            pfam=args.pfam,
            tool_dir=args.tool_dir,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        raise SystemExit(f'boltzscan doctor: {exc}') from None
    if summary.failed_required:
        raise SystemExit(1)


def _cmd_run(args):
    from boltzscan.tfdna import run_tf_dna_workflow

    try:
        summary = run_tf_dna_workflow(
            tf_fasta=args.tf,
            promoters=args.promoters,
            out_dir=args.run,
            refs=args.refs,
            pfam=args.pfam,
            exclude_species=args.exclude_species,
            collapse_clusters=not args.no_pwm_cluster,
            pvalue=args.pvalue,
            overlap_thresh=args.overlap_thresh,
            dna_flank=args.dna_flank,
            crop=args.crop,
            model=args.model,
            seed=args.seed,
            stage=args.stage,
            reuse=args.reuse,
        )
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(f"boltzscan run: {exc}") from None
    print(f"Run directory: {summary.out_dir}")


def _cmd_wash(args):
    from boltzscan.predict.wash import wash_prediction_outputs

    try:
        summary = wash_prediction_outputs(
            args.run,
            model=args.model,
            mode=args.mode,
            arms=args.arms,
        )
    except (FileExistsError, FileNotFoundError, OSError, RuntimeError, ValueError) as exc:
        raise SystemExit(f"boltzscan wash: {exc}") from None
    print(
        f"Washed {summary.files} files from {', '.join(summary.arms)} "
        f"with mode={summary.mode} -> {summary.method_root}"
    )
    print(f"Tracking manifest: {summary.manifest}")


def _cmd_fimo_scan(args):
    from boltzscan.fimocistarget.runner import run_fimo_scan

    try:
        run_fimo_scan(
            target_fasta=args.fasta,
            motif_dir=args.motif_dir,
            tf2pwms_path=args.tf2pwms,
            tf_fasta=args.tf_fasta,
            output_dir=args.output,
            pvalue_thresh=args.pvalue_thresh,
            overlap_thresh=args.overlap_thresh,
            dna_flank=args.dna_flank,
            no_pwm_cluster=args.no_pwm_cluster,
            reuse=args.reuse,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(f'boltzscan fimo-scan: {exc}') from None


def _cmd_fimo2yaml(args):
    try:
        full_dir, _crop_summary = write_fimo_model_yamls(
            model_table=args.fimo,
            msa_dir=args.msa_dir,
            domtbl_path=None,
            output_dir=args.output,
            crop=None,
        )
    except (FileExistsError, FileNotFoundError, ValueError) as exc:
        raise SystemExit(f'boltzscan fimo2yaml: {exc}') from None
    print(f'Wrote full-length model YAMLs: {full_dir}')


def _cmd_find_interface(args):
    from boltzscan.interface import run_interface_stage

    try:
        summary = run_interface_stage(
            args.run,
            flank=args.flank,
            cutoff=args.cutoff,
            model=args.model,
            seed=args.seed,
        )
    except (FileNotFoundError, OSError, RuntimeError, ValueError) as exc:
        raise SystemExit(f'boltzscan find-interface: {exc}') from None
    print(
        f'Found reusable interfaces for {summary.n_tfs} TF(s); wrote '
        f'{summary.n_model_inputs} cropped inputs -> {summary.crop_input_dir}'
    )
    print(f'Interface coordinates: {summary.boundaries}')


def _cmd_candidate_tf_fasta(args):
    from boltzscan.utils.tf_candidates import build_candidate_tf_fasta

    try:
        summary = build_candidate_tf_fasta(
            hit_statistics=args.hit_statistics,
            tf2pwms_path=args.tf2pwms,
            protein_fasta=args.protein_fasta,
            output=args.output,
            motif_col=args.motif_col,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f'boltzscan hit2fasta: {exc}') from None
    print(
        f"Wrote {summary.output} "
        f"({summary.written_tfs}/{summary.candidate_tfs} candidate TFs with protein sequence)"
    )
    print(
        f"Motifs: {summary.matched_motifs}/{summary.hit_motifs} hit motifs matched tf2pwms"
    )
    if summary.missing_tfs:
        print(
            f"Warning: {len(summary.missing_tfs)} candidate TFs missing from protein FASTA "
            f"(first few: {summary.missing_tfs[:5]})",
            file=sys.stderr,
        )
    if summary.unmatched_motifs:
        print(
            f"Warning: {len(summary.unmatched_motifs)} hit motifs not found in tf2pwms "
            f"(first few: {summary.unmatched_motifs[:5]})",
            file=sys.stderr,
        )


def _cmd_find_tf(args):
    from boltzscan.utils.find_tf import (
        DEFAULT_PFAM,
        find_transcription_factors,
        print_report,
    )

    summary = find_transcription_factors(
        proteins=args.proteins,
        out_dir=args.output,
        pfam=args.pfam or DEFAULT_PFAM,
        cpu=args.cpu,
        force=args.force,
        domtbl=args.domtbl,
    )
    print_report(summary)
    print(
        f"Wrote {summary.annotation_tsv} "
        f"({summary.n_tf}/{summary.n_total} proteins are TFs across "
        f"{len(summary.family_counts)} families; "
        f"{summary.n_te_excluded} transposon-derived dropped)"
    )


def _cmd_build_pwm_refs(args):
    from boltzscan.pwmmap.archive import pack_reference_store
    from boltzscan.pwmmap.cluster import cluster_reference_motifs
    from boltzscan.pwmmap.refs import build_reference_db
    from boltzscan.pwmmap.pfam import runtime_pfam_paths, validate_runtime_pfam

    store = build_reference_db(
        refs_dir=args.refs,
        pfam=args.pfam,
        cpu=args.cpu,
        refresh=args.refresh,
        include_cisbp=not args.no_cisbp,
        include_jaspar=not args.no_jaspar,
    )
    print(f"Reference store ready at {store.root}")
    compact_pfam = validate_runtime_pfam(*runtime_pfam_paths(store.root))
    print(
        f"Runtime Pfam: {compact_pfam.n_profiles} profiles, "
        f"SHA256 {compact_pfam.sha256}"
    )
    clusters_path = Path(args.refs) / 'motif_clusters.tsv'
    if args.refresh or not clusters_path.is_file():
        clusters = cluster_reference_motifs(
            refs_dir=args.refs,
            tomtom=args.tomtom,
            qthresh=args.qthresh,
            cpu=args.cpu,
        )
        print(
            f"Representative PWMs: {clusters.n_clustered} motifs -> "
            f"{clusters.n_clusters} clusters"
        )
    else:
        print(f"Representative PWM map reused: {clusters_path}")

    if not args.no_archive:
        archive = Path(args.archive) if args.archive else (
            Path(args.refs).parent / 'boltzscan_pwm_refs.tar.gz'
        )
        release = pack_reference_store(args.refs, archive)
        print(
            f"Release archive: {release.archive} "
            f"({release.size_bytes / 1024 / 1024:.1f} MiB)"
        )
        print(f"SHA256: {release.sha256}")
        print(f"Checksum file: {release.checksum_file}")


def _cmd_install_pwm_refs(args):
    from boltzscan.pwmmap.archive import install_reference_store

    print(f"PWM reference release: {args.url}")
    print(f"Expected SHA256: {args.sha256}")
    try:
        summary = install_reference_store(
            args.url,
            args.refs,
            args.sha256,
            replace=args.replace,
            progress=lambda message: print(message, flush=True),
        )
    except (FileExistsError, FileNotFoundError, OSError, RuntimeError, ValueError) as exc:
        raise SystemExit(f'boltzscan install-pwm-refs: {exc}') from None
    print(
        f"Installed PWM refs: {summary.n_dbd_rows} DBD rows, "
        f"{summary.n_meme_motifs} scan-ready motifs, "
        f"{summary.n_pfam_profiles} Pfam profiles -> {summary.refs_dir}"
    )
    print(f"Verified SHA256: {summary.archive_sha256}")
    if summary.backup_dir is not None:
        print(f"Previous store moved to: {summary.backup_dir}")


def _cmd_map_pwm(args):
    from boltzscan.pwmmap.mapper import map_species

    try:
        refs_for_mapping = Path(args.refs)
        excluded_species = tuple(
            species.strip() for species in args.exclude_species if species.strip()
        )
        if excluded_species:
            from boltzscan.pwmmap.leaveout import create_reference_subset

            refs_for_mapping = Path(args.output) / 'refs_loso'
            loso = create_reference_subset(
                args.refs,
                refs_for_mapping,
                exclude_species=excluded_species,
            )
            print(
                f"LOSO reference subset: retained {loso.n_retained_dbd_rows}/"
                f"{loso.n_input_dbd_rows} DBD records after excluding "
                f"{', '.join(excluded_species)}"
            )
            if not args.no_pwm_cluster:
                from boltzscan.pwmmap.cluster import cluster_reference_motifs

                clusters = cluster_reference_motifs(
                    refs_dir=refs_for_mapping,
                    cpu=args.cpu,
                )
                print(
                    f"LOSO representative PWMs: {clusters.n_clustered} motifs -> "
                    f"{clusters.n_clusters} clusters"
                )

        s = map_species(species_fasta=args.proteins, out_dir=args.output,
                        refs_dir=refs_for_mapping,
                        domtbl=args.domtbl, threshold_mode=args.threshold_mode,
                        threshold=args.threshold, min_cov=args.min_cov,
                        blastp=args.blastp, makeblastdb=args.makeblastdb,
                        pfam=args.pfam, cpu=args.cpu,
                        collapse_clusters=not args.no_pwm_cluster)
    except (FileNotFoundError, ModuleNotFoundError, RuntimeError, ValueError) as exc:
        raise SystemExit(f'boltzscan map-pwm: {exc}') from None
    print(f"Wrote {s.out_dir}/tf2pwms.json "
          f"({s.n_mapped}/{s.n_species_tfs} TFs mapped, {s.n_motifs} motifs)")


def _cmd_promoter(args):
    """Unified promoter extraction for BED, FASTA, or both."""
    print(f"Extracting promoter {args.format.upper()} using GFF: {args.gff}")
    extract_promoter_regions(args)


def _cmd_predict(args):
    from boltzscan.predict.runners import native_prediction_dir, run_prediction

    try:
        rc = run_prediction(
            input_dir=args.input_dir,
            out_dir=args.output,
            model=args.model,
            seed=args.seed,
        )
        native_output = native_prediction_dir(args.input_dir, args.output, args.model)
    except ValueError as exc:
        raise SystemExit(f'boltzscan predict: {exc}') from None
    if rc != 0:
        raise SystemExit(f"{args.model} prediction exited with code {rc}")
    print(f"{args.model} native prediction done -> {native_output}")
    print("Run `boltzscan wash --help` to publish the shared prediction layout.")


def _cmd_ipsae(args):
    from boltzscan.tfdna import score_tf_dna_run

    processes = available_cpu_count()
    try:
        score_tf_dna_run(
            out_dir=args.run,
            model=args.model,
            processes=processes,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f'boltzscan ipsae: {exc}') from None


def _cmd_txt2meme(args):
    from boltzscan.utils.io_utils import calculate_fasta_background, txt_to_meme

    if args.background is not None and args.bg_fasta is not None:
        raise SystemExit("Use either --background or --bg-fasta, not both")
    if args.bg_fasta is not None:
        background = calculate_fasta_background(args.bg_fasta)
    else:
        background = _parse_background(args.background)
    output_files = txt_to_meme(
        input_path=args.input,
        output_path=args.output,
        nsites=args.nsites,
        background=background,
        force=args.force,
    )
    bg_text = ','.join(f'{value:.4f}' for value in background)
    print(f"Background A,C,G,T: {bg_text}")
    print(f"Wrote {len(output_files)} MEME file(s) to {args.output}")


def _cmd_valid(args):
    from boltzscan.utils.ipsae_score import (
        print_ipsae_warnings,
        validate_ipsae_iptm,
    )

    summary = validate_ipsae_iptm(
        res_dir=args.res_dir,
        output=args.output,
        tolerance=args.tolerance,
        force=args.force,
        processes=args.processes,
    )
    print_ipsae_warnings(summary.warnings)
    if summary.output is not None:
        print(f"Wrote {summary.output}")
    incomplete = summary.total_predictions - summary.compared_predictions
    print(
        f"Validated {summary.compared_predictions}/{summary.total_predictions} predictions "
        f"({summary.valid_predictions} valid, {summary.invalid_predictions} invalid, "
        f"{incomplete} incomplete, "
        f"tolerance={args.tolerance})"
    )

    if summary.total_predictions == 0:
        print("No prediction directories found.", file=sys.stderr)
        raise SystemExit(1)

    failed = summary.result[summary.result['valid'] != True]
    if not failed.empty:
        cols = ['name', 'boltz_iptm', 'ipsae_iptm', 'iptm_diff', 'tolerance', 'valid']
        print("Failed or incomplete validation rows:", file=sys.stderr)
        print(failed.loc[:, cols].head(20).to_string(index=False), file=sys.stderr)
        if len(failed) > 20:
            print(f"... {len(failed) - 20} more", file=sys.stderr)
        raise SystemExit(1)


_DISPATCH = {
    'doctor': _cmd_doctor,
    'run': _cmd_run,
    'wash': _cmd_wash,
    'msa': _cmd_msa,
    'promoter': _cmd_promoter,
    'fimo-scan': _cmd_fimo_scan,
    'fimo2yaml': _cmd_fimo2yaml,
    'find-interface': _cmd_find_interface,
    'hit2fasta': _cmd_candidate_tf_fasta,
    'find-tf': _cmd_find_tf,
    'build-pwm-refs': _cmd_build_pwm_refs,
    'install-pwm-refs': _cmd_install_pwm_refs,
    'map-pwm': _cmd_map_pwm,
    'predict': _cmd_predict,
    'ipsae': _cmd_ipsae,
    'txt2meme': _cmd_txt2meme,
    'valid': _cmd_valid,
}


def _parse_background(value):
    if value is None:
        return [0.25, 0.25, 0.25, 0.25]
    try:
        background = [float(v.strip()) for v in str(value).split(',')]
    except ValueError as exc:
        raise SystemExit(f"Invalid --background value: {value}") from exc
    if len(background) != 4:
        raise SystemExit("--background must contain four comma-separated values for A,C,G,T")
    total = sum(background)
    if total <= 0:
        raise SystemExit("--background frequencies must sum to a positive value")
    if not 0.99 <= total <= 1.01:
        raise SystemExit(f"--background frequencies must sum to 1.0, got {total:.6f}")
    return background


def _command_log_path(args):
    """Choose the run-root log for a command.

    ``boltzscan run`` is the production interface. Lower-level commands infer
    the same root when their output follows the standard RUN layout.
    """
    command = args.command
    if command == 'doctor':
        from boltzscan.toolchain import managed_tool_root

        return managed_tool_root(args.tool_dir) / 'doctor.log'
    if command == 'run':
        return Path(args.run) / 'run.log'
    if command == 'build-pwm-refs':
        return Path(args.refs) / f'{command}.log'
    if command == 'install-pwm-refs':
        return Path(args.refs).parent / 'install-pwm-refs.log'
    if command == 'wash':
        return Path(args.run) / f'{command}.log'
    if command == 'find-interface':
        return Path(args.run) / 'find-interface.log'
    if command == 'msa':
        output = Path(args.output_dir)
        root = output.parent if output.name == 'msa' else output
        return root / 'msa.log'
    if command == 'promoter':
        output = Path(args.output)
        return output.parent / 'promoter.log'
    if command == 'ipsae':
        return Path(args.run) / 'ipsae.log'
    if command == 'valid':
        if getattr(args, 'output', None) is None:
            return Path(args.res_dir) / f'{command}.log'
        output = Path(args.output)
        root = output.parent.parent if output.parent.name == 'results' else output.parent
        return root / 'valid.log'
    output = getattr(args, 'output', None)
    if output is not None:
        output = Path(output)
        if command == 'hit2fasta' or output.suffix == '.meme':
            return output.parent / f'{output.name}.log'
        standard_component = {
            'find-tf': {'find_tf', 'tf'},
            'map-pwm': {'pwm'},
            'fimo-scan': {'scan'},
            'fimo2yaml': {'inputs'},
        }.get(command, set())
        root = output.parent if output.name in standard_component else output
        return root / f'{command}.log'
    return Path.cwd() / f'boltzscan-{command}.log'


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)
    handler = _DISPATCH.get(args.command)
    if handler is None:
        parser.error(f"Unknown command: {args.command}")
    from boltzscan.runlog import CommandLog

    command_argv = list(sys.argv[1:] if argv is None else argv)
    command_log = _command_log_path(args)
    with CommandLog(
        command_log,
        ['boltzscan', *command_argv],
        args,
    ):
        handler(args)


if __name__ == '__main__':
    main()
