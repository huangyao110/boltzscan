"""CLI entry point for the boltzscan package.

Exposed as the `boltzscan` console script via pyproject.toml. The thin
`bscan.py` shim at the repo root and `boltzscan/__main__.py` both delegate to
:func:`main` here.
"""
import argparse
import subprocess
import sys
from pathlib import Path

from boltzscan import __version__
from boltzscan.fimocistarget.expro import (
    extract_promoters,
    extract_bed_regions,
    extract_promoters_both,
)
from boltzscan.fimocistarget.runner import run_cistarg_with_tf_aggregation
from boltzscan.utils.boltz_input import build_boltz_input_from_fimo

__author__ = 'huangyao'
__date__ = '2026-05-16'

# Vendored Boltz CLI lives inside the package tree but is not a Python
# submodule (it has its own pyproject); we shell out to its main.py.
_PKG_ROOT = Path(__file__).resolve().parent
_BOLTZ_MAIN = _PKG_ROOT / 'boltz' / 'src' / 'boltz' / 'main.py'


def _build_parser():
    parser = argparse.ArgumentParser(
        prog='boltzscan',
        description='BoltzScan: GRN inference via structural modeling of TFs and their targets.',
        epilog=(
            f'Author : {__author__}\n'
            f'Date   : {__date__}\n'
            f'Version: {__version__}'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-V', '--version', action='version',
        version=f'boltzscan {__version__} (by {__author__}, {__date__})')
    sub_parsers = parser.add_subparsers(dest='command')
    sub_parsers.required = True

    # extract_pd_from_pdb
    p = sub_parsers.add_parser('extract_pd_from_pdb', help='Extract protein domains from PDB files')
    p.add_argument('pdb_file', help='PDB file')

    # msa
    p = sub_parsers.add_parser('msa', help='Run MSA on FASTA files using run_msa_mutil')
    p.add_argument('-f', '--fasta_file', required=True, help='Input FASTA file path')
    p.add_argument('-o', '--output_dir', required=True, help='Output directory for MSA results')
    p.add_argument('-c', '--cpus', type=int, default=1, help='Number of CPU cores to use (default: 1)')

    # promoter
    p = sub_parsers.add_parser('promoter', help='Extract promoter regions from genome using GFF annotations')
    p.add_argument('-gff', '--gff', required=False,
        default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.60.gff3',
        help='Input GFF3 annotation file (can be .gz). Assumes ### gene separator or groups implicitly.')
    p.add_argument('-g', '--genome', required=False,
        default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.dna_sm.toplevel.fa',
        help='Input Genome FASTA file (can be .gz).')
    p.add_argument('-o', '--output', required=False,
        default='./Rosa_chinensis_promoter.fasta',
        help='Output FASTA file name (all promoters written here).')
    p.add_argument('-u', '--upstream', type=int, default=2000, help='Bases UPSTREAM of CDS anchor (will be uppercase).')
    p.add_argument('-d', '--downstream', type=int, default=200, help='Bases DOWNSTREAM of CDS anchor, including anchor position (will be lowercase).')

    # extract-bed
    p = sub_parsers.add_parser('extract-bed', help='Extract genomic regions from GFF and generate BED file')
    p.add_argument('-i', '--input', required=False,
        default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.60.gff3',
        help='Input GFF file path')
    p.add_argument('-o', '--output', required=False,
        default='./Rosa_chinensis_regions.bed',
        help='Output BED file path')
    p.add_argument('-t', '--type', default='gene', help='Feature type to extract (default: gene)')
    p.add_argument('-n', '--name', default='ID', help='Attribute field to use as name (default: ID)')

    # promoter-both
    p = sub_parsers.add_parser('promoter-both', help='Extract promoter regions to both BED and FASTA files')
    p.add_argument('-gff', '--gff', required=False,
        default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.60.gff3',
        help='Input GFF3 annotation file (can be .gz). Assumes ### gene separator or groups implicitly.')
    p.add_argument('-g', '--genome', required=False,
        default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.dna_sm.toplevel.fa',
        help='Input Genome FASTA file (can be .gz).')
    p.add_argument('-o', '--output', required=False,
        default='./Rosa_chinensis_promoter',
        help='Output file prefix, will generate .bed and .fasta files')
    p.add_argument('-u', '--upstream', type=int, default=2000, help='Bases UPSTREAM of CDS anchor (will be uppercase).')
    p.add_argument('-d', '--downstream', type=int, default=200, help='Bases DOWNSTREAM of CDS anchor, including anchor position (will be lowercase).')

    # boltzscan
    p = sub_parsers.add_parser('boltzscan', help='Run BoltzScan on FASTA files')
    p.add_argument('-fd', '--fasta_dir', required=True, help='Input FASTA file path')
    p.add_argument('-od', '--output_dir', required=True, help='Output directory for BoltzScan results')
    p.add_argument('--sampling_steps', type=int, default=2, help='Sampling steps for BoltzScan')
    p.add_argument('-m', '--model_name', type=str, default='boltz_ode', required=True, help='Model name for BoltzScan')
    p.add_argument('-or', '--override', type=str, default='',
        help='Pass --override to boltz (any non-empty value enables it; rerun existing predictions)')
    p.add_argument('-pt', '--preprocessing-threads', type=int, default=1,
        help='CPU threads for the boltz preprocessing step (default: 1)')
    p.add_argument('-nw', '--num-workers', type=int, default=2,
        help='DataLoader workers during predict. High values (>8) often deadlock with CUDA. (default: 2)')

    # cistarg
    p = sub_parsers.add_parser(
        'cistarg',
        help='Build cisTarget DB (FIMO + score matrix) and roll motif scores up to TF level',
    )
    p.add_argument('-f', '--fasta', required=True,
        help='Target promoter FASTA (also used to compute ACGT background frequencies)')
    p.add_argument('-m', '--motif-dir', required=True,
        help='Directory of per-motif .meme files (e.g. data/pwms/tair/pwms_all_memes)')
    p.add_argument('-t', '--tf2pwms', required=True,
        help='JSON mapping TF_id -> [motif_id, ...] (e.g. data/pwms/tair/tf2pwms_dct.json)')
    p.add_argument('-o', '--output', required=True,
        help='Output directory for cisTarget DB and TF-level matrices')
    p.add_argument('--tf-fasta', default=None,
        help='Optional TF FASTA; if given, only motifs mapped to these TF IDs (via tf2pwms) are scanned')
    p.add_argument('--motif-agg', default='max',
        choices=['max', 'sum', 'mean', 'count', 'top_n_mean', 'weighted_sum', 'max_with_count'],
        help='Per-motif score aggregation across multiple FIMO hits in a sequence (default: max)')
    p.add_argument('--tf-agg', default='max', choices=['max', 'sum', 'mean'],
        help='How to combine multiple motifs that map to the same TF (default: max)')
    p.add_argument('--pvalue-thresh', type=float, default=1e-4, help='FIMO p-value threshold (default: 1e-4)')
    p.add_argument('--max-stored-scores', type=int, default=500000, help='FIMO max_stored_scores (default: 500000)')
    p.add_argument('--motif-pseudo', type=float, default=0.1, help='Motif pseudocount (default: 0.1)')
    p.add_argument('--top-n', type=int, default=3, help='top_n for top_n_mean aggregation (default: 3, ignored otherwise)')
    p.add_argument('--no-analyze', action='store_true', help='Skip multi-hit statistics step to save time')
    p.add_argument('--reuse', action='store_true',
        help='Reuse <output>/score_matrix.feather if present (skip FIMO rebuild)')

    # fimo2boltz
    p = sub_parsers.add_parser(
        'fimo2boltz',
        help='Filter FIMO hits, attach TF/promoter/MSA metadata, emit per-row Boltz input FASTAs',
    )
    p.add_argument('-f', '--fimo', required=True,
        help='Raw FIMO results CSV (e.g. tasks/llx/cisdb/fimo_results_raw.csv.gz)')
    p.add_argument('-t', '--tf-fasta', required=True,
        help='TF protein FASTA; record IDs must match TF ids in tf2pwms')
    p.add_argument('-j', '--tf2pwms', required=True,
        help='JSON mapping TF_id -> [motif_id, ...]')
    p.add_argument('-r', '--promoter-fasta', required=True,
        help='Promoter FASTA whose record IDs match sequence_name in the FIMO results')
    p.add_argument('-m', '--msa-dir', required=True,
        help='Directory containing per-TF MSA subdirs (each with 0.a3m); subdir names use "_" not "."')
    p.add_argument('-o', '--output', required=True, help='Output directory; FASTAs land in <output>/boltzscan_input/')
    p.add_argument('-p', '--pvalue', type=float, default=1e-5, help='Drop FIMO hits with pvalue > this (default: 1e-5)')
    p.add_argument('-O', '--overlap-thresh', type=float, default=0.0,
        help='Overlap-ratio threshold for filter_fimo_by_seq_and_overlap (default: 0.0)')
    p.add_argument('-e', '--extend-range', type=int, default=5,
        help='Bases added on each side when extracting the motif sub-sequence (default: 5)')
    p.add_argument('-L', '--layout', default='auto', choices=['sub', 'all', 'auto'],
        help='Output layout: "sub"=split BOLTZ_INPUT/ into subdir_NNN/, "all"=flat, '
             '"auto"=sub if > --auto-threshold files else all (default: auto)')
    p.add_argument('-a', '--auto-threshold', type=int, default=20000,
        help='File-count cutoff used by --layout auto (default: 20000)')
    p.add_argument('-n', '--num-subdirs', type=int, default=100,
        help='Number of subdirectories to split BOLTZ_INPUT/ into (default: 100)')
    p.add_argument('-w', '--max-workers', type=int, default=8, help='Threads for parallel FASTA writes (default: 8)')

    # hit2fasta
    p = sub_parsers.add_parser(
        'hit2fasta',
        help='Build a candidate TF protein FASTA from motif hit statistics',
    )
    p.add_argument('-s', '--hit-statistics', required=True,
        help='hit_statistics.csv from cisTarget/FIMO output')
    p.add_argument('-j', '--tf2pwms', required=True,
        help='JSON mapping TF_id -> [motif_id, ...]')
    p.add_argument('-p', '--protein-fasta', required=True,
        help='Protein FASTA containing TF sequences')
    p.add_argument('-o', '--output', required=True,
        help='Output candidate TF protein FASTA')
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
        help='hmmpress-ed Pfam-A.hmm (default: data/pfam/Pfam-A.hmm)')
    p.add_argument('-c', '--cpu', type=int, default=8,
        help='Threads for hmmsearch (default: 8)')
    p.add_argument('--domtbl', default=None,
        help='Reuse an existing hmmsearch domtblout instead of running hmmsearch')
    p.add_argument('--force', action='store_true',
        help='Re-run hmmsearch even if <output>/pfam.domtbl already exists')

    # build-pwm-refs
    p = sub_parsers.add_parser('build-pwm-refs',
        help='Stage A: download cisBP+JASPAR, extract reference DBDs, cache the reference store')
    p.add_argument('--refs', default='data/pwms/_refs', help='Reference store dir (default: data/pwms/_refs)')
    p.add_argument('-m', '--pfam', default=None)
    p.add_argument('-c', '--cpu', type=int, default=8)
    p.add_argument('--no-cisbp', action='store_true')
    p.add_argument('--no-jaspar', action='store_true')
    p.add_argument('--refresh', action='store_true', help='Re-download even if present')

    # map-pwm
    p = sub_parsers.add_parser('map-pwm',
        help='Stage B: map a species TF FASTA to reference motifs via DBD %%ID (no network)')
    p.add_argument('-f', '--proteins', required=True, help='Species TF protein FASTA')
    p.add_argument('-o', '--output', required=True, help='Output dir, e.g. data/pwms/cm_pwms')
    p.add_argument('--refs', default='data/pwms/_refs')
    p.add_argument('--domtbl', default=None, help='Reuse a find-tf pfam.domtbl (else run hmmsearch)')
    p.add_argument('--threshold-mode', default='family', choices=['family', 'global'])
    p.add_argument('--threshold', type=float, default=0.70, help='Global %%ID cutoff (global mode)')
    p.add_argument('--min-cov', type=float, default=0.8, help='Min DBD coverage for a blast hit')
    p.add_argument('--blastp', default=None)
    p.add_argument('--makeblastdb', default=None)
    p.add_argument('-m', '--pfam', default=None)
    p.add_argument('-c', '--cpu', type=int, default=8)
    p.add_argument('--collapse-clusters', action='store_true',
        help='Collapse transferred motifs to cluster representatives (needs cluster-motifs first)')

    # cluster-motifs
    p = sub_parsers.add_parser('cluster-motifs',
        help='Cluster reference motifs per family (Tomtom) into a non-redundant set -> motif_clusters.tsv')
    p.add_argument('--refs', default='data/pwms/_refs', help='Reference store dir')
    p.add_argument('--qthresh', type=float, default=0.05,
        help='Tomtom q-value cutoff for merging motifs (default: 0.05)')
    p.add_argument('--tomtom', default=None, help='Path to tomtom (default: PATH or the meme conda env)')
    p.add_argument('-c', '--cpu', type=int, default=8, help='Families clustered in parallel (default: 8)')

    # cpg-islands
    p = sub_parsers.add_parser('cpg-islands',
        help='Scan a promoter FASTA for CpG islands (cisdb for CpG-binding/CXXC factors, e.g. DDR1)')
    p.add_argument('-f', '--fasta', required=True, help='Promoter FASTA')
    p.add_argument('-o', '--output', required=True, help='Output dir (cpg_per_promoter.tsv, cpg_islands.bed, cpg_target_genes.txt)')
    p.add_argument('--window', type=int, default=200, help='Min island length / window bp (default: 200)')
    p.add_argument('--min-gc', type=float, default=0.5, help='Min GC fraction (default: 0.5)')
    p.add_argument('--min-oe', type=float, default=0.6, help='Min observed/expected CpG ratio (default: 0.6)')

    # ipsae
    p = sub_parsers.add_parser(
        'ipsae',
        help='Calculate TF-DNA ipSAE ranking scores and merge them into a score table',
    )
    p.add_argument('-r', '--res-dir', required=True,
        help='Boltz result directory or its predictions/ subdirectory')
    p.add_argument('-s', '--score-file', default=None,
        help='Optional input score table (.csv, .csv.gz, .tsv, or .tsv.gz)')
    p.add_argument('-o', '--output', default=None,
        help='Output CSV path for the scored table (default: <res-dir>/ipsae_scored.csv)')
    p.add_argument('-i', '--id-col', default=None,
        help='Optional score table column matching Boltz prediction directory names')
    p.add_argument('-f', '--force', action='store_true',
        help='Recalculate IPSAE even when an existing *_10_10.txt file is present')
    p.add_argument('-p', '--processes', type=int, default=1,
        help='Number of concurrent IPSAE calculations (default: 1)')

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

    # esm-embed
    p = sub_parsers.add_parser(
        'esm-embed',
        help='Embed protein sequences with ESMC-6B (runs in the esmfold2 conda env)',
    )
    p.add_argument('-f', '--fasta', required=True, help='Input protein FASTA')
    p.add_argument('-o', '--output', required=True, help='Output .npz path for embeddings')
    p.add_argument('-w', '--weights', default=None,
        help='ESMC-6B weights directory (default: $BOLTZSCAN_ESMC_WEIGHTS or the local esm_weights/ESMC-6B)')
    p.add_argument('-p', '--pooling', default='mean', choices=['mean', 'cls', 'none'],
        help='mean/cls -> [N,D]; none -> per-residue [L,D] arrays (default: mean)')
    p.add_argument('-d', '--device', default='cuda', help='torch device (default: cuda)')
    p.add_argument('-b', '--batch-size', type=int, default=8, help='Sequences per batch (default: 8)')
    p.add_argument('-L', '--max-len', type=int, default=None,
        help='Truncate sequences to this many residues (default: no truncation)')
    p.add_argument('--dtype', default='bfloat16', choices=['bfloat16', 'float16', 'float32'],
        help='Model/compute dtype (default: bfloat16; forced to float32 on CPU)')
    p.add_argument('--python-exe', default=None,
        help='Python interpreter to run the worker (default: $BOLTZSCAN_ESM_PYTHON or the esmfold2 env)')

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
    from boltzscan.cropd2p.run_msa import run_msa_mutil
    from boltzscan.cropd2p.disopred import batch_run_disoreder

    print(f"Running MSA on FASTA file: {args.fasta_file}")
    batch_run_disoreder.split_fasta(args.fasta_file, args.output_dir)
    files_lst = [str(i) for i in Path(args.output_dir).iterdir() if i.suffix == '.fasta']
    result = run_msa_mutil(fasta_files=files_lst, save_dir=args.output_dir, n_cpu=args.cpus)
    if result:
        print(f"MSA completed successfully. Results saved to: {args.output_dir}")
    else:
        print("MSA processing failed or no results generated.")


def _cmd_boltzscan(args):
    print(f"Running BoltzScan on FASTA files: {args.fasta_dir}")
    cmd = [
        sys.executable, str(_BOLTZ_MAIN), 'predict', args.fasta_dir,
        '--out_dir', args.output_dir,
        '--sampling_steps', str(args.sampling_steps),
        '--seed', '42',
        '--write_full_pae',
        '--step_scale', '1.0',
        '--preprocessing-threads', str(args.preprocessing_threads),
        '--num_workers', str(args.num_workers),
        '--model', args.model_name,
    ]
    if args.override:
        cmd.append('--override')
    subprocess.run(cmd)


def _cmd_cistarg(args):
    run_cistarg_with_tf_aggregation(
        target_fasta=args.fasta,
        motif_dir=args.motif_dir,
        tf2pwms_path=args.tf2pwms,
        output_dir=args.output,
        motif_agg=args.motif_agg,
        tf_agg=args.tf_agg,
        pvalue_thresh=args.pvalue_thresh,
        max_stored_scores=args.max_stored_scores,
        motif_pseudo=args.motif_pseudo,
        top_n=args.top_n,
        analyze_hits=not args.no_analyze,
        tf_fasta=args.tf_fasta,
        reuse=args.reuse,
    )


def _cmd_fimo2boltz(args):
    build_boltz_input_from_fimo(
        fimo_csv=args.fimo,
        tf_fasta=args.tf_fasta,
        tf2pwms_path=args.tf2pwms,
        promoter_fasta=args.promoter_fasta,
        msa_dir=args.msa_dir,
        output_dir=args.output,
        pvalue_thresh=args.pvalue,
        overlap_thresh=args.overlap_thresh,
        extend_range=args.extend_range,
        layout=args.layout,
        auto_threshold=args.auto_threshold,
        num_subdirs=args.num_subdirs,
        max_workers=args.max_workers,
    )


def _cmd_candidate_tf_fasta(args):
    from boltzscan.utils.tf_candidates import build_candidate_tf_fasta

    summary = build_candidate_tf_fasta(
        hit_statistics=args.hit_statistics,
        tf2pwms_path=args.tf2pwms,
        protein_fasta=args.protein_fasta,
        output=args.output,
        motif_col=args.motif_col,
    )
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
    from boltzscan.pwmmap.refs import build_reference_db
    store = build_reference_db(refs_dir=args.refs, pfam=args.pfam, cpu=args.cpu,
                               refresh=args.refresh,
                               include_cisbp=not args.no_cisbp,
                               include_jaspar=not args.no_jaspar)
    print(f"Reference store ready at {store.root}")


def _cmd_map_pwm(args):
    from boltzscan.pwmmap.mapper import map_species
    s = map_species(species_fasta=args.proteins, out_dir=args.output, refs_dir=args.refs,
                    domtbl=args.domtbl, threshold_mode=args.threshold_mode,
                    threshold=args.threshold, min_cov=args.min_cov,
                    blastp=args.blastp, makeblastdb=args.makeblastdb,
                    pfam=args.pfam, cpu=args.cpu,
                    collapse_clusters=args.collapse_clusters)
    print(f"Wrote {s.out_dir}/tf2pwms.json "
          f"({s.n_mapped}/{s.n_species_tfs} TFs mapped, {s.n_motifs} motifs)")


def _cmd_cluster_motifs(args):
    from boltzscan.pwmmap.cluster import cluster_reference_motifs
    s = cluster_reference_motifs(refs_dir=args.refs, tomtom=args.tomtom,
                                 qthresh=args.qthresh, cpu=args.cpu)
    print(f"Wrote {s.clusters_tsv} "
          f"({s.n_clustered} motifs with PWM -> {s.n_clusters} clusters "
          f"across {s.n_families} families; {s.n_motifs} total in map)")


def _cmd_cpg_islands(args):
    from boltzscan.fimocistarget.cpg import scan_promoters_for_cpg
    s = scan_promoters_for_cpg(args.fasta, args.output, window=args.window,
                               min_gc=args.min_gc, min_oe=args.min_oe)
    print(f"Wrote {s.out_dir} "
          f"({s.n_with_island}/{s.n_promoters} promoters with a CpG island, "
          f"{s.n_islands} islands)")


def _cmd_ipsae(args):
    from boltzscan.utils.ipsae_score import print_ipsae_warnings, score_ipsae_table

    summary = score_ipsae_table(
        res_dir=args.res_dir,
        score_file=args.score_file,
        output=args.output,
        id_col=args.id_col,
        force=args.force,
        processes=args.processes,
    )
    print_ipsae_warnings(summary.warnings)
    print(
        f"Wrote {summary.output} "
        f"({summary.matched_rows} matched rows, "
        f"{summary.scored_predictions}/{summary.total_predictions} predictions scored)"
    )


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


def _cmd_esm_embed(args):
    from boltzscan.esm_embed.runner import run_esmc_embed

    out = run_esmc_embed(
        fasta=args.fasta,
        output=args.output,
        weights=args.weights,
        pooling=args.pooling,
        device=args.device,
        batch_size=args.batch_size,
        max_len=args.max_len,
        dtype=args.dtype,
        python_exe=args.python_exe,
    )
    print(f"Wrote ESMC-6B embeddings to {out}")


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
    'msa': _cmd_msa,
    'promoter': lambda a: (print(f"Extracting promoter regions using GFF: {a.gff}"), extract_promoters(a)),
    'extract-bed': lambda a: (print(f"Extracting genomic regions from GFF: {a.input}"),
                              extract_bed_regions(a.input, a.output, a.type, a.name)),
    'promoter-both': lambda a: (print(f"Extracting promoter regions to both BED and FASTA using GFF: {a.gff}"),
                                extract_promoters_both(a)),
    'boltzscan': _cmd_boltzscan,
    'cistarg': _cmd_cistarg,
    'fimo2boltz': _cmd_fimo2boltz,
    'hit2fasta': _cmd_candidate_tf_fasta,
    'find-tf': _cmd_find_tf,
    'build-pwm-refs': _cmd_build_pwm_refs,
    'map-pwm': _cmd_map_pwm,
    'cluster-motifs': _cmd_cluster_motifs,
    'cpg-islands': _cmd_cpg_islands,
    'ipsae': _cmd_ipsae,
    'txt2meme': _cmd_txt2meme,
    'esm-embed': _cmd_esm_embed,
    'valid': _cmd_valid,
    'extract_pd_from_pdb': lambda a: print(f"Extracting protein domains from PDB file: {a.pdb_file}"),
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


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)
    handler = _DISPATCH.get(args.command)
    if handler is None:
        parser.error(f"Unknown command: {args.command}")
    handler(args)


if __name__ == '__main__':
    main()
