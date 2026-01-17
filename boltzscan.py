from cropd2p.cropdock import interaction
from cropd2p.cropdock.disopred import batch_run_disoreder
import argparse
from pathlib import Path
import subprocess
import sys
from fimocistarget.expro import extract_promoters, extract_bed_regions, extract_promoters_both



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='BoltzScan: A tool for construct GRN')
    sub_parsers = parser.add_subparsers(dest='command')
    sub_parsers.required = True

    # subparser for extract_pd_from_pdb
    extract_pd_from_pdb_parser = sub_parsers.add_parser('extract_pd_from_pdb', help='Extract protein domains from PDB files')
    extract_pd_from_pdb_parser.add_argument('pdb_file', help='PDB file')

    # subparser for disopred
    disopred_parser = sub_parsers.add_parser('disopred', help='Run DISOPRED on multiple FASTA files')
    disopred_parser.add_argument('-t', '--algo_type', required=True, help='disopred algorithm type: disopred3 or iupred2')
    disopred_parser.add_argument('-i', '--fasta_file', required=True, help='Input FASTA file containing sequences')
    disopred_parser.add_argument('-s', '--script_path', help='Path to the run_disopred.pl script')
    disopred_parser.add_argument('-c', '--cpus', type=int, default=4, help='Number of CPU cores to use (default: 4)')
    disopred_parser.add_argument('-o', '--output_dir', required=True, default='disopred_res', help='Output directory for results (default: disopred_results)')

    # subparser for i (run_msa_mutil)
    i_parser = sub_parsers.add_parser('msa', help='Run MSA on FASTA files using run_msa_mutil')
    i_parser.add_argument('-f', '--fasta_file', required=True, help='Input FASTA file path')
    i_parser.add_argument('-o', '--output_dir', required=True, help='Output directory for MSA results')
    i_parser.add_argument('-c', '--cpus', type=int, default=1, help='Number of CPU cores to use (default: 1)')

    # subparser for promoter extraction
    promoter_parser = sub_parsers.add_parser('promoter', help='Extract promoter regions from genome using GFF annotations')
    promoter_parser.add_argument('-gff', '--gff', required=False,
                        default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.60.gff3',
                        help='Input GFF3 annotation file (can be .gz). Assumes ### gene separator or groups implicitly.')
    promoter_parser.add_argument('-g', '--genome', required=False,
                        default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.dna_sm.toplevel.fa',
                        help='Input Genome FASTA file (can be .gz).')
    promoter_parser.add_argument('-o', '--output', required=False,
                        default = './Rosa_chinensis_promoter.fasta',
                        help='Output FASTA file name (all promoters written here).')
    promoter_parser.add_argument('-u', '--upstream', type=int, default=2000, help='Bases UPSTREAM of CDS anchor (will be uppercase).')
    promoter_parser.add_argument('-d', '--downstream', type=int, default=200, help='Bases DOWNSTREAM of CDS anchor, including anchor position (will be lowercase).')
    
    # subparser for BED extraction
    bed_parser = sub_parsers.add_parser('extract-bed', help='Extract genomic regions from GFF and generate BED file')
    bed_parser.add_argument('-i', '--input', required=False,
                      default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.60.gff3',
                      help='Input GFF file path')
    bed_parser.add_argument('-o', '--output', required=False,
                      default='./Rosa_chinensis_regions.bed',
                      help='Output BED file path')
    bed_parser.add_argument('-t', '--type', default='gene',
                      help='Feature type to extract (default: gene)')
    bed_parser.add_argument('-n', '--name', default='ID',
                      help='Attribute field to use as name (default: ID)')
    
    # subparser for promoter extraction to both BED and FASTA
    promoter_both_parser = sub_parsers.add_parser('promoter-both', help='Extract promoter regions to both BED and FASTA files')
    promoter_both_parser.add_argument('-gff', '--gff', required=False,
                              default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.60.gff3',
                              help='Input GFF3 annotation file (can be .gz). Assumes ### gene separator or groups implicitly.')
    promoter_both_parser.add_argument('-g', '--genome', required=False,
                              default='./odata/Rosa_chinensis/Rosa_chinensis.RchiOBHm-V2.dna_sm.toplevel.fa',
                              help='Input Genome FASTA file (can be .gz).')
    promoter_both_parser.add_argument('-o', '--output', required=False,
                              default='./Rosa_chinensis_promoter',
                              help='Output file prefix, will generate .bed and .fasta files')
    promoter_both_parser.add_argument('-u', '--upstream', type=int, default=2000, help='Bases UPSTREAM of CDS anchor (will be uppercase).')
    promoter_both_parser.add_argument('-d', '--downstream', type=int, default=200, help='Bases DOWNSTREAM of CDS anchor, including anchor position (will be lowercase).')

    # Parse arguments
    args = parser.parse_args()

    # Handle subcommands
    if args.command == 'extract_pd_from_pdb':
        # Handle extract_pd_from_pdb subcommand
        print(f"Extracting protein domains from PDB file: {args.pdb_file}")
        # Add your extract_pd_from_pdb logic here
        pass
    elif args.command == 'disopred':
        if args.algo_type  == 'disopred3':
        # Handle disopred subcommand
            print(f"Running DISOPRED on FASTA file: {args.fasta_file}")
            import sys


            # Create a mock sys.argv to pass to batch_run_disoreder.main()
            original_argv = sys.argv
            sys.argv = [
                'batch_run_disoreder',
                '-i', args.fasta_file,
                '-s', args.script_path,
                '-c', str(args.cpus),
                '-o', args.output_dir
            ]

            try:
                batch_run_disoreder.main()
            finally:
                sys.argv = original_argv
        elif args.algo_type  == 'iupred2':
            print(f"Running IUPRED2 on FASTA file: {args.fasta_file}")
            # Add your iupred2 logic here
            iupred_script = 'cropd2p/iupred3/iupred3.py'
            batch_run_disoreder.split_fasta(args.fasta_file, args.output_dir)
            files_lst = [str(i) for i in Path(args.output_dir).iterdir() if i.suffix == '.fasta']
            for i in files_lst:
                cmd = f"python {iupred_script} {i} short > {args.output_dir}/{Path(i).stem}_iupred_res.txt"
                print(f"Executing command: {cmd}")
                subprocess.run(cmd, shell=True, check=True)

    elif args.command == 'msa':
        # Handle i subcommand (run_msa_mutil)
        print(f"Running MSA on FASTA file: {args.fasta_file}")
        from cropd2p.cropdock.run_msa import run_msa_mutil
        batch_run_disoreder.split_fasta(args.fasta_file, args.output_dir)
        files_lst = [str(i) for i in Path(args.output_dir).iterdir() if i.suffix == '.fasta']

        # Call run_msa_mutil with the provided arguments
        result = run_msa_mutil(
            fasta_files=files_lst,
            save_dir=args.output_dir,
            n_cpu=args.cpus
        )

        if result:
            print(f"MSA completed successfully. Results saved to: {args.output_dir}")
        else:
            print("MSA processing failed or no results generated.")

    elif args.command == 'promoter':
        # Handle promoter extraction subcommand
        print(f"Extracting promoter regions using GFF: {args.gff}")
        extract_promoters(args)
    
    elif args.command == 'extract-bed':
        # Handle BED extraction subcommand
        print(f"Extracting genomic regions from GFF: {args.input}")
        extract_bed_regions(args.input, args.output, args.type, args.name)
    
    elif args.command == 'promoter-both':
        # Handle promoter extraction to both BED and FASTA subcommand
        print(f"Extracting promoter regions to both BED and FASTA using GFF: {args.gff}")
        extract_promoters_both(args)
