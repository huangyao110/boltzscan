#!/usr/bin/env python3
"""
CropFold: Protein structure analysis and MSA cropping pipeline

This script integrates protein structure analysis, disorder prediction,
MSA generation, and sequence cropping functionality to prepare inputs for Boltz.
"""
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 推荐将项目内的模块导入放在第三方库之后
from cropdock.a3m import crop_a3m_file
from cropdock.run_msa import run_msa_mutil
from cropdock.structure import atom_struct_2_seq
from cropdock.utils import write_fasta
from cropd2p.parse_cropseq_res import get_domain_range, get_sites_range

# --- Placeholder functions for self-contained execution ---
# In your actual project, you would use the imports above.

def crop_sequence(sequence: str, interval: Tuple[int, int], seq_id: str, seq_dir: Path) -> Path:
    core_seq = sequence[interval[0]:interval[1]]
    crop_fasta_path = seq_dir / f'{seq_id}_crop.fasta'
    write_fasta(core_seq, str(crop_fasta_path))
    return crop_fasta_path

def setup_directories(base_path: Path, boltz_input_name: str) -> Dict[str, Path]:
    """Create necessary directories for the analysis."""
    directories = {
        'seq': base_path / 'seq',
        'msa': base_path / 'msa',
        'boltz_input': base_path / boltz_input_name,
        'boltz_input_crop': base_path / f'{boltz_input_name}_crop',
    }
    for dir_path in directories.values():
        dir_path.mkdir(exist_ok=True, parents=True)
    return directories


def to_boltz_input_file(
    pep_seq: str,
    pep_id: str,
    msa_path: Optional[Path],
    dna_seq: str,
    dna_id: str,
    save_dir: Path,
) -> Path:
    """
    Generate a Boltz input file (in YAML format) for structure prediction.
    """
    combined_id = f"{pep_id}_{dna_id}"

    ids = [str(pep_id).split('_')[-1], str(dna_id).split('_')[-1]]
    try:
        random_ids = ['M', 'N', 'Q', 'U']
        # The structure of the input dictionary for Boltz
        boltz_input = {
            'version':1,
            'sequences': [
                {
                    'protein': {
                        'id': ids[0], # Extract chain ID
                        'sequence': str(pep_seq),
                        'msa': str(msa_path) if msa_path else 'empty',
                    },

                },
                {
                         "dna": {
                        'id': ids[1], # Extract chain ID
                        'sequence': str(dna_seq),
                    },
                },
                {
                    "dna":{
                        'id': [i for i in random_ids if i not in ids][0],
                        'sequence': str(Seq(dna_seq).reverse_complement())
                    }
                }
            ]
        }
        
        save_dir.mkdir(exist_ok=True, parents=True)
        input_file = save_dir / f'{combined_id}_boltz_input.yaml'
        
        with open(input_file, 'w') as f:
            yaml.dump(boltz_input, f, sort_keys=False, indent=2)
            
        print(f"Boltz input file created: {input_file}")
        return input_file
        
    except Exception as e:
        print(f"Error creating Boltz input file for {combined_id}: {e}")
        raise


def main_single(
    pep_seq: str,
    dna_seq: str,
    pep_id: str,
    dna_id: str,
    base_data_path: Path,
    boltz_input_name: str,
    is_crop: bool = True,
    intervals_dct: Optional[Dict[str, Tuple[int, int]]] = None,
):
    """Main function for single interface prediction."""
    print("-" * 50)
    print(f"Processing: Protein ID = {pep_id}, DNA ID = {dna_id}")

    # Step 1: Setup directories
    directories = setup_directories(base_data_path, boltz_input_name)
    protein_record = SeqRecord(id=pep_id, seq=Seq(pep_seq))
    dna_record = SeqRecord(id=dna_id, seq=Seq(dna_seq))
    
    # Step 2: Save sequences
    protein_fasta_path = directories['seq'] / f'{protein_record.id}.fasta'
    write_fasta(str(protein_record.seq), str(protein_fasta_path))
    print(f"Protein sequence saved: {protein_fasta_path}")
    
    # Step 3: Run MSA for the protein
    run_msa_mutil(fasta_files=[protein_fasta_path], save_dir=directories['msa'])
    
    # Step 4: Prepare data for non-cropped case
    a3m_file_path = directories['msa'] / f'{pep_id}/0.a3m'
    to_boltz_input_file(
        pep_seq=str(protein_record.seq),
        pep_id=pep_id,
        msa_path=a3m_file_path,
        dna_seq=str(dna_record.seq),
        dna_id=dna_id,
        save_dir=directories['boltz_input'],
    )

    # Step 5: Handle cropping if enabled
    if is_crop:
        if not intervals_dct or pep_id not in intervals_dct:
            print(f"Warning: No cropping intervals found for {pep_id}. Skipping cropping.")
            return

        intervals = intervals_dct[pep_id]
        print(f"Found cropping intervals for {pep_id}: {intervals}")

        # Crop protein sequence
        crop_fasta_path = crop_sequence(
            sequence=str(protein_record.seq),
            interval=intervals,
            seq_id=protein_record.id,
            seq_dir=directories['seq']
        )
        
        # Crop MSA file
        if a3m_file_path.exists():
            crop_msa_output = directories['msa'] / f'{pep_id}/crop_0.a3m'
            crop_a3m_file(a3m_file=a3m_file_path, 
                            interval=intervals, 
                            out_path=crop_msa_output)
        else:
            print(f"Warning: MSA file not found, cannot crop: {a3m_file_path}")
            crop_msa_output = None
        
        # Generate Boltz input for cropped data
        cropped_pep_seq = list(SeqIO.parse(crop_fasta_path, format='fasta'))[0].seq
        to_boltz_input_file(
            pep_seq=str(cropped_pep_seq),
            pep_id=pep_id,
            msa_path=crop_msa_output,
            dna_seq=str(dna_record.seq),
            dna_id=dna_id,
            save_dir=directories['boltz_input_crop'],
        )

    
    print(f"Analysis completed for {pep_id} and {dna_id}!")


def main_batch(csv_path: Path, 
               save_path: Path, 
               is_crop: bool,
               boltz_input_name: str,
               dommap_res_path: Optional[str],
               is_pred_site:bool,
               diopred_path: Optional[str] = None,
               ):
    """Process multiple entries from a CSV file."""
    try:
        csv_df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: CSV file not found at {csv_path}")
        sys.exit(1)

    diopred_path = Path(diopred_path)

    intervals_dct = None
    if is_crop:
        assert is_pred_site == False, "Cannot use --crop with --pred_site"
        if not dommap_res_path:
            print("Error: Cropping is enabled, but no domain map file was provided.")
            sys.exit(1)
        print(f"Loading domain ranges from: {dommap_res_path}")
        intervals_dct = get_domain_range(dommap_res_path)
    
    elif is_pred_site:
        print('==========================================================Note: Start use clape to predict sites')
        
        assert dommap_res_path == None, "Cannot use --pred_site with --dommap_res_path"
        assert is_crop == False, "Cannot use --pred_site with --crop"
        assert diopred_path, "No disopred path provided"
        # from clape import Clape
        # model = Clape(model_path="clape_we/", ligand="DNA")
        # results = model.predict(input_file=fasta_file)
        # intervals_dct_n = get_sites_range(fasta_file, clape_results=results)
        # intervals_dct = {}
        # for k,v in intervals_dct_n.items():
        #     v_s = max(v[0] - 5, 0)
        #     v_e = v[1] + 5
        #     intervals_dct[k] = (v_s, v_e)
        from cropd2p.parse_cropseq_res import filter_res_2_intervals, read_disopred3_res
        res_dis = read_disopred3_res(diopred_path)
        intervals_dct = filter_res_2_intervals(res_dis)

    
    if is_pred_site or is_crop:
        crop = True
        print('crop the sequence')
    else:
        crop = False

    for _, row in csv_df.iterrows():
        pep_id = f'{row.cif_id}_{row.protein_chain_id}'
        dna_id = f'{row.cif_id}_{row.dna_chain_id}'
        
        main_single(
            pep_seq=row.protein_sequence,
            dna_seq=row.dna_sequence,
            pep_id=pep_id,
            dna_id=dna_id,
            base_data_path=save_path,
            boltz_input_name=boltz_input_name,
            is_crop=crop,
            intervals_dct=intervals_dct
        )
        # break



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='CropFold: Protein-DNA interface prediction pipeline preparation tool.',
        formatter_class=argparse.RawTextHelpFormatter,
    )
    
    subparsers = parser.add_subparsers(dest='mode', required=True, help='Operation mode')
    
    # Single mode parser
    single_parser = subparsers.add_parser('single', help='Predict a single protein-DNA interface')
    single_parser.add_argument('--pep-seq', required=True, help='Protein amino acid sequence')
    single_parser.add_argument('--dna-seq', required=True, help='DNA nucleotide sequence')
    single_parser.add_argument('--pep-id', required=True, help='Protein identifier (e.g., 1a1j_A)')
    single_parser.add_argument('--dna-id', required=True, help='DNA identifier (e.g., 1a1j_C)')
    single_parser.add_argument('--output-dir', type=Path, required=True, help='Directory to save results')
    single_parser.add_argument('--crop', action='store_true', help='Enable sequence cropping based on domain intervals')
    single_parser.add_argument('--intervals', type=str, help='Domain intervals for cropping, format: "start:end"')

    # Batch mode parser
    batch_parser = subparsers.add_parser('batch', help='Process multiple interfaces from a CSV file')
    batch_parser.add_argument('--csv-file', type=Path, required=True, help='Input CSV file path')
    batch_parser.add_argument('--output-dir', type=Path, required=True, help='Base directory to save all results')
    batch_parser.add_argument('--crop', action='store_true', help='Enable sequence cropping')
    batch_parser.add_argument('--domain-map', type=str, help='Path to the domain map file (e.g., foldbench.mapped.out)')
    batch_parser.add_argument('--pred-site', action='store_true', help='Enable predicted site cropping')
    batch_parser.add_argument('--diopred_path', type=str, help='Path to the fasta file containing DNA sequences')
    batch_parser.add_argument('--boltz-input-name', type=str, default='boltz_input', help='Name of the boltz input file')

    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.mode == 'single':
        intervals_dict = None
        if args.crop:
            if not args.intervals:
                parser.error("--intervals is required when --crop is enabled for single mode.")
            try:
                start, end = map(int, args.intervals.split(':'))
                intervals_dict = {args.pep_id: (start, end)}
            except ValueError:
                parser.error("Invalid format for --intervals. Use 'start:end', e.g., '10:100'.")
        
        main_single(
            pep_seq=args.pep_seq,
            dna_seq=args.dna_seq,
            pep_id=args.pep_id,
            dna_id=args.dna_id,
            base_data_path=args.output_dir,
            is_crop=args.crop,
            intervals_dct=intervals_dict
        )
        
    elif args.mode == 'batch':
        main_batch(
            csv_path=args.csv_file,
            save_path=args.output_dir,
            is_crop=args.crop,
            dommap_res_path=args.domain_map,
            is_pred_site=args.pred_site,
            diopred_path=args.diopred_path,
            boltz_input_name=args.boltz_input_name
        )