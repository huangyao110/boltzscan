#!/usr/bin/env python3
"""
Parser for foldbench.mapped.out files from DOMAIN MAPPER v3.0.2
This script parses domain mapping results with metadata and domain annotations.
"""

import pandas as pd
from typing import Dict, List, Union
import re
from dataclasses import dataclass
import io
from Bio import SeqIO
import numpy as np
from pathlib import Path


@dataclass
class DomainEntry:
    """Represents a single domain mapping entry."""
    accession: str
    e_value: float
    residue_range: str
    property: str
    architecture: str
    x_group: str
    t_group: str
    f_group: str
    f_id: str


@dataclass
class DomainMapperMetadata:
    """Represents metadata from the domain mapper output header."""
    version: str
    execution_date: str
    input_hmm: str
    output_file: str
    intra_domain_gap: int
    inter_domain_gap: int
    overlap: int
    e_value_cutoff: float
    total_proteins: int
    total_domains: int
    nc_domains: int
    cp_domains: int
    is_domains: int


def parse_foldbench_mapped(file_path: str) -> Dict[str, Union[DomainMapperMetadata, List[DomainEntry], pd.DataFrame]]:
    """
    Parse a foldbench.mapped.out file and return structured data.
    
    Args:
        file_path: Path to the foldbench.mapped.out file
        
    Returns:
        Dictionary containing:
            - 'metadata': DomainMapperMetadata object with header information
            - 'domains': List of DomainEntry objects with domain annotations
            - 'raw_data': Raw pandas DataFrame of the domain data
            
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid or cannot be parsed
    """
    
    metadata = DomainMapperMetadata(
        version="",
        execution_date="",
        input_hmm="",
        output_file="",
        intra_domain_gap=0,
        inter_domain_gap=0,
        overlap=0,
        e_value_cutoff=0.0,
        total_proteins=0,
        total_domains=0,
        nc_domains=0,
        cp_domains=0,
        is_domains=0
    )
    
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Parse header metadata from commented lines
    for i, line in enumerate(lines):
        line = line.strip()
        
        if not line.startswith('#'):
            # Stop parsing metadata if we reach the data section
            break
            
        if 'DOMAIN MAPPER' in line:
            match = re.search(r'v(\d+\.\d+\.\d+)', line)
            if match:
                metadata.version = match.group(1)
        
        # Handle keys that might have values on the next line
        elif 'Excecuted on:' in line:
            parts = line.split(':', 1)
            if len(parts) > 1 and parts[1].strip().lstrip('#').strip():
                metadata.execution_date = parts[1].strip().lstrip('#').strip()
            elif i + 1 < len(lines):
                metadata.execution_date = lines[i+1].lstrip('#').strip()

        elif 'Input HMM:' in line:
            parts = line.split(':', 1)
            if len(parts) > 1 and parts[1].strip().lstrip('#').strip():
                metadata.input_hmm = parts[1].strip().lstrip('#').strip()
            elif i + 1 < len(lines):
                metadata.input_hmm = lines[i+1].lstrip('#').strip()

        elif 'Output:' in line:
            parts = line.split(':', 1)
            if len(parts) > 1 and parts[1].strip().lstrip('#').strip():
                metadata.output_file = parts[1].strip().lstrip('#').strip()
            elif i + 1 < len(lines):
                metadata.output_file = lines[i+1].lstrip('#').strip()
        
        # Handle options, which are always on a single line
        elif 'Intra domain gap' in line:
            match = re.search(r'Intra domain gap\s*=\s*(\d+)', line)
            if match: metadata.intra_domain_gap = int(match.group(1))
        
        elif 'Inter domain gap' in line:
            match = re.search(r'Inter domain gap\s*=\s*(\d+)', line)
            if match: metadata.inter_domain_gap = int(match.group(1))
        
        elif 'overlap' in line and '=' in line:
            match = re.search(r'overlap\s*=\s*(\d+)', line)
            if match: metadata.overlap = int(match.group(1))
        
        elif 'E-value cutoff' in line:
            match = re.search(r'E-value cutoff\s*=\s*([0-9.eE+-]+)', line)
            if match: metadata.e_value_cutoff = float(match.group(1))

        # Handle Domain Counts / Summary Statistics (supports old and new formats)
        elif 'Total Proteins:' in line: # New format
            match = re.search(r'Total Proteins:\s*(\d+)\s+Total Domains:\s*(\d+)', line)
            if match:
                metadata.total_proteins = int(match.group(1))
                metadata.total_domains = int(match.group(2))
        elif 'Total no of proteins' in line: # Old format support
            match = re.search(r'Total no of proteins\s*=\s*(\d+)', line)
            if match: metadata.total_proteins = int(match.group(1))
        
        elif 'Total no of domains' in line: # Old format support
            match = re.search(r'Total no of domains\s*=\s*(\d+)', line)
            if match: metadata.total_domains = int(match.group(1))
        
        elif 'NC :' in line or 'NC-domains' in line:
            match = re.search(r'(\d+)', line) # Simple number extraction is sufficient here
            if match: metadata.nc_domains = int(match.group(1))
        
        elif 'CP :' in line or 'CP-domains' in line:
            match = re.search(r'(\d+)', line)
            if match: metadata.cp_domains = int(match.group(1))
        
        elif 'IS :' in line or 'IS-domains' in line:
            match = re.search(r'(\d+)', line)
            if match: metadata.is_domains = int(match.group(1))

    # Define column names based on the expected format
    column_names = [
        'accession', 'e_value', 'residue_range', 'property', 
        'architecture', 'x_group', 't_group', 'f_group', 'f_id'
    ]
    
    # Use pandas to read the data table, skipping comment lines
    try:
        data_io = io.StringIO("".join(lines))
        raw_df = pd.read_csv(
            data_io,
            sep='\t',                         # CRITICAL: Use tab as the separator
            comment='#',
            header=None,
            names=column_names,
            usecols=range(len(column_names)), # Safely read only the expected number of columns
            na_filter=False,                  # Treat empty fields as empty strings, not NaN
            engine='python'                   # Python engine is more robust for complex cases
        )
    except Exception as e:
        raise ValueError(f"Failed to parse the data table. Error: {e}")

    # Convert e_value to float, as na_filter=False may cause it to be read as object type
    raw_df['e_value'] = pd.to_numeric(raw_df['e_value'])
    
    # Create a list of DomainEntry objects from the DataFrame
    domains = [DomainEntry(**row) for row in raw_df.to_dict('records')]
    
    return {
        'metadata': metadata,
        'domains': domains,
        'raw_data': raw_df
    }


if __name__ == '__main__':
    # Create a dummy file based on the complex example provided
    file_content = """#===========================================================================================
#  DOMAIN MAPPER v3.0.2
#  Excecuted on:
#               2025-09-11 20:04:10.964888
#  Input HMM: 
#               ./test/foldbench.hmm.out
#  Output:
#               ./test/foldbench.mapped.out
#  Options:
#               Intra domain gap = 30
#               Inter domain gap = 30
#               overlap = 40
#               E-value cutoff = 1.00e-05
#  Domain Counts:
#               Total Proteins:    254         Total Domains:     395
#                                                        NC :  28 (7.09%)
#                                                        CP :   0 (0.00%)
#                                                        IS :  24 (6.08%)
#===========================================================================================
# Accession	E-Value	Residue Range	Property	Architecture	X-group	T-group	F-group	F-id
8cuc-assembly1_F	1.40e-16	1-30		few secondary structure elements	beta-beta-alpha zinc fingers	beta-beta-alpha zinc fingers	zf-C2H2_16	386.1.1.223	
8g8j-assembly1_A	3.50e-42	6-27,88-128,184-235	NC	a+b two layers	Alpha-beta plaits	Adenylyl and guanylyl cyclase catalytic domain-like	IMS	304.48.1.46	
8g8j-assembly1_A	7.80e-25	27-76	IS	a+b three layers	Hypothetical protein Ta1206-like	Hypothetical protein Ta1206-like	IMS_1	850.1.1.4	
"""
    
    dummy_file_path = "foldbench.mapped.out.tmp"
    with open(dummy_file_path, "w") as f:
        f.write(file_content)

    try:
        # Parse the dummy file
        parsed_data = parse_foldbench_mapped(dummy_file_path)
        
        # Print metadata
        print("--- Metadata ---")
        print(parsed_data['metadata'])
        assert parsed_data['metadata'].version == "3.0.2"
        assert parsed_data['metadata'].total_proteins == 254
        assert parsed_data['metadata'].nc_domains == 28
        assert parsed_data['metadata'].input_hmm == "./test/foldbench.hmm.out"
        print("\nMetadata parsed successfully!")
        
        # Print first domain entry
        print("\n--- First Domain Entry ---")
        first_domain = parsed_data['domains'][0]
        print(first_domain)
        assert first_domain.accession == "8cuc-assembly1_F"
        assert first_domain.property == ""
        assert first_domain.architecture == "few secondary structure elements"
        
        # Print second domain entry
        print("\n--- Second Domain Entry ---")
        second_domain = parsed_data['domains'][1]
        print(second_domain)
        assert second_domain.property == "NC"
        assert second_domain.residue_range == "6-27,88-128,184-235"
        
        # Print raw DataFrame
        print("\n--- Raw DataFrame ---")
        print(parsed_data['raw_data'])
        
    except (FileNotFoundError, ValueError) as e:
        print(f"An error occurred: {e}")
    finally:
        # Clean up dummy files
        import os
        if os.path.exists(dummy_file_path):
            os.remove(dummy_file_path)


def get_domain_range(dommap_res_path='foldbench.mapped.out'):
    intervals_dct = {}
    dat = parse_foldbench_mapped(dommap_res_path)['raw_data']
    dat_lst = list(dat.groupby(by='accession'))
    for da in dat_lst:
        id, df = da[0], da[1]
        resi_range_lst = df['residue_range'].tolist()
        if len(resi_range_lst) == 1:
            r_l = resi_range_lst[0].split('-')
            s = int(r_l[0])
            e = int(r_l[-1])
        elif len(resi_range_lst) > 1:
            r_l_s = resi_range_lst[0]
            r_l_e = resi_range_lst[-1]
            s = int(r_l_s.split('-')[0])
            e = int(r_l_e.split('-')[-1])
        intervals_dct[id] = (s, e)
    return intervals_dct


def get_sites_range(fasta_file, clape_results):
    seq_dct = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    clape_res = {}

    for idx, (k, v) in enumerate(seq_dct.items()):
        res =  np.array(list(clape_results[idx]), dtype=np.uint8)  
        assert len(v) == len(res), f'length mismatch for {k}'
        range_intervals = np.where(np.array(res) == 1)[0] + 1
        if range_intervals.shape[0] == 0:
            print(k)
            clape_res[k] = (0, len(v))
        else:
            clape_res[k] = (range_intervals[0], range_intervals[-1])
        # clape_res[k] = (range_intervals[0], range_intervals[-1])
    return clape_res


def read_disopred3_res(res_dir):
    res_dir = Path(res_dir)
    res_files = list(res_dir.glob('*.pbdat'))
    res_dct = {}
    for res_file in res_files:
        res_dct[res_file.stem] = _read_disopred3_res(res_file)
    return res_dct


def _read_disopred3_res(file):
    ids, resis, res_lst = [],[],[]
    with open(file) as f:
        r = [i.strip() for i in f.readlines()]
        r = [i for i in r if '#' not in i]
        r = [i.split(' ') for i in r]
    for j in r:
        ids.append(j[0])
        resis.append(j[1])
        res_lst.append(j[2])
    df = pd.DataFrame([ids, res_lst, resis],index=[0,1,2]).T
    df.columns = ['cif_id', 'res', 'residue']
    return df

def filter_res_2_intervals(res_df):
    intervl_dct = {}
    for name, v in res_df.items():
        v = v[v['res'] != '-']['cif_id'].tolist()
        v = list(map(int, v))
        interval = (v[0], v[-1])
        names = name.split('_')
        n = names[0] + '-' + names[1] + '_' + names[2]
        intervl_dct[n] = interval
    return intervl_dct

def _read_iupred2_res(res_file, dis_threshold=0.5):
    ids, resis, res_lst = [],[],[]
    with open(res_file) as f:
        r = [i.strip() for i in f.readlines()]
        r = [i for i in r if '#' not in i]
        r = [i.split('\t') for i in r]
    for j in r:
        ids.append(j[0])
        resis.append(j[1])
        res_lst.append(j[2])
    df = pd.DataFrame([ids, res_lst, resis],index=[0,1,2]).T
    df.columns = ['cif_id', 'res', 'residue']
    df = df[df['res'].astype(float) <= dis_threshold]
    resi_range = df['cif_id'].values.astype(int)
    s = np.min(resi_range)
    e = np.max(resi_range)
    return (s, e)


def read_iupred2_res(res_dir, dis_threshold):
    res_dir = Path(res_dir)
    res_files = list(res_dir.glob('*.txt'))
    res_dct = {}
    for res_file in res_files:
        res_dct[res_file.stem.split('_iupred')[0]] = _read_iupred2_res(res_file, dis_threshold)
    return res_dct