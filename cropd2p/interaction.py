import argparse
import os
import sys
import csv
import warnings
from Bio.PDB import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SeqUtils import seq1
from Bio.Data import PDBData
import numpy as np
from scipy.spatial import KDTree
from pathlib import Path

# --- 常量定义 ---
DISTANCE_THRESHOLD = 5.0  # 交互距离阈值（单位：埃）
DNA_RESIDUES = {"DA", "DC", "DG", "DT", "DI"}
PROTEIN_RESIDUES = set(PDBData.protein_letters_3to1.keys())

# --- 辅助函数 ---

def get_chain_info(chain):
    """
    分析一个链，判断其类型（蛋白质/DNA）并提取序列。
    """
    dna_count = 0
    protein_count = 0
    sequence_map = {}

    dna_3to1 = {
        "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DI": "I"
    }

    for residue in chain:
        if residue.id[0] != ' ':
            continue
            
        res_id = residue.get_id()[1]
        res_name = residue.get_resname().strip()

        if res_name in PROTEIN_RESIDUES:
            protein_count += 1
            try:
                sequence_map[res_id] = seq1(res_name)
            except KeyError:
                sequence_map[res_id] = 'X'
        elif res_name in DNA_RESIDUES:
            dna_count += 1
            sequence_map[res_id] = dna_3to1.get(res_name, 'N')

    if protein_count > dna_count and protein_count > 0:
        chain_type = "protein"
    elif dna_count > protein_count and dna_count > 0:
        chain_type = "dna"
    else:
        chain_type = "other"

    sorted_res_ids = sorted(sequence_map.keys())
    sequence = "".join([sequence_map[res_id] for res_id in sorted_res_ids])

    return chain_type, sequence

def check_interaction(chain1, chain2, distance_threshold):
    """
    使用KDTree高效检查两个链之间是否存在原子间距小于阈值的交互。
    """
    coords1 = np.array([atom.get_coord() for atom in chain1.get_atoms()])
    coords2 = np.array([atom.get_coord() for atom in chain2.get_atoms()])

    if coords1.size == 0 or coords2.size == 0:
        return False

    tree = KDTree(coords1)
    nearby_indices = tree.query_ball_point(coords2, r=distance_threshold)
    
    for indices in nearby_indices:
        if indices:
            return True
            
    return False

def identify_dna_pairs(dna_chains, distance_threshold):
    """
    预先识别哪些DNA链彼此形成双链（或紧密复合物）。
    返回一个字典，将每个DNA链ID映射到一个代表其复合物的唯一ID（例如，排序后的链ID元组）。
    """
    dna_ids = sorted(dna_chains.keys())
    # 使用并查集（Union-Find）或简单的图遍历来分组相互作用的DNA链会更严谨，
    # 这里采用简化的两两配对方法。如果链B和链C互作，它们属于同一组。
    
    # 初始时，每个链自成一组
    complex_map = {dna_id: dna_id for dna_id in dna_ids} 
    
    # 简单的并查集思想来合并有交互的DNA链
    def find_root(d_id):
        if complex_map[d_id] == d_id:
            return d_id
        complex_map[d_id] = find_root(complex_map[d_id]) # 路径压缩
        return complex_map[d_id]

    def union_sets(id1, id2):
        root1 = find_root(id1)
        root2 = find_root(id2)
        if root1 != root2:
            # 简单的按字母顺序合并，保证ID稳定
            if root1 < root2:
                complex_map[root2] = root1
            else:
                complex_map[root1] = root2

    print("  Checking for DNA-DNA pairings (dsDNA detection)...")
    pair_count = 0
    for i in range(len(dna_ids)):
        for j in range(i + 1, len(dna_ids)):
            id1, id2 = dna_ids[i], dna_ids[j]
            # 这里使用相同的距离阈值来判断DNA之间是否形成双链
            if check_interaction(dna_chains[id1]['obj'], dna_chains[id2]['obj'], distance_threshold):
                union_sets(id1, id2)
                pair_count += 1
                # print(f"    - DNA {id1} seems paired with DNA {id2}")

    # 整理最终的复合物ID映射
    final_complex_map = {}
    for dna_id in dna_ids:
        final_complex_map[dna_id] = find_root(dna_id)
        
    if pair_count > 0:
        print(f"  Found {pair_count} DNA-DNA interaction(s). grouped into complexes.")
    
    return final_complex_map

def process_cif_file(cif_path, csv_writer, unique_dsdna=False):
    """
    处理单个CIF文件，找到蛋白质-DNA交互并写入CSV。
    """
    cif_id = os.path.splitext(os.path.basename(cif_path))[0]
    print(f"--- Processing {cif_id} ---")

    parser = MMCIFParser()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            structure = parser.get_structure(cif_id, cif_path)
    except Exception as e:
        print(f"Error parsing {cif_path}: {e}", file=sys.stderr)
        return

    model = structure[0]

    protein_chains = {}
    dna_chains = {}
    for chain in model:
        chain_type, sequence = get_chain_info(chain)
        if not sequence:
            continue
        
        if chain_type == "protein":
            protein_chains[chain.id] = {'obj': chain, 'seq': sequence}
        elif chain_type == "dna":
            dna_chains[chain.id] = {'obj': chain, 'seq': sequence}
    
    print(f"Found {len(protein_chains)} protein chain(s) and {len(dna_chains)} DNA chain(s).")

    if not protein_chains or not dna_chains:
        print("No protein-DNA pairs to compare.")
        return

    # --- 新增步骤：如果需要去冗余，先识别DNA复合物 ---
    dna_complex_map = {}
    if unique_dsdna and len(dna_chains) > 1:
        dna_complex_map = identify_dna_pairs(dna_chains, DISTANCE_THRESHOLD)
    elif unique_dsdna:
        # 只有一条DNA链，它自己就是一个复合物
        for dna_id in dna_chains:
             dna_complex_map[dna_id] = dna_id

    found_interactions = 0
    for prot_id, prot_data in protein_chains.items():
        # 用于记录当前蛋白质已经报告过的DNA复合物ID
        reported_complexes_for_this_protein = set()
        
        for dna_id, dna_data in dna_chains.items():
            
            # --- 核心去冗余逻辑 ---
            if unique_dsdna:
                complex_id = dna_complex_map.get(dna_id)
                # 如果这个DNA链所属的复合物（双链）已经被当前蛋白质报告过了，则跳过
                if complex_id in reported_complexes_for_this_protein:
                    print(f"  [Skipping] Protein {prot_id} vs DNA {dna_id} (already reported its partner in complex {complex_id})")
                    continue

            print(f"Comparing Protein chain {prot_id} with DNA chain {dna_id}...")
            
            if check_interaction(prot_data['obj'], dna_data['obj'], DISTANCE_THRESHOLD):
                print(f"  > Interaction FOUND between Protein {prot_id} and DNA {dna_id}!")
                
                csv_writer.writerow([
                    cif_id,
                    prot_id,
                    prot_data['seq'],
                    dna_id,
                    dna_data['seq']
                ])
                found_interactions += 1
                
                # 如果发现了交互，且开启了去冗余，将该复合物ID标记为已报告
                if unique_dsdna:
                    reported_complexes_for_this_protein.add(complex_id)
                    
            else:
                # print(f"  - No interaction found.")
                pass
    
    print(f"Finished processing {cif_id}. Found {found_interactions} interacting pair(s).")


# --- 主程序入口 ---
def main():
    parser = argparse.ArgumentParser(
        description=f"Extract interacting protein and DNA sequences from CIF files based on a distance threshold of {DISTANCE_THRESHOLD} Å."
    )
    parser.add_argument("--cif_dir", required=True, help='Directory containing input CIF files.')
    parser.add_argument("--output_csv", required=True, help="Path to the output CSV file.")
    # --- 新增参数 ---
    parser.add_argument("--unique_dsdna", action='store_true', 
                        help="If set, for a given protein, only one DNA chain from a double-stranded DNA complex will be reported.")
    
    args = parser.parse_args()

    try:
        cif_dir = Path(args.cif_dir)
        if not cif_dir.is_dir():
            print(f"Error: Provided path '{args.cif_dir}' is not a directory.", file=sys.stderr)
            sys.exit(1)

        with open(args.output_csv, 'w', newline='', encoding='utf-8') as f_out:
            writer = csv.writer(f_out)
            writer.writerow([
                "cif_id", 
                "protein_chain_id", 
                "protein_sequence", 
                "dna_chain_id", 
                "dna_sequence"
            ])
            
            cif_files = list(cif_dir.glob('*.cif'))
            print(f"Found {len(cif_files)} CIF file(s) in '{cif_dir}'.")
            if args.unique_dsdna:
                print("Feature enabled: Unique dsDNA reporting (reporting only one chain per dsDNA complex for each protein).")

            for cif_file in cif_files:
                # 将新参数传递给处理函数
                process_cif_file(cif_file, writer, unique_dsdna=args.unique_dsdna)
            
    except IOError as e:
        print(f"Error writing to output file {args.output_csv}: {e}", file=sys.stderr)
        sys.exit(1)
        
    print(f"\nDone. Results saved to {args.output_csv}")

if __name__ == '__main__':
    main()