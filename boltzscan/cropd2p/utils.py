from dataclasses import dataclass
from typing import Optional
import pandas as pd
from pathlib import Path
from .structure import read_file_to_atomarray, atom_struct_2_seq


def get_domain_indices(crop_domain_file):
    # 读取包含结构域信息的表格文件
    dat = pd.read_table(crop_domain_file)
    # 获取结构域的索引信息
    indices =  dat['chopping'].values.item()
    print(indices)
    # 检查索引是否为浮点数（可能是空值或NaN）
    if type(indices) == float:
        print('ooooo')
        # 如果是浮点数，则创建一个包含整个残基范围的元组
        sites = [(1, dat['nres'].values.item())]
        # 返回链ID和整个残基范围
        return dat['chain_id'].values.item(), sites
    # 初始化sites列表
    sites = []
    print(indices)
    # 检查索引是否包含逗号，判断是否为多个结构域
    if ',' in indices: # 多个domian的情况
        # 将索引字符串分割成列表
        indices_lst = indices.split(',')
        # 遍历每个结构域区间
        for interval in indices_lst:
            # 提取每个区间的起始和结束位置
            s,e = int(interval.split('-')[0]), int(interval.split('-')[-1])
            print(s,e)
            # 将区间添加到sites列表中
            sites.append((s,e))
    else: # 单个domian的情况
        print(', not in indices')
        # 处理单个结构域区间
        interval = indices
        s,e = int(interval.split('-')[0]), int(interval.split('-')[-1])
        print(s,e)
        # 将单个区间添加到sites列表中
        sites.append((s,e))
    # 返回链ID和结构域区间列表
    return dat['chain_id'].values.item(), sites

def extend_core_seq(core_seq, full_seq, extend_range):
    core_seq_start_pos = full_seq.find(core_seq)
    if core_seq_start_pos == -1:
        raise(f'{core_seq} not find in full_seq')
    core_seq_len = len(core_seq)
    core_seq_end_pos= core_seq_start_pos + core_seq_len
    extend_core_range = (core_seq_start_pos-extend_range, core_seq_end_pos+extend_range)
    return full_seq[extend_core_range[0]: extend_core_range[-1]+1], extend_core_range



def filter_short_assembly(assembly_dir, interaction_df):
    assembly_dir = Path(assembly_dir)
    filter_lst = []
    for row in interaction_df.itertuples():
        aseq_id = row.interface_chain_id_1
        bseq_id = row.interface_chain_id_2
        assembly_id = row.pdb_id
        assembly_files = list(assembly_dir.glob(f'{assembly_id}*'))
        if len(assembly_files) == 0:
            print(f'No assembly file found for {assembly_id}')
            continue
        else:
            assembly_file = next(iter(assembly_files))
            assembly,_ = read_file_to_atomarray(assembly_file)
            seq_size = []
            for chain_id in [aseq_id, bseq_id]:
                try:
                    chain = assembly[assembly.chain_id == chain_id]
                    seq = atom_struct_2_seq(cif_file=None, atom_array=chain)[0][0]
                    seq_size.append(len(seq))
                except Exception as e:
                    print(f'Error processing {assembly_id} {chain_id}: {str(e)}')
            print('-='*20)
            if len(seq_size) > 1:
                print(f'{assembly_id} {aseq_id} {bseq_id} {max(seq_size)}')
                if max(seq_size) > 100:
                    filter_lst.append(pd.DataFrame(row))
    return filter_lst



def write_fasta(seq, save_path):
    save_path = Path(save_path)
    id = save_path.stem
    with open(save_path, 'w') as f:
        f.write(f'>{id}\n')
        f.write(seq)


def select_core_structured_indices(diso_pred, sequence, max_gap=4):
    if len(diso_pred) != len(sequence):
        raise ValueError("diso_pred 和 sequence 的长度必须相同。")
    seq_len = len(sequence)
    n_break_point = seq_len
    gap_count = 0
    for i, pred in enumerate(diso_pred):
        if pred == 0:
            gap_count += 1
        else:
            gap_count = 0
        if gap_count > max_gap:
            n_break_point = i - max_gap
            break
    if n_break_point == seq_len:
         if gap_count > 0 and diso_pred[-1] == 0:
             is_all_disorder = True
             for i in range(seq_len - gap_count, seq_len):
                 if diso_pred[i] != 0:
                     is_all_disorder = False
                     break
             if is_all_disorder:
                 return []
    c_break_point = -1
    gap_count = 0
    for i in range(seq_len - 1, -1, -1):
        pred = diso_pred[i]
        if pred == 0:
            gap_count += 1
        else:
            gap_count = 0
        if gap_count > max_gap:
            c_break_point = i + max_gap
            break
    core_zero_indices = [
        i for i in range(n_break_point, c_break_point + 1)
        if diso_pred[i] == 0
    ]
    intervals = (core_zero_indices[0], core_zero_indices[-1])
    core_seq  = sequence[intervals[0]:intervals[1]+1]
    return intervals, core_seq