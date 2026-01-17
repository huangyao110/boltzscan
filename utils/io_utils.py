from pathlib import Path
import json
from Bio import SeqIO
import numpy as np
from Bio.Seq import Seq
from Bio import motifs
from typing import List, Optional, Dict, Any
import pandas as pd
import logging
from tqdm import tqdm
from cropd2p.cropdock.a3m import read_a3m
from Bio.SeqRecord import SeqRecord
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_data_path():
    """
    Get the path to the data directory.
    """
    return Path(__file__).parent.parent / "data"


def read_fasta(path):
    """
    Read a fasta file and return a dictionary of sequences.
    """
    return {record.id: str(record.seq) for record in SeqIO.parse(path, "fasta")}


def read_fasta_from_dir(fasta_dir):
    """
    Read all fasta files in a directory and return a dictionary of sequences.
    """
    id2seq = {}
    fasta_dir = Path(fasta_dir)
    for fasta_file in fasta_dir.glob("*.fasta"):
        for record in SeqIO.parse(fasta_file, "fasta"):
            id2seq[record.id] = str(record.seq)
    return id2seq


def write_fasta(path, sequences):
    """
    Write a dictionary of sequences to a fasta file.
    """
    with open(path, "w") as f:
        for id, seq in sequences.items():
            f.write(f">{id}\n{seq}\n")
    return


def return_boltz_res_path(boltz_path):
    res_dct = {}
    boltz_path = Path(boltz_path)
    pred_path = list(boltz_path.joinpath('predictions').iterdir())[0]
    try:
        msa_path = [i for i in list(boltz_path.joinpath('msa').iterdir()) if i.is_dir()][0]
        msa_file = msa_path.joinpath('uniref.a3m')
        res_dct['a3m'] = msa_file
    except:
        # print(f'NOTE: Not find the msa file')
        res_dct['a3m'] = ''
    for f in pred_path.iterdir():
        if f.suffix == '.json':
            res_dct['json'] = f
        elif f.suffix == '.cif':
            res_dct['cif'] = f
        elif 'pae' in f.stem:
            res_dct['pae'] = f
        elif f.suffix == '.txt' and 'byres' not in f.stem:
            res_dct['ipsae'] = f
    return res_dct


def return_boltz_res_path_v2(boltz_path):
    res_dct = {}
    boltz_path = Path(boltz_path)
    try:
        msa_path = [i for i in list(boltz_path.joinpath('msa').iterdir()) if i.is_dir()][0]
        msa_file = msa_path.joinpath('uniref.a3m')
        res_dct['a3m'] = msa_file
    except:
        # print(f'NOTE: Not find the msa file')
        res_dct['a3m'] = ''
    for f in boltz_path.iterdir():
        if f.suffix == '.json':
            res_dct['json'] = f
        elif f.suffix == '.cif':
            res_dct['cif'] = f
        elif 'pae' in f.stem:
            res_dct['pae'] = f
        elif f.suffix == '.txt' and 'byres' not in f.stem:
            res_dct['ipsae'] = f
    return res_dct


def load_boltz_json(boltz_json):
    boltz_json = Path(boltz_json)
    with open(boltz_json, 'r') as f:
        data = json.load(f)
    return data


def get_boltz_res(boltz_res_dir, res_type='cif'):
    if res_type == 'iptm':
        res_dct = {}
        for i in boltz_res_dir.iterdir():
            res = return_boltz_res_path_v2(i)['json']
            res_dct[i.name] = load_boltz_json(res)['iptm']
        res_df = pd.DataFrame(res_dct, index=['boltzscan score']).T
        return res_df
    elif res_type == 'cif':
        res_dct = {}
        for i in boltz_res_dir.iterdir():
            res = return_boltz_res_path_v2(i)['cif']
            res_dct[i.name] = res
    elif res_type == 'a3m':
        res_dct = {}
        for i in boltz_res_dir.iterdir():
            res = return_boltz_res_path_v2(i)['a3m']
            res_dct[i.name] = res
    elif res_type == 'pae':
        res_dct = {}
        for i in boltz_res_dir.iterdir():
            res = return_boltz_res_path_v2(i)['pae']
            res_dct[i.name] = res
    elif res_type == 'ipsae':
        res_dct = {}
        for i in boltz_res_dir.iterdir():
            res = return_boltz_res_path_v2(i)['ipsae']
            res_dct[i.name] = res
    return res_dct





def read_motif(motif_path: str) -> Dict[str, Any]:
    """
    从文件中读取基序数据。

    Args:
        motif_path (str): 基序文件路径

    Returns:
        Dict[str, Any]: 包含基序矩阵和名称的字典

    Raises:
        FileNotFoundError: 当文件不存在时
        ValueError: 当文件格式不正确时
    """
    try:
        motif_path = Path(motif_path)
        if not motif_path.exists():
            raise FileNotFoundError(f"基序文件不存在: {motif_path}")

        name = motif_path.stem
        with open(motif_path, 'r') as handle:
            m = motifs.read(handle, 'pfm-four-columns')

        return {'motif': m.pwm, 'name': name}
    except Exception as e:
        logger.error(f"读取基序文件时出错: {str(e)}")
        raise

def pwm_to_meme(
    pwm: List[List[float]],
    output_file: str,
    motif_id: str,
    nsites: int,
    motif_name: Optional[str] = None,
    e_value: float = 0.0,
    meme_version: str = "5.5.0",
    alphabet: str = "ACGT",
    background: List[float] = [0.25, 0.25, 0.25, 0.25],
    url: Optional[str] = None
) -> None:
    """
    将一个PWM/PPM矩阵转换为MEME格式的文件。

    Args:
        pwm (List[List[float]]):
            位置概率矩阵。一个列表的列表，外层列表代表基序的位置，
            内层列表代表在那个位置上A, C, G, T的概率。
            例如: [[0.1, 0.4, 0.4, 0.1], [0.8, 0.1, 0.05, 0.05], ...]
        output_file (str):
            要保存的输出文件名 (例如 'motif.meme')。
        motif_id (str):
            基序的唯一标识符 (例如 'MP00028')。
        nsites (int):
            用于构建此基序的序列位点数量。这是MEME格式的必需字段。
        motif_name (Optional[str]):
            基序的可读名称 (例如 'AT1G09770')。如果未提供，将使用 motif_id。
        e_value (float):
            基序的统计显著性 E-value。默认为 0.0。
        meme_version (str):
            写入文件头的MEME版本号。默认为 '5.5.0'。
        alphabet (str):
            字母表。对于DNA，应为 'ACGT'。默认为 'ACGT'。
        background (List[float]):
            背景碱基频率 [A, C, G, T]。默认为均一分布 [0.25, 0.25, 0.25, 0.25]。
        url (Optional[str]):
            一个可选的URL，提供关于该基序的更多信息。

    Raises:
        ValueError: 当输入参数不合法时
        IOError: 当文件写入失败时
    """
    try:
        if not pwm:
            raise ValueError("PWM 不能为空")
        if nsites <= 0:
            raise ValueError("nsites 必须是正整数")


        # 转换为DataFrame以便处理
        pwm_df = pd.DataFrame([pwm[i] for i in alphabet]).T
        width = len(pwm_df)
        alength = len(alphabet)

        if len(background) != alength:
            raise ValueError(f"背景频率列表的长度 ({len(background)}) 必须与字母表长度 ({alength}) 匹配。")

        # 检查每行的概率和
        for i, row in enumerate(pwm_df.itertuples(index=False)):
            row_sum = sum(row)
            if not (0.99 <= row_sum <= 1.01):
                logger.warning(f"PWM 第 {i+1} 行的概率之和为 {row_sum:.3f}，不接近 1.0")


        # --- 2. 写入文件 ---
        try:
            with open(output_file, 'w') as f:
                # 写入文件头
                f.write(f"MEME version {meme_version}\n\n")
                f.write(f"ALPHABET= {alphabet}\n\n")
                f.write("strands: + -\n\n")

                f.write("Background letter frequencies (from uniform background):\n")
                bg_line = " ".join([f"{char} {freq:.4f}" for char, freq in zip(alphabet, background)])
                f.write(bg_line + "\n\n")

                # 写入基序定义
                f.write(f"MOTIF {motif_id}\n\n")

                # 写入概率矩阵头信息
                f.write(f"letter-probability matrix: alength= {alength} w= {width} nsites= {nsites} E= {e_value}\n")

                # 写入概率矩阵
                for _, row in pwm_df.iterrows():
                    line = "\t".join([f"{prob:0.6f}" for prob in row])
                    f.write(line + "\n")

                # 写入可选的URL
                if url:
                    f.write(f"\nURL {url}\n")

            logger.info(f"成功将基序 '{motif_id}' 写入文件: {output_file}")

        except IOError as e:
            logger.error(f"写入文件时出错: {str(e)}")
            raise

    except Exception as e:
        logger.error(f"处理PWM时出错: {str(e)}")
        raise


def clean_msa_inplace(msa_path):
    """
    1. 检查第一条序列末尾是否为 '-', 'X', '*'。
    2. 若是，检查所有序列该位置是否也是如此。
    3. 同步删除所有序列的该位置字符。
    4. 手动写入文件：确保一个序列占满一行。
    """
    msa_path = Path(msa_path)
    if not msa_path.exists():
        return False, f"MSA file not found: {msa_path}"

    try:
        # 读取 MSA
        records = list(SeqIO.parse(msa_path, "fasta"))
        if not records:
            return False, "MSA is empty"

        modified = False
        while True:
            first_seq_str = str(records[0].seq)
            if not first_seq_str:
                break
            
            last_char = first_seq_str[-1].upper()
            
            if last_char in ['-', 'X', '*']:
                # 严格检查：确保所有序列在这一列都是垃圾字符
                for i, rec in enumerate(records):
                    current_last = str(rec.seq)[-1].upper()
                    if current_last not in ['-', 'X', '*']:
                        raise ValueError(
                            f"Consistency Error: Query ends with '{last_char}', "
                            f"but Record {i} ({rec.id}) has '{current_last}'."
                        )
                
                # 全列删除最后一个字符
                for rec in records:
                    rec.seq = rec.seq[:-1]
                modified = True
                continue 
            else:
                break
        
        # --- 核心修改：手动写入文件，确保一个序列一行 ---
        if modified:
            with open(msa_path, 'w') as f:
                for rec in records:
                    # 写入 ID/描述行
                    f.write(f">{rec.description}\n")
                    # 写入完整序列并换行（不进行 60 字符切割）
                    f.write(f"{str(rec.seq)}\n")
            return True, "MSA cleaned and overwritten (flat format)"
        
        return True, "No cleaning needed"

    except Exception as e:
        return False, str(e)
# ==========================================
# 2. 生成 Boltz FASTA 主函数 (集成严格校验逻辑)
# ==========================================
def to_boltz_fasta(fimo_df, save_dir, tf_seqs_dict=None, tf_msa_dct=None):
    """
    最终生成 Boltz 输入文件。
    """
    INPUT_SUBDIR_NAME = 'BOLTZ_INPUT'
    df = fimo_df.copy()
    
    # 数据补全
    if 'tf_seq' not in df.columns and tf_seqs_dict:
        df['tf_seq'] = df['motif_id'].map(lambda x: str(tf_seqs_dict[x].seq) if x in tf_seqs_dict else None)
    if 'msa_path' not in df.columns and tf_msa_dct:
        df['msa_path'] = df['mid'].map(tf_msa_dct)
    
    save_path = Path(save_dir) / INPUT_SUBDIR_NAME
    save_path.mkdir(exist_ok=True, parents=True)
    
    success_count = 0
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc='Generating Boltz FASTAs'):
        msa_p = row['msa_path']
        tf_seq = row['tf_seq']
        
        # --- 步骤 1: 物理清洗 MSA ---
        is_clean_ok, clean_msg = clean_msa_inplace(msa_p)
        if not is_clean_ok:
            print(f"Skip {row.boltz_name} due to MSA Error: {clean_msg}")
            continue

        # --- 步骤 2: 读取清洗后的 MSA 第一条序列 ---
        try:
            msa_records = list(SeqIO.parse(msa_p, "fasta"))
            # 去除 Gap，获取纯蛋白质序列字符串
            msa_query_seq = str(msa_records[0].seq).replace('-', '').upper()
        except Exception as e:
            print(f"Skip {row.boltz_name}: Failed to read cleaned MSA - {e}")
            continue

        # --- 步骤 3: 同样清洗 tf_seq 并进行严格布尔检查 ---
        # 按照相同逻辑，删去 tf_seq 末尾的 '-', '*', 'X'
        cleaned_tf_seq = str(tf_seq).upper().strip()
        while cleaned_tf_seq and cleaned_tf_seq[-1] in ['-', '*', 'X']:
            cleaned_tf_seq = cleaned_tf_seq[:-1]
        
        # 核心布尔检查：要求完全一致
        if msa_query_seq != cleaned_tf_seq:
            print(f"Skip {row.boltz_name}: Sequence mismatch! "
                  f"MSA_query({len(msa_query_seq)}) != TF_seq({len(cleaned_tf_seq)})")
            continue

        # --- 步骤 4: 准备 DNA 并写入文件 ---
        core_dna = str(row['core_seq']).upper()
        rev_dna = str(Seq(core_dna).reverse_complement())

        safe_fn = re.sub(r'[\\/*?:"<>|]', "_", str(row.boltz_name))
        save_file = save_path / f'{safe_fn}.fasta'
        
        try:
            with open(save_file, 'w') as f:
                # 此时 msa_query_seq 的长度与物理 MSA 文件里的有效列数完全对齐
                # 这将彻底解决 ipsae.py 中的 IndexError
                f.write(f">A|protein|{msa_p}\n{msa_query_seq}\n")
                f.write(f">B|dna\n{core_dna}\n")
                f.write(f">C|dna\n{rev_dna}\n")
            success_count += 1
        except Exception as e:
            print(f"File write error {row.boltz_name}: {e}")

    print(f"\nComplete: {success_count} Boltz FASTA files generated at {save_path}")

def get_core_seq(row, extend_range=5):
    seq = row['promoter_seq']
    s, e = row['start'], row['stop']
    left = min(s, e) - 1 - extend_range
    right = max(s, e) + extend_range
    left = max(left, 0)
    right = min(right, len(seq))
    return seq[left:right]