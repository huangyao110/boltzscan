from pathlib import Path
import json
from Bio import SeqIO, motifs
from typing import List, Optional, Dict, Any
import pandas as pd
import logging
import re
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import os
import shutil
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

def pre_process_msas(unique_msa_map):
    """
    针对唯一的 MSA 列表进行清洗和读取。
    unique_msa_map: List of (msa_path, expected_tf_seq) tuples
    返回: dict {msa_path: clean_protein_seq or None}
    """
    msa_cache = {}
    
    # 这里也可以并行，但考虑到 MSA 文件数量通常远少于 Motif 数量，单线程通常足够
    # 如果 TF 种类特别多 (>1000)，也可以改为并行
    for msa_p, raw_tf_seq in tqdm(unique_msa_map, desc="Pre-processing MSAs"):
        if not msa_p: 
            continue
            
        # 1. 清洗 MSA (IO 操作)
        is_clean_ok, clean_msg = clean_msa_inplace(msa_p)
        if not is_clean_ok:
            print(f"MSA Error [{msa_p}]: {clean_msg}")
            msa_cache[msa_p] = None
            continue

        # 2. 读取序列 (IO 操作)
        try:
            # 只读第一条，获取纯序列
            with open(msa_p, 'r') as f:
                # 简单解析第一行，比 SeqIO.parse 快
                # 假设是标准 FASTA: 第一行 >, 第二行序列(可能多行)
                # 为了稳妥还是用 SeqIO，但只读一条
                record = next(SeqIO.parse(f, "fasta"))
                msa_query_seq = str(record.seq).replace('-', '').upper()
        except Exception as e:
            print(f"Read Error [{msa_p}]: {e}")
            msa_cache[msa_p] = None
            continue

        # 3. 校验序列 (CPU 操作)
        # 清洗 raw_tf_seq
        if raw_tf_seq:
            cleaned_tf_seq = str(raw_tf_seq).upper().strip()
            # 移除末尾特殊字符
            cleaned_tf_seq = re.sub(r'[\-\*X]+$', '', cleaned_tf_seq)
            
            if msa_query_seq != cleaned_tf_seq:
                print(f"Mismatch [{msa_p}]: MSA({len(msa_query_seq)}) != TF({len(cleaned_tf_seq)})")
                msa_cache[msa_p] = None # 标记为无效
                continue
        
        # 校验通过
        msa_cache[msa_p] = msa_query_seq
        
    return msa_cache

def write_single_file(args):
    """
    写入单个文件的原子操作，用于线程池调用
    """
    save_path, content = args
    try:
        with open(save_path, 'w') as f:
            f.write(content)
        return True
    except OSError as e:
        logger.error(f"File write error {save_path}: {e}")
        return False

def to_boltz_fasta(fimo_df, save_dir, tf_seqs_dict=None, tf_msa_dct=None, max_workers=8):
    """
    高性能版 Boltz 输入生成器
    """
    INPUT_SUBDIR_NAME = 'BOLTZ_INPUT'
    save_path = Path(save_dir) / INPUT_SUBDIR_NAME
    save_path.mkdir(exist_ok=True, parents=True)

    df = fimo_df.copy()
    
    print("Step 1: 补全数据...")
    if 'tf_seq' not in df.columns and tf_seqs_dict:
        # 使用 map 比 apply 快
        df['tf_seq'] = df['motif_id'].map(lambda x: str(tf_seqs_dict[x].seq) if x in tf_seqs_dict else None)
    if 'msa_path' not in df.columns and tf_msa_dct:
        df['msa_path'] = df['mid'].map(tf_msa_dct)

    # 过滤掉没有 MSA 路径的数据
    df = df.dropna(subset=['msa_path'])

    print("Step 2: 集中处理 MSA (Deduplication)...")
    # 提取唯一的 (msa_path, tf_seq) 对，避免重复处理同一个文件
    unique_pairs = df[['msa_path', 'tf_seq']].drop_duplicates().values.tolist()
    
    # 获得 {msa_path: valid_protein_seq} 字典
    # 如果校验失败，值为 None
    valid_msa_dict = pre_process_msas(unique_pairs)
    
    # 将清洗后的序列映射回主表
    df['clean_protein_seq'] = df['msa_path'].map(valid_msa_dict)
    
    # 剔除无效 MSA 的行 (None)
    initial_len = len(df)
    df = df.dropna(subset=['clean_protein_seq'])
    print(f"MSA 校验完成: 保留 {len(df)}/{initial_len} 行")

    if len(df) == 0:
        return

    print("Step 3: 向量化准备 DNA 序列...")
    # 反向互补使用 str.translate (上一行已 upper, 所以只需处理 ACGTN)
    dna_seqs = df['extracted_motif_seq'].astype(str).str.upper()
    trans_table = str.maketrans("ACGTN", "TGCAN")
    rev_dna_seqs = dna_seqs.str.translate(trans_table).str[::-1]

    # 文件名替换非法字符 + 拼后缀
    filenames = (
        df['boltz_name'].astype(str).str.replace(r'[\\/*?:"<>|]', "_", regex=True)
        + ".fasta"
    )
    full_paths = [str(save_path / fn) for fn in filenames]

    print(f"Step 4: 并行写入 {len(df)} 个文件 (使用 {max_workers} 线程)...")
    msa_paths = df['msa_path'].values
    prot_seqs = df['clean_protein_seq'].values
    dna_list = dna_seqs.values
    rev_dna_list = rev_dna_seqs.values

    def gen_tasks():
        for path, mp, ps, d, rd in zip(full_paths, msa_paths, prot_seqs, dna_list, rev_dna_list):
            yield (path, f">A|protein|{mp}\n{ps}\n>B|dna\n{d}\n>C|dna\n{rd}\n")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(tqdm(
            executor.map(write_single_file, gen_tasks()),
            total=len(df),
            unit="file"
        ))

    print(f"完成! 成功生成 {sum(results)} 个文件。")


def distribute_files(source_dir, save_boltz_dir, num_subdirs=100):
    """
    将源目录中的文件均匀分配到目标目录的子文件夹中
    
    参数:
        source_dir (str): 源目录路径
        save_boltz_dir (str): 目标目录路径
        num_subdirs (int): 要创建的子文件夹数量，默认为100
    """
    # 确保目标目录存在
    os.makedirs(save_boltz_dir, exist_ok=True)
    
    # 获取源目录中的所有文件
    files = [f for f in os.listdir(source_dir) 
             if os.path.isfile(os.path.join(source_dir, f))]
    
    if not files:
        print(f"在 {source_dir} 中没有找到任何文件")
        return
    
    # 计算每个子文件夹应该包含的文件数
    files_per_dir = len(files) // num_subdirs
    remainder = len(files) % num_subdirs
    
    # 创建子文件夹并分配文件
    file_index = 0
    for i in range(num_subdirs):
        # 创建子文件夹
        subdir_path = os.path.join(save_boltz_dir, f"subdir_{i:03d}")
        os.makedirs(subdir_path, exist_ok=True)
        
        # 确定当前子文件夹应该包含的文件数
        current_files_count = files_per_dir + (1 if i < remainder else 0)
        
        # 将文件移动到子文件夹 (跨文件系统时由 shutil.move 自动 fallback 到 copy+remove)
        for _ in range(current_files_count):
            if file_index < len(files):
                src_file = os.path.join(source_dir, files[file_index])
                dst_file = os.path.join(subdir_path, files[file_index])
                shutil.move(src_file, dst_file)
                file_index += 1

    print(f"成功将 {len(files)} 个文件分配到 {num_subdirs} 个子文件夹中")
