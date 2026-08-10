from pathlib import Path
import gzip
from Bio import SeqIO
from typing import List, Optional
import pandas as pd
import logging
from tqdm import tqdm
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def calculate_fasta_background(fasta_file: str) -> List[float]:
    """Calculate A,C,G,T background frequencies from a FASTA file."""
    fasta_file = Path(fasta_file)
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    opener = gzip.open if fasta_file.suffix == '.gz' else open
    with opener(fasta_file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            seq = str(record.seq).upper().replace('U', 'T')
            for base in counts:
                counts[base] += seq.count(base)

    total = sum(counts.values())
    if total == 0:
        raise ValueError(f"No A/C/G/T bases found in FASTA: {fasta_file}")
    return [counts[base] / total for base in ['A', 'C', 'G', 'T']]


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

            logger.debug(f"成功将基序 '{motif_id}' 写入文件: {output_file}")

        except IOError as e:
            logger.error(f"写入文件时出错: {str(e)}")
            raise

    except Exception as e:
        logger.error(f"处理PWM时出错: {str(e)}")
        raise


def txt_to_meme(
    input_path: str,
    output_path: str,
    nsites: int = 100,
    background: Optional[List[float]] = None,
    force: bool = False,
) -> List[Path]:
    """
    Convert one PWM .txt file or a directory of .txt files to MEME files.

    The input txt format is the four-column PFM/PWM format used in data/pwms,
    for example columns Pos, A, C, G, T. When input_path is a directory,
    output_path is treated as the output directory.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    if background is None:
        background = [0.25, 0.25, 0.25, 0.25]

    if input_path.is_file():
        output_file = output_path
        if output_path.suffix != '.meme':
            output_file = output_path / f'{input_path.stem}.meme'
        return [_convert_one_txt_to_meme(input_path, output_file, nsites, background, force)]

    if input_path.is_dir():
        if output_path.suffix == '.meme':
            raise ValueError("When input is a directory, output must be a directory")
        output_path.mkdir(parents=True, exist_ok=True)
        output_files = []
        failed = 0
        for txt_file in tqdm(sorted(input_path.glob('*.txt')), desc='txt2meme'):
            output_file = output_path / f'{txt_file.stem}.meme'
            try:
                output_files.append(_convert_one_txt_to_meme(txt_file, output_file, nsites, background, force))
            except Exception as exc:
                failed += 1
                logger.warning(f"跳过 {txt_file}: {exc}")
        if failed:
            logger.warning(f"跳过 {failed} 个无法转换的 PWM txt 文件")
        return output_files

    raise FileNotFoundError(f"Input path does not exist: {input_path}")


def _convert_one_txt_to_meme(txt_file, output_file, nsites, background, force):
    output_file = Path(output_file)
    if output_file.exists() and not force:
        return output_file

    motif_info = _read_pwm_txt(txt_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    pwm_to_meme(
        pwm=motif_info['motif'],
        output_file=output_file,
        motif_id=motif_info['name'],
        nsites=nsites,
        background=background,
    )
    return output_file


def _read_pwm_txt(txt_file):
    txt_file = Path(txt_file)
    text = txt_file.read_text()
    if text.lstrip().startswith('MEME version'):
        return _read_meme_txt_matrix(txt_file, text)

    df = pd.read_csv(txt_file, sep=r'\s+', engine='python')
    df.columns = [str(col).strip().upper() for col in df.columns]
    required = ['A', 'C', 'G', 'T']
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"missing PWM columns: {', '.join(missing)}")
    df = df.loc[:, required].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(how='all')
    if df.empty:
        raise ValueError("no PWM rows")
    if df.isna().any().any():
        raise ValueError("PWM contains non-numeric values")
    return {
        'motif': {base: df[base].tolist() for base in required},
        'name': Path(txt_file).stem,
    }


def _read_meme_txt_matrix(txt_file, text):
    matrix = []
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith('letter-probability matrix'):
            for matrix_line in lines[i + 1:]:
                matrix_line = matrix_line.strip()
                if not matrix_line:
                    break
                parts = matrix_line.split()
                if len(parts) < 4:
                    break
                try:
                    matrix.append([float(v) for v in parts[:4]])
                except ValueError:
                    break
            break

    if not matrix:
        raise ValueError("no MEME letter-probability matrix rows")
    df = pd.DataFrame(matrix, columns=['A', 'C', 'G', 'T'])
    return {
        'motif': {base: df[base].tolist() for base in ['A', 'C', 'G', 'T']},
        'name': Path(txt_file).stem,
    }
