import multiprocessing
from pathlib import Path
import subprocess
import logging


def _run_protenix_msa_sig(args):
    """运行单个MSA任务的内部函数"""
    fasta_file, save_dir = args
    fasta_file = Path(fasta_file)
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    
    # 构建输出文件路径
    save_dir = save_dir
    fasta_name = fasta_file.stem
    save_dir = save_dir / fasta_name
    
    # 如果输出文件已存在，则跳过
    if save_dir.joinpath('0.a3m').exists():
        print(f'the {save_dir} file has been existed')
        return None
    try:
        cmd = f"/home/zlab/miniconda3/envs/protenix/bin/protenix msa --input {fasta_file} --out_dir {save_dir}"
        logging.info(f"Running command: {cmd}")
        result = subprocess.run(cmd, shell=True, check=True, 
                              capture_output=True, text=True)
        if result.returncode == 0:
            logging.info(f"Successfully processed {fasta_file}")
            return save_dir
        else:
            logging.error(f"Error processing {fasta_file}: {result.stderr}")
            return None
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to process {fasta_file}: {str(e)}")
        return None
    except Exception as e:
        logging.error(f"Unexpected error processing {fasta_file}: {str(e)}")
        return None

def run_msa_mutil(fasta_files, save_dir, n_cpu=1):
    """
    并行运行多个MSA任务
    
    Args:
        fasta_files: 输入的FASTA文件路径列表或单个路径
        save_dir: 输出目录路径
        n_cpu: 使用的CPU核心数
    
    Returns:
        成功处理的输出文件路径列表
    """
    # 设置日志
    logging.basicConfig(level=logging.INFO)
    
    # 确保输入是列表形式
    if isinstance(fasta_files, (str, Path)):
        fasta_files = [fasta_files]
    
    # 转换为Path对象
    fasta_files = [Path(f) for f in fasta_files]
    save_dir = Path(save_dir)
    
    # 验证输入文件
    valid_files = []
    for fasta_file in fasta_files:
        if fasta_file.exists():
            valid_files.append(fasta_file)
        else:
            logging.warning(f"Input file {fasta_file} does not exist, skipping...")
    
    if not valid_files:
        logging.error("No valid input files found")
        return []
    
    # 创建参数元组列表
    args_list = [(fasta_file, save_dir) for fasta_file in valid_files]
    
    # 创建进程池并执行任务
    pool = multiprocessing.Pool(min(n_cpu, len(valid_files)))
    results = pool.map(_run_protenix_msa_sig, args_list)
    pool.close()
    pool.join()
    
    # 过滤成功的结果
    successful_results = [r for r in results if r is not None]
    
    logging.info(f"Successfully processed {len(successful_results)}/{len(valid_files)} files")
    return successful_results
