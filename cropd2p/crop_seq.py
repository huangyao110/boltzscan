import subprocess
import os
import logging
from pathlib import Path

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def _run_hmmer(fasta_file, output_file, hmm_file):
    cmd = f'hmmscan -o  {output_file} {hmm_file} {fasta_file}'
    try:
        subprocess.run(cmd, shell=True, check=True)
        print("HMMER search completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmscan: {e}")
    except FileNotFoundError:
        print("Error: hmmscan not found")

def run_dommap(fasta_file, 
               output_file, 
               dom_def=None, 
               intra_gap=30, 
               inter_gap=30, 
               overlap=40, 
               frac_overlap=0.7, 
               eval_cutoff=1e-5, 
               update=False):
    """
    运行dommap.py脚本的包装函数
    
    参数:
        input_file: hmmscan输入文件路径
        output_file: 输出文件路径
        dom_def: ECOD域定义文件路径(可选)
        intra_gap: 域内间隙容忍度(默认30)
        inter_gap: 域间间隙容忍度(默认30)
        overlap: 重叠容忍度(默认40)
        frac_overlap: 分数重叠容忍度(默认0.7)
        eval_cutoff: E值阈值(默认1e-5)
        update: 是否更新ECOD域定义(默认False)
    """
    logger.info("Running dommap.py...")
    logger.info(f"STEP 1: Running hmmscan...")
    output_file = Path(output_file)
    _run_hmmer(fasta_file, output_file, "/home/zlab/boltzscan/cropd2p/cropdock/domaper/ecodf/ecodf.hmm")


    script_path = "cropd2p/cropdock/domaper/src/DomainMapper/dommap.py"
    
    # 构建命令
    dommap_out_file = output_file.with_suffix(".dommap")
    cmd = ["python", script_path]
    cmd.extend(["-f", output_file])
    cmd.extend(["-o", dommap_out_file])
    
    if dom_def is not None:
        cmd.extend(["--dom_def", dom_def])
    
    cmd.extend(["--intra_gap", str(intra_gap)])
    cmd.extend(["--inter_gap", str(inter_gap)])
    cmd.extend(["--overlap", str(overlap)])
    cmd.extend(["--frac_overlap", str(frac_overlap)])
    cmd.extend(["--eval_cutoff", str(eval_cutoff)])
    
    if update:
        cmd.append("--update")
    
    # 运行命令
    try:
        subprocess.run(cmd, check=True)
        print("Domain mapping completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error running dommap: {e}")
    except FileNotFoundError:
        print("Error: dommap.py script not found")
