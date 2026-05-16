import pandas as pd

def calc_ipsae(cif_file, pae_file, ipsae_script=None):
    import subprocess
    if ipsae_script is None:
        try:
            from pathlib import Path
            script_path = Path(__file__).resolve()
            script_dir = script_path.parent
            ipsae_script = script_dir / "ipsae.py"
            # print(f"Using ipsae script at {ipsae_script}")
        except:
            raise Exception("ipsae script not found")
    cmd = f"python {ipsae_script} {pae_file} {cif_file} 10 10"
    try:
        result = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"ipsae script execution failed: {e}")
    except:
        raise Exception("ipsae script execution failed")

def read_ipsae(ipsae_file):
    df = pd.read_csv(ipsae_file, skipinitialspace=True, sep=' ')  # 跳过列名后的空格
    # 替换数据中的多余空格
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    return df