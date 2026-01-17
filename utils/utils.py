# from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
from cropd2p.cropdock.structure import pdb_cif_convert_template
import pandas as pd


# def calc_dockq(pred_struct, native_struct, map_lst=['A', 'B']):

#     pdb_cif_convert_template(pred_struct)

#     model = load_PDB(str(pred_struct.with_suffix('.pdb')))
#     native = load_PDB(str(native_struct))

#     # native:model chain map dictionary for two interfaces
#     chain_map = {map_lst[0]:map_lst[0], map_lst[1]:map_lst[1]}
#     # returns a dictionary containing the results and the total DockQ score
#     RES = run_on_all_native_interfaces(model, native, chain_map=chain_map)
#     return pd.DataFrame(RES[0]).T


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