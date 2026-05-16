from pathlib import Path
import json


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
