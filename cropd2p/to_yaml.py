import yaml
from pathlib import Path


def generate_boltz_yaml_adv(seqs: str, save_path, concat,templates_with_chains: dict[str, str], a3m_path = None) -> str:
    """
    Generates a YAML config where template paths are mapped to specific chain IDs.

    Args:
        seq (str): The amino acid sequence of the protein.
        a3m_path (str): The file path to the a3m MSA file.
        templates_with_chains (dict[str, str]): A dictionary mapping template file paths
                                                 to the chain ID to use from that file.
                                                 e.g., {"path/to/7l4u.cif": "B"}
    """
    
    templates_list = [
        {'cif': str(Path(path)), 'chain_id': chain}
        for chain, path in templates_with_chains.items()
    ] if templates_with_chains else None

    number2alphalet = {1:'A', 2:'B', 3:'C', 4:'D', 5:'E', 6:'F', 7:'G', 8:'H',}

    if templates_list:
        config_dict = {
            'sequences': [
                        {'protein': {'id': f'{number2alphalet[i+1]}', 
                                       'sequence': str(seqs[i]), 
                                       'msa': str(Path(a3m_path[i])) if a3m_path else 'empty'}
                        }  for i in range(len(seqs))
                        ],
            'templates': templates_list,
            'pocket': [
                {'binder': 'B'},
                {'concat': concat}
            ]
        }
    else:
        config_dict = {
            'sequences': [
                        {'protein': {'id': f'{number2alphalet[i+1]}', 
                                     'sequence': str(seqs[i]), 
                                     'msa': str(Path(a3m_path[i])) if a3m_path else 'empty'}
                        }  for i in range(len(seqs))
                        ],
        }       
    with open(save_path, "w") as yaml_file:
        yaml.dump(config_dict, yaml_file, default_flow_style=False)


def add_a3m(yaml_file, a3m_path_dct):
    yaml_file = Path(yaml_file)
    with open(yaml_file, "r") as ya:
        data = yaml.safe_load(ya)
    for sequence in data['sequences']:
        sequence['protein']['msa'] = a3m_path_dct[sequence['protein']['id']]
    with open(f'add_a3m_{yaml_file.name}', 'w') as f:
        yaml.dump(data, f, default_flow_style=False)


def add_templates(yaml_file, templates_dict):
    yaml_file = Path(yaml_file)
    with open(yaml_file, "r") as ya:
        data = yaml.safe_load(ya)
    data['templates'] = [
        {'cif': str(Path(path)), 'chain_id': chain}
        for chain, path in templates_dict.items() if chain in [da['protein']['id'] for da in data['sequences']]
    ]
    with open(f'add_template_{yaml_file.name}', 'w') as f:
        yaml.dump(data, f, default_flow_style=False)    


def add_a3m_template(yaml_file, a3m_path_dct, templates_dict):
    yaml_file = Path(yaml_file)
    add_a3m(yaml_file=yaml_file, a3m_path_dct=a3m_path_dct)
    out_a3m_path = f'add_a3m_{yaml_file.name}'
    add_templates(out_a3m_path, templates_dict=templates_dict)
