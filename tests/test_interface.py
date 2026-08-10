import numpy as np
import pandas as pd
from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure

from boltzscan.interface import protein_dna_contacts, select_interface_models


def _atom(name, coordinate, serial, element):
    return Atom.Atom(
        name,
        np.asarray(coordinate, dtype=float),
        1.0,
        1.0,
        ' ',
        name,
        serial,
        element=element,
    )


def test_contact_finder_uses_protein_a_and_both_dna_chains(tmp_path):
    structure = Structure.Structure('prediction')
    model = Model.Model(0)
    structure.add(model)
    protein = Chain.Chain('A')
    dna_b = Chain.Chain('B')
    dna_c = Chain.Chain('C')
    model.add(protein)
    model.add(dna_b)
    model.add(dna_c)
    for index, (name, x) in enumerate((('ALA', 0.0), ('GLY', 20.0)), start=1):
        residue = Residue.Residue((' ', index, ' '), name, ' ')
        residue.add(_atom('CA', (x, 0, 0), index, 'C'))
        protein.add(residue)
    for index, (chain, x) in enumerate(((dna_b, 4.0), (dna_c, 40.0)), start=3):
        residue = Residue.Residue((f'H_DA', 1, ' '), 'DA', ' ')
        residue.add(_atom('P', (x, 0, 0), index, 'P'))
        chain.add(residue)
    path = tmp_path / 'prediction.cif'
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(path))

    assert protein_dna_contacts(path, 'AG', cutoff=5.0) == [1]


def test_one_highest_support_localizer_is_selected_per_tf():
    models = pd.DataFrame([
        {'model_id': 'a', 'tf_name': 'TF1', 'tf_seq': 'AG',
         'canonical_dsdna': 'AAAA', 'n_promoter_candidates': 1},
        {'model_id': 'b', 'tf_name': 'TF1', 'tf_seq': 'AG',
         'canonical_dsdna': 'CCCC', 'n_promoter_candidates': 5},
        {'model_id': 'c', 'tf_name': 'TF2', 'tf_seq': 'AA',
         'canonical_dsdna': 'GGGG', 'n_promoter_candidates': 2},
    ])

    selected = select_interface_models(models)

    assert selected[['tf_name', 'model_id']].values.tolist() == [
        ['TF1', 'b'], ['TF2', 'c'],
    ]
