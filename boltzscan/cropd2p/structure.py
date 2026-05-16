from Bio.PDB import PDBParser, PDBIO, Select, Selection, NeighborSearch, MMCIFParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.SeqUtils import seq1

from pathlib import Path
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb
import numpy as np
import logging


def read_file_to_biostruc(struct_path, 
                            structure_name="structure"):
    """
    Given the path to a pdb file and the name of the structure, returns a biopython
    structure object.
    """
    parser = PDBParser()
    cif_parser = MMCIFParser()
    struct_path = Path(struct_path)
    if struct_path.suffix == '.pdb':
        structure = parser.get_structure(structure_name, str(struct_path))
    elif struct_path.suffix == '.cif':
        structure = cif_parser.get_structure(structure_name, str(struct_path))

    return structure


def struc_to_seq(structure):
    """
    Takes in a biopython structure object and returns the amino acid sequence as a
    string.
    """

    chains = {
        chain.id: seq1("".join(residue.resname for residue in chain))
        for chain in structure.get_chains()
    }
    if len(chains) > 1:
        msg = (
            "This function is designed for AF2 or Colabfold-generated structures with"
            " a single chain. The input has multiple chains!"
        )
        raise ValueError(msg)

    for chain, seq in chains.items():
        return seq


def structure_to_pLDDT(structure, format="d"):
    """
    Takes in a biopython structure object and returns pLDDT.

    If format is specified as 'd', will return a dictionary of structure pos:pLDDT,
    where pos is the 1-indexed residue position.

    If format is specified as 'l', will return a list with numeric pLDDTs.
    """

    if format not in set(["l", "d"]):
        msg = f"The 'format' parameter must be 'd' or 'l'. Received {format}."
        raise ValueError(msg)

    def get_bfactor_from_residue(residue):
        for atom in residue:
            return atom.get_bfactor()

    pLDDTs = dict()
    for residue in structure.get_residues():

        # 1-indexed position
        pos = residue.get_full_id()[3][1]

        # Get pLDDT, which is stored under bfactor of each atom
        pLDDT = get_bfactor_from_residue(residue)

        pLDDTs[pos] = pLDDT

    if format == "d":
        return pLDDTs
    elif format == "l":
        return list(pLDDTs.values())


def write_structure_to_pdb(structure, path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(path)


def write_structure_subset(structure, residues_to_keep, outfile):
    """
    Writes a pdb outfile that is just a subset of the input structure from the
    residue start to the residue end. Note that the pdb file is 1-indexed, so the
    coordinates are expected to be 1-indexed as well.

    - structure: biopython structure object
    - residues_to_keep: some iterable/list of numbers, where each number is the position
        of one of the residues to keep.
    """

    class ResSelect(Select):
        def accept_residue(self, res):
            if res.id[1] in residues_to_keep:
                return True
            else:
                return False

    io = PDBIO()
    io.set_structure(structure)
    io.save(outfile, ResSelect())


def struc_rebase(structure):
    """
    Input is a biopython structure object. This function renumbers all residues such
    that the first residue starts at 1 and all residues are sequential - so, it takes
    out number gaps.

    Note that this function acts in-place. It doesn't return anything - it just edits
    the structure that is input.
    """
    for i, residue in enumerate(structure.get_residues()):
        i += 1
        res_id = list(residue.id)
        res_id[1] = i
        residue.id = tuple(res_id)


def compare_structures(structure1, structure2):
    """
    This makes sure that all residues are the same in the two input structures. The
    inputs should be biopython structure objects. This is used for testing.
    """
    s1_residues = [r for r in structure1.get_residues()]
    s2_residues = [r for r in structure2.get_residues()]

    assert len(s1_residues) == len(s2_residues)
    for i in range(len(s1_residues)):
        r1 = s1_residues[i]
        r2 = s2_residues[i]
        assert r1 == r2


def extract_chains(structure, chains):
    """
    This function makes a new structure object containing the chain or chains that
    are specified.

    structure - a structure object
    chains - an iterable of chain IDs you want to keep in the new structure object
    """

    # Create a new structure and model
    new_structure = Structure("new_structure")
    new_model = Model(0)
    new_structure.add(new_model)

    # Get the chain from the original structure
    for model in structure:
        for chain in model:
            if chain.get_id() in chains:
                # Create a new chain and copy residues from the original chain
                new_chain = Chain(chain.get_id())
                for residue in chain:
                    new_chain.add(residue.copy())
                # Add the new chain to the new model
                new_model.add(new_chain)
                break

    # Check that all desired chains are present in the final structure
    final_chains = [chain.get_id() for chain in new_structure.get_chains()]
    for chain in chains:
        if chain not in final_chains:
            msg = (
                f"Cannot find the chain, {chain}. The input 'chains' iterable was "
                f"{chains}. Make sure this chain is in the input structure."
            )
            raise ValueError(msg)

    return new_structure


def read_file_to_atomarray(struct_file):
    struct_file = Path(struct_file)
    if struct_file.suffix == '.cif':
        f = pdbx.CIFFile.read(str(struct_file))
        struct = pdbx.get_structure(f, model=1)
    elif struct_file.suffix == '.pdb':
        f = pdb.PDBFile.read(str(struct_file))
        struct = pdb.get_structure(f, model=1)
    return struct, f


def change_chain_id(atom_array):
    chain_id = np.unique(atom_array.chain_id)
    assert len(chain_id) == 1, '.....'
    new_chain_id = 'B'
    new_chain = []
    for atom in atom_array:
        atom.chain_id = new_chain_id
        new_chain.append(atom)
    atom_array = struc.array(new_chain)
    return atom_array


def atom_struct_write(atom_array, out_file):
    out_file = Path(out_file)
    if out_file.suffix == '.cif':
        w_p = pdbx.CIFFile()
        pdbx.set_structure(w_p, atom_array)
        w_p.write(file=out_file)
    elif out_file.suffix == '.pdb':
        w_p = pdb.PDBFile()
        try:
            pdb.set_structure(w_p, atom_array)
        except:
            atom_array = change_chain_id(atom_array)
            pdb.set_structure(w_p, atom_array)
        w_p.write(file=out_file)


def atom_struct_2_seq(cif_file=None, atom_array=None, chain_ids=None):
    if cif_file:
        atom_array, _ = read_file_to_atomarray(cif_file)
    else:
        atom_array = atom_array
    atom_array = atom_array[~atom_array.hetero]
    # atom_array = atom_array[atom_array.chain_id == chain_id]
    if chain_ids is None:
        chain_ids = np.unique(atom_array.chain_id)
    seq_dct = {chain_id: struc.to_sequence(atom_array[atom_array.chain_id == chain_id])[0][0] for chain_id in chain_ids}
    return seq_dct


def write_atom_structure_subset(
    structure: struc.AtomArray, 
    interval,
    chain_id,
    out_file
) -> struc.AtomArray:
    chain_a = structure[structure.chain_id == chain_id]
    chain_a = chain_a[~chain_a.hetero]
    domain = []
    interval_range = range(interval[0], interval[-1]+1)
    for site in interval_range:
        atoms = chain_a[chain_a.res_id == site]
        domain.extend(atoms)
    if out_file:
        atom_struct_write(struc.array(domain), out_file=out_file)
        print(f"成功生成新的CIF文件: {out_file}")
    return struc.array(domain)


def find_interface_atom_residues(
    structure_file: str,
    chain_id1: str,
    chain_id2: str,
    distance_cutoff: float = 2,
) -> dict:
    filepath = Path(structure_file)
    if not filepath.is_file():
        logging.error(f"Structure file not found: {filepath}")
        return {}
    
    atom_array,_ = read_file_to_atomarray(struct_file=structure_file)

    # --- 2. Isolate the two chains of interest ---
    chain1_atoms = atom_array[atom_array.chain_id == chain_id1]
    chain1_atoms = chain1_atoms[~chain1_atoms.hetero]
    chain2_atoms = atom_array[atom_array.chain_id == chain_id2]
    chain2_atoms = chain2_atoms[~chain2_atoms.hetero]

    if chain1_atoms.array_length() == 0:
        logging.error(f"Chain '{chain_id1}' not found in the structure.")
        return {}
    if chain2_atoms.array_length() == 0:
        logging.error(f"Chain '{chain_id2}' not found in the structure.")
        return {}

    # --- 3. Find contacts using CellList ---
    logging.info(f"Using CellList to find contacts with a {distance_cutoff} Å cutoff.")
    
    # Initialize the CellList with atoms from the first chain.
    # The cell_size must be >= distance_cutoff for correctness.
    cell_list = struc.CellList(chain1_atoms, cell_size=distance_cutoff)
    
    # Query the CellList to find which atoms in chain1 are close to atoms in chain2.
    # The result is an array where each row corresponds to an atom in chain2.
    # The columns contain the indices of contacting atoms from chain1.
    # A value of -1 indicates no contact in that slot.
    contacts = cell_list.get_atoms(chain2_atoms.coord, radius=distance_cutoff)

    # --- 4. Identify the interacting residues for both chains ---
    
    # For chain 2: Find which atoms in chain2 had at least one contact.
    # These are the rows in the 'contacts' array that are not all -1.
    contact_indices_in_chain2 = np.where((contacts != -1).any(axis=1))[0]

    # For chain 1: Find which atoms in chain1 were contacted.
    # These are all the non-(-1) values in the 'contacts' array.
    contact_indices_in_chain1 = contacts[contacts != -1]
    
    if contact_indices_in_chain1.size == 0:
        logging.warning("No interactions found between the chains within the cutoff distance.")
        return {
            f'chain_{chain_id1}_residues': [], f'chain_{chain_id1}_range': None,
            f'chain_{chain_id2}_residues': [], f'chain_{chain_id2}_range': None,
        }

    # Get the unique residue IDs for the contacting atoms of each chain.
    # Using auth_res_id is robust as it refers to the original file's numbering.
    interface_res_ids1 = np.unique(chain1_atoms[contact_indices_in_chain1].res_id)
    interface_res_ids2 = np.unique(chain2_atoms[contact_indices_in_chain2].res_id)
    
    logging.info(f"Found {len(interface_res_ids1)} interface residues in chain '{chain_id1}' "
                 f"and {len(interface_res_ids2)} interface residues in chain '{chain_id2}'.")

    # --- 5. Format the output ---
    residues_list1 = sorted(list(interface_res_ids1))
    residues_list2 = sorted(list(interface_res_ids2))

    result = {
        f'{chain_id1}_resi': residues_list1,
        f'{chain_id1}': (min(residues_list1), max(residues_list1)) if residues_list1 else None,
        f'{chain_id2}_resi': residues_list2,
        f'{chain_id2}': (min(residues_list2), max(residues_list2)) if residues_list2 else None,
    }

    return result

def atomarray_reresi(struct_file, outpath):
    path_struct = Path(struct_file)
    outpath = Path(outpath)
    def _res_single_chain(chain):
        old_resid = np.unique(chain.res_id)
        num_resi = len(np.unique(chain.res_id))
        cor_resid = np.array(range(1, num_resi+1))
        old_new_resid_dct = {o:n for o,n in zip(old_resid, cor_resid)}
        reresi_chain = []
        for atom in chain:
            o_atom_id = atom.res_id
            atom.res_id = old_new_resid_dct[o_atom_id]
            # print(atom.red_id)
            reresi_chain.append(atom)
        # print(reresi_chain)
        return struc.array(reresi_chain)
    struct, f = read_file_to_atomarray(struct_file=struct_file)
    chain_ids = np.unique(struct.chain_id)
    chains = []
    for chain_id in chain_ids:
        new_chain = _res_single_chain(struct[struct.chain_id == chain_id])
        chains.append(new_chain)
    new_struct = struc.concatenate(chains)
    if outpath:
        struct_name = path_struct.stem
        if path_struct.suffix == '.cif':
            pdbx.set_structure(f, new_struct)
        elif path_struct.suffix == '.pdb':
            pdb.set_structure(f, new_struct)
        f.write(outpath.joinpath(f'{struct_name}_reresi.cif'))


def crop_assembly(chian_interval_dct, struct_file, outpath):
    struct, f = read_file_to_atomarray(struct_file=struct_file)
    chain_ids = np.unique(struct.chain_id)
    need_chain_ids = list(chian_interval_dct.keys())
    for chain_id in need_chain_ids:
        assert chain_id in chain_ids, f'.......the need_chain_ids:{chain_id} not in chain_ids:{chain_ids}'
    chains = []
    for chain_id in list(chian_interval_dct.keys()):
        crop_chain = write_atom_structure_subset(
            structure=struct[struct.chain_id == chain_id],
            interval=chian_interval_dct[chain_id],
            chain_id=chain_id,
            out_file=None)
        chains.append(crop_chain)
    new_struct = struc.concatenate(chains)
    if outpath:
        outpath = Path(outpath)
        outpath.mkdir(parents=True, exist_ok=True)
        struct_name = Path(struct_file).stem
        if Path(struct_file).suffix == '.cif':
            pdbx.set_structure(f, new_struct)
            f.write(outpath.joinpath(f'{struct_name}_crop.cif'))
        elif Path(struct_file).suffix == '.pdb':
            pdb.set_structure(f, new_struct)
            f.write(outpath.joinpath(f'{struct_name}_crop.pdb'))
    return new_struct
        
        

def pdb_cif_convert_template(sturct_file):
    import subprocess
    sturct_file = Path(sturct_file)
    if sturct_file.suffix == '.pdb':
        model = 1
        out_file = sturct_file.with_suffix('.cif')
    elif sturct_file.suffix == '.cif':
        model = 2
        out_file = sturct_file.with_suffix('.pdb')

    cmd = f"/home/linuxbrew/.linuxbrew/bin/maxit -input {sturct_file} -output {out_file} -o {model}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while predicting for {sturct_file}: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)