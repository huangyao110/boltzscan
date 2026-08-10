"""Find one reusable TF--dsDNA interface and build cropped model inputs."""

from dataclasses import dataclass
import json
from pathlib import Path

import pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch, is_aa
from Bio.SeqUtils import seq1

from boltzscan.msa import read_a3m_query
from boltzscan.predict.runners import (
    engine_for_model,
    native_prediction_dir,
    prediction_root_name,
    run_prediction,
)
from boltzscan.utils.boltz_input import (
    _build_msa_dct,
    _resolve_msa_path,
    write_coordinate_cropped_model_inputs,
    write_full_model_inputs,
)


@dataclass(frozen=True)
class InterfaceSummary:
    boundaries: Path
    crop_input_dir: Path
    crop_manifest: Path
    n_tfs: int
    n_model_inputs: int


def select_interface_models(model_inputs):
    """Select one deterministic, highest-support dsDNA localizer per TF."""
    required = {'model_id', 'tf_name', 'tf_seq', 'canonical_dsdna'}
    missing = sorted(required.difference(model_inputs.columns))
    if missing:
        raise ValueError('Model-level scan table is missing: ' + ', '.join(missing))
    order = ['tf_name']
    ascending = [True]
    if 'n_promoter_candidates' in model_inputs.columns:
        order.append('n_promoter_candidates')
        ascending.append(False)
    order.append('model_id')
    ascending.append(True)
    return (
        model_inputs.sort_values(order, ascending=ascending, kind='mergesort')
        .drop_duplicates('tf_name', keep='first')
        .reset_index(drop=True)
    )


def _tf_msas(model_inputs, msa_dir):
    available = _build_msa_dct(msa_dir)
    result = {}
    for tf_name, group in model_inputs.groupby('tf_name', sort=False):
        sequences = set(group['tf_seq'].astype(str).str.upper())
        if len(sequences) != 1:
            raise ValueError(f'TF {tf_name!r} has inconsistent protein sequences')
        path = _resolve_msa_path(str(tf_name), available)
        if path is None:
            raise FileNotFoundError(f'MSA not found under {msa_dir} for TF {tf_name}')
        path = Path(path).resolve()
        if read_a3m_query(path) != sequences.pop():
            raise ValueError(f'Protein/MSA query mismatch for TF {tf_name}: {path}')
        result[tf_name] = path
    return result


def write_interface_localizer_inputs(model_inputs, input_dir, msa_dir=None):
    """Write exactly one full-length TF--dsDNA YAML per TF."""
    selected = select_interface_models(model_inputs)
    msa_paths = _tf_msas(selected, msa_dir) if msa_dir is not None else None
    input_dir = Path(input_dir)
    expected = {f'{model_id}.yaml' for model_id in selected['model_id'].astype(str)}
    if input_dir.is_dir():
        stale = [path for path in input_dir.glob('*.yaml') if path.name not in expected]
        if stale:
            raise FileExistsError(
                f'{input_dir} contains {len(stale)} stale localizer YAML file(s)'
            )
    write_full_model_inputs(selected, input_dir, msa_paths=msa_paths)
    selected.to_csv(input_dir / 'interface_manifest.csv', index=False)
    return input_dir


def _heavy_atoms(residues):
    return [
        atom
        for residue in residues
        for atom in residue.get_atoms()
        if str(getattr(atom, 'element', '')).strip().upper() != 'H'
    ]


def protein_dna_contacts(cif_path, expected_sequence, cutoff=5.0):
    """Return 1-based protein residues contacting either dsDNA strand."""
    if cutoff <= 0:
        raise ValueError('Interface contact cutoff must be positive')
    structure = MMCIFParser(QUIET=True).get_structure('prediction', str(cif_path))
    model = next(structure.get_models())
    protein = model.child_dict.get('A')
    dna = [model.child_dict.get(chain) for chain in ('B', 'C')]
    if protein is None or any(chain is None for chain in dna):
        raise ValueError('Prediction must contain protein A and dsDNA chains B/C')
    residues = [residue for residue in protein.get_residues() if is_aa(residue, standard=False)]
    observed = ''.join(seq1(residue.resname, undef_code='X') for residue in residues)
    if observed != str(expected_sequence).upper():
        raise ValueError(
            f'Predicted/input protein sequence mismatch: {len(observed)} != '
            f'{len(str(expected_sequence))}'
        )
    dna_atoms = _heavy_atoms(
        residue for chain in dna for residue in chain.get_residues()
    )
    if not dna_atoms:
        raise ValueError('Prediction contains no dsDNA heavy atoms')
    neighbours = NeighborSearch(dna_atoms)
    return [
        index
        for index, residue in enumerate(residues, start=1)
        if any(
            neighbours.search(atom.coord, cutoff, level='A')
            for atom in _heavy_atoms([residue])
        )
    ]


def _prediction_directory(run_dir, model):
    run_dir = Path(run_dir)
    for arm in ('interface', 'full'):
        published = run_dir / prediction_root_name(model) / arm
        if published.is_dir():
            return arm, published
        native = native_prediction_dir(run_dir / 'inputs' / arm, run_dir, model)
        if native.is_dir():
            return arm, native
    raise FileNotFoundError(
        f'No interface/full predictions found in {run_dir}; run full interface localization first'
    )


def _prediction_cif(prediction_dir, model_id):
    directory = Path(prediction_dir) / str(model_id)
    preferred = sorted(directory.glob('*_model_0.cif'))
    candidates = preferred or sorted(directory.glob('*.cif'))
    return candidates[0] if candidates else None


def run_interface_stage(run_dir, flank=None, cutoff=5.0, model=None, seed=None):
    """Run or reuse one full localizer prediction per TF, then build crop inputs."""
    run_dir = Path(run_dir)
    config_path = run_dir / 'run_config.json'
    config = json.loads(config_path.read_text()) if config_path.is_file() else {}
    model = model or config.get('model')
    seed = config.get('seed') if seed is None else seed
    if model is None:
        raise ValueError('Prediction model is unknown; pass --model or use a BoltzScan run')
    input_dir = run_dir / 'inputs' / 'interface'
    if not input_dir.is_dir() or not any(input_dir.glob('*.yaml')):
        raise FileNotFoundError(
            f'Interface localizer inputs not found: {input_dir}. Run prepare with --crop first.'
        )
    table = pd.read_csv(run_dir / 'scan' / 'filtered_model_scan_level_res.csv')
    selected = select_interface_models(table)
    native = native_prediction_dir(input_dir, run_dir, model)
    published = run_dir / prediction_root_name(model) / 'interface'
    existing = published if published.is_dir() else native
    complete = existing.is_dir() and all(
        _prediction_cif(existing, model_id) is not None
        for model_id in selected['model_id']
    )
    if complete:
        print(f'Reusing {len(selected)} completed interface localizer prediction(s)')
    else:
        result = run_prediction(
            input_dir=input_dir,
            out_dir=run_dir,
            model=model,
            seed=seed,
        )
        if result != 0:
            raise RuntimeError(f'{model} interface localization failed with exit code {result}')
    return find_interface_crop(
        run_dir,
        flank=flank,
        cutoff=cutoff,
        model=model,
    )


def find_interface_crop(run_dir, flank=None, cutoff=5.0, model=None):
    """Infer one interface per TF and materialize all cropped TF--DNA YAMLs."""
    run_dir = Path(run_dir)
    config_path = run_dir / 'run_config.json'
    config = json.loads(config_path.read_text()) if config_path.is_file() else {}
    model = model or config.get('model')
    if model is None:
        raise ValueError('Prediction model is unknown; pass --model or use a BoltzScan run')
    flank = config.get('crop') if flank is None else flank
    if flank is None:
        raise ValueError('Crop flank is unknown; pass --flank or prepare RUN with --crop')
    if flank < 0:
        raise ValueError('Crop flank must be non-negative')

    table_path = run_dir / 'scan' / 'filtered_model_scan_level_res.csv'
    if not table_path.is_file():
        raise FileNotFoundError(f'Model-level FIMO table not found: {table_path}')
    model_inputs = pd.read_csv(table_path)
    selected = select_interface_models(model_inputs)
    arm, prediction_dir = _prediction_directory(run_dir, model)

    intervals = {}
    evidence = {}
    rows = []
    failures = []
    for row in selected.itertuples(index=False):
        cif = _prediction_cif(prediction_dir, row.model_id)
        if cif is None:
            failures.append(f'{row.tf_name}: missing prediction for {row.model_id}')
            continue
        try:
            contacts = protein_dna_contacts(cif, row.tf_seq, cutoff=cutoff)
        except ValueError as exc:
            failures.append(f'{row.tf_name}: {exc}')
            continue
        if not contacts:
            failures.append(f'{row.tf_name}: no protein-DNA contact within {cutoff:g} A')
            continue
        start, stop = min(contacts), max(contacts)
        intervals[row.tf_name] = (start, stop)
        evidence[row.tf_name] = cif.resolve()
        rows.append({
            'tf_name': row.tf_name,
            'localizer_model_id': row.model_id,
            'localizer_arm': arm,
            'prediction_cif': str(cif.resolve()),
            'contact_cutoff_angstrom': cutoff,
            'contact_residue_count': len(contacts),
            'contact_start_1based_inclusive': start,
            'contact_stop_1based_inclusive': stop,
            'contact_positions_1based': ';'.join(map(str, contacts)),
        })
    if failures:
        raise ValueError('Interface localization failed: ' + '; '.join(failures[:10]))

    msa_paths = (
        _tf_msas(model_inputs, run_dir / 'msa')
        if engine_for_model(model) == 'boltz'
        else None
    )
    crop_dir = run_dir / 'inputs' / f'crop{flank}'
    summary = write_coordinate_cropped_model_inputs(
        model_inputs,
        intervals_by_tf=intervals,
        crop_input_dir=crop_dir,
        crop_manifest_path=crop_dir / 'crop_manifest.csv',
        flank=flank,
        msa_paths=msa_paths,
        boundary_source='predicted_interface',
        evidence_by_tf=evidence,
    )
    boundaries = run_dir / 'inputs' / 'interface_boundaries.csv'
    pd.DataFrame(rows).sort_values('tf_name', kind='mergesort').to_csv(
        boundaries, index=False
    )
    return InterfaceSummary(
        boundaries=boundaries,
        crop_input_dir=summary.crop_input_dir,
        crop_manifest=summary.crop_manifest,
        n_tfs=len(rows),
        n_model_inputs=summary.n_model_inputs,
    )
