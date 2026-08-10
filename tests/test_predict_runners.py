import json
from pathlib import Path
import sys
from types import ModuleType
from types import SimpleNamespace

import pytest
import yaml

from boltzscan.predict import runners
from boltzscan.predict import _boltz_worker


def _write_prediction_yaml(input_dir, *, protein='MPEPT', msa='empty', dna='ACGT'):
    input_dir.mkdir(parents=True, exist_ok=True)
    path = input_dir / 'model.yaml'
    path.write_text(yaml.safe_dump({
        'version': 1,
        'sequences': [
            {'protein': {'id': ['A'], 'sequence': protein, 'msa': str(msa)}},
            {'dna': {'id': ['B'], 'sequence': dna}},
            {'dna': {'id': ['C'], 'sequence': dna}},
        ],
    }, sort_keys=False))
    return path


def test_preflight_reports_esmfold_no_msa(tmp_path, capsys):
    input_dir = tmp_path / 'inputs'
    _write_prediction_yaml(input_dir)

    assert runners.preflight_prediction_inputs(input_dir, 'esmfold2') == 'no-MSA'
    assert 'MSA mode=no-MSA (0/1)' in capsys.readouterr().out


def test_preflight_requires_msa_for_every_boltz_model(tmp_path):
    input_dir = tmp_path / 'inputs'
    _write_prediction_yaml(input_dir)

    for model in ('boltz1', 'boltz2', 'boltz1_ode', 'boltz2_ode'):
        with pytest.raises(ValueError, match='requires a local MSA'):
            runners.preflight_prediction_inputs(input_dir, model)


def test_preflight_validates_msa_query_against_protein(tmp_path):
    msa = tmp_path / 'query.a3m'
    msa.write_text('>query\nMPeEPT\n>hit\nMP-EPT\n')
    input_dir = tmp_path / 'inputs'
    yaml_path = _write_prediction_yaml(input_dir, msa=msa)

    assert runners.preflight_prediction_inputs(input_dir, 'boltz2_ode') == 'MSA'
    spec = yaml.safe_load(yaml_path.read_text())
    spec['sequences'][0]['protein']['sequence'] = 'MPEPA'
    yaml_path.write_text(yaml.safe_dump(spec, sort_keys=False))
    with pytest.raises(ValueError, match='Protein/MSA query mismatch'):
        runners.preflight_prediction_inputs(input_dir, 'boltz2_ode')


def test_preflight_rejects_illegal_polymer_characters(tmp_path):
    input_dir = tmp_path / 'inputs'
    _write_prediction_yaml(input_dir, dna='ACGU')

    with pytest.raises(ValueError, match='Illegal dna character'):
        runners.preflight_prediction_inputs(input_dir, 'esmfold2')


def test_preflight_rejects_noncomplementary_dna(tmp_path):
    input_dir = tmp_path / 'inputs'
    _write_prediction_yaml(input_dir, dna='AACG')

    with pytest.raises(ValueError, match='not reverse complements'):
        runners.preflight_prediction_inputs(input_dir, 'esmfold2')


@pytest.mark.parametrize(
    ('public_model', 'native_model', 'patched_family'),
    [
        ('boltz1_ode', 'boltz1', 'boltz1'),
        ('boltz2_ode', 'boltz2', 'boltz2'),
    ],
)
def test_ode_worker_patches_only_the_selected_base_family(
    monkeypatch, public_model, native_model, patched_family
):
    base1 = type('Base1Params', (), {})
    base2 = type('Base2Params', (), {})
    fake_main = SimpleNamespace(
        BoltzDiffusionParams=base1,
        Boltz2DiffusionParams=base2,
    )
    fake_boltz = ModuleType('boltz')
    fake_boltz.main = fake_main
    monkeypatch.setitem(sys.modules, 'boltz', fake_boltz)

    configured_main, argv = _boltz_worker._configure_model([
        'predict', 'inputs', '--model', public_model,
    ])

    assert configured_main is fake_main
    assert argv[argv.index('--model') + 1] == native_model
    if patched_family == 'boltz1':
        assert fake_main.BoltzDiffusionParams is _boltz_worker.BoltzODEDiffusionParams
        assert fake_main.Boltz2DiffusionParams is base2
    else:
        assert fake_main.BoltzDiffusionParams is base1
        assert fake_main.Boltz2DiffusionParams is _boltz_worker.BoltzODEDiffusionParams


def test_boltz_prediction_layout_and_parameters(tmp_path, monkeypatch):
    input_dir = tmp_path / 'inputs' / 'full'
    input_dir.mkdir(parents=True)
    output_root = tmp_path / 'run'
    calls = []

    monkeypatch.setattr(
        runners.subprocess,
        'run',
        lambda cmd, **kwargs: calls.append((cmd, kwargs)) or SimpleNamespace(returncode=0),
    )
    monkeypatch.setattr(runners, 'available_cpu_count', lambda: 32)

    result = runners.run_boltz(
        input_dir,
        output_root,
        model='boltz1_ode',
        seed=42,
        python_exe=sys.executable,
    )

    native_root = output_root / 'boltz_results_full'
    assert result == 0
    assert not (output_root / 'boltz2_prediction').exists()
    metadata = json.loads((native_root / 'inference_parameters.json').read_text())
    assert metadata['model'] == 'boltz1_ode'
    assert metadata['output_dir'] == str(native_root / 'predictions')
    assert metadata['parameters']['sampling_steps'] == 2
    assert metadata['parameters']['preprocessing_threads'] == 32
    assert calls[0][0][1].endswith('_boltz_worker.py')
    assert calls[0][0][calls[0][0].index('--model') + 1] == 'boltz1_ode'


def test_boltz_receives_only_yaml_when_input_dir_contains_manifest(tmp_path, monkeypatch):
    input_dir = tmp_path / 'inputs' / 'crop20'
    input_dir.mkdir(parents=True)
    (input_dir / 'model.yaml').write_text('version: 1\n')
    (input_dir / 'crop_manifest.csv').write_text('model_id,crop_start\n')
    observed = {}

    def fake_run(cmd, **_kwargs):
        staged = Path(cmd[3])
        observed['name'] = staged.name
        observed['files'] = sorted(path.name for path in staged.iterdir())
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runners.subprocess, 'run', fake_run)
    runners.run_boltz(
        input_dir,
        tmp_path / 'run',
        model='boltz2',
        python_exe=sys.executable,
    )

    assert observed == {'name': 'crop20', 'files': ['model.yaml']}


@pytest.mark.parametrize('model', ['boltz1', 'boltz2'])
def test_standard_boltz_models_use_native_inference_defaults(tmp_path, monkeypatch, model):
    input_dir = tmp_path / 'inputs' / 'full'
    input_dir.mkdir(parents=True)
    output_root = tmp_path / 'run'
    calls = []

    monkeypatch.setattr(
        runners.subprocess,
        'run',
        lambda cmd, **kwargs: calls.append((cmd, kwargs)) or SimpleNamespace(returncode=0),
    )

    result = runners.run_boltz(
        input_dir,
        output_root,
        model=model,
        python_exe=sys.executable,
    )

    assert result == 0
    command = calls[0][0]
    assert '--sampling_steps' not in command
    assert '--step_scale' not in command
    assert '--seed' not in command
    metadata = json.loads(
        (output_root / 'boltz_results_full' / 'inference_parameters.json').read_text()
    )
    assert metadata['parameters']['configuration_source'] == 'boltz_native_default'
    assert metadata['parameters']['sampling_steps'] is None
    assert metadata['parameters']['sampling_steps_source'] == 'boltz_native_default'
    assert metadata['parameters']['step_scale'] is None
    assert metadata['parameters']['step_scale_source'] == 'boltz_native_default'
    assert metadata['parameters']['seed'] is None
    assert metadata['parameters']['seed_source'] == 'native_default'


@pytest.mark.parametrize('model', ['boltz1_ode', 'boltz2_ode'])
def test_ode_variants_default_to_two_sampling_steps(
    tmp_path, monkeypatch, model
):
    input_dir = tmp_path / 'inputs' / 'full'
    input_dir.mkdir(parents=True)
    output_root = tmp_path / 'run'
    calls = []

    monkeypatch.setattr(
        runners.subprocess,
        'run',
        lambda cmd, **kwargs: calls.append((cmd, kwargs)) or SimpleNamespace(returncode=0),
    )

    result = runners.run_boltz(
        input_dir,
        output_root,
        model=model,
        seed=42,
        python_exe=sys.executable,
    )

    assert result == 0
    command = calls[0][0]
    assert command[command.index('--model') + 1] == model
    assert command[command.index('--sampling_steps') + 1] == '2'
    assert command[command.index('--step_scale') + 1] == '1.0'
    metadata = json.loads(
        (output_root / 'boltz_results_full' / 'inference_parameters.json').read_text()
    )
    assert metadata['parameters']['sampling_steps'] == 2
    assert metadata['parameters']['sampling_steps_source'] == f'boltzscan_{model}'
    assert metadata['parameters']['step_scale'] == 1.0


def test_esmfold_accepts_only_the_shared_seed_override(tmp_path, monkeypatch):
    input_dir = tmp_path / 'inputs' / 'crop20'
    input_dir.mkdir(parents=True)
    native_root = tmp_path / 'run' / 'esmfold_results_crop20'
    calls = []
    weights = tmp_path / 'weights' / 'ESMFold2'
    weights.mkdir(parents=True)
    (weights / 'ccd.pkl').write_bytes(b'ccd')
    esmc_weights = tmp_path / 'weights' / 'ESMC-6B'
    monkeypatch.setenv('BOLTZSCAN_ESMFOLD2_WEIGHTS', str(weights))
    monkeypatch.setenv('BOLTZSCAN_ESMC_WEIGHTS', str(esmc_weights))
    monkeypatch.delenv('ESMCFOLD_CCD_PATH', raising=False)

    monkeypatch.setattr(
        runners.subprocess,
        'run',
        lambda cmd, **kwargs: calls.append((cmd, kwargs)) or SimpleNamespace(returncode=0),
    )

    result = runners.run_esmfold(
        input_dir,
        native_root,
        seed=42,
        python_exe=sys.executable,
    )
    assert result == 0
    command = calls[0][0]
    assert calls[0][1]['env']['PYTHONPATH'].split(':')[0] == str(runners._ESM_SRC)
    assert calls[0][1]['env']['ESMCFOLD_CCD_PATH'] == str(weights / 'ccd.pkl')
    assert command[command.index('--seed') + 1] == '42'
    assert command[command.index('--weights') + 1] == str(weights)
    assert command[command.index('--esmc-weights') + 1] == str(esmc_weights)
    assert '--sampling-steps' not in command
    metadata = json.loads((native_root / 'inference_parameters.json').read_text())
    assert metadata['model'] == 'esmfold2'
    assert metadata['parameters']['seed'] == 42
    assert metadata['parameters']['seed_source'] == 'explicit'
    assert metadata['parameters']['deterministic'] is True
    assert metadata['parameters']['deterministic_source'] == 'explicit_seed'
    assert metadata['parameters']['esmfold2_weights'] == str(weights)
    assert metadata['parameters']['esmc_weights'] == str(esmc_weights)
    assert metadata['parameters']['ccd_path'] == str(weights / 'ccd.pkl')


def test_esmfold_uses_public_local_inference_defaults(tmp_path, monkeypatch):
    input_dir = tmp_path / 'inputs' / 'full'
    input_dir.mkdir(parents=True)
    native_root = tmp_path / 'run' / 'esmfold_results_full'
    calls = []

    monkeypatch.setattr(
        runners.subprocess,
        'run',
        lambda cmd, **kwargs: calls.append((cmd, kwargs)) or SimpleNamespace(returncode=0),
    )

    result = runners.run_esmfold(
        input_dir,
        native_root,
        python_exe=sys.executable,
    )

    assert result == 0
    command = calls[0][0]
    for option in ('--num-loops', '--sampling-steps', '--diffusion-samples', '--seed'):
        assert option not in command
    metadata = json.loads((native_root / 'inference_parameters.json').read_text())
    assert metadata['parameters']['configuration_source'] == 'esmfold2_native_default'
    assert metadata['parameters']['seed'] is None
    assert metadata['parameters']['seed_source'] == 'native_default'
    assert metadata['parameters']['deterministic'] is False
    assert metadata['parameters']['deterministic_source'] == 'native_default'


@pytest.mark.parametrize('model', ['boltz1', 'boltz2'])
def test_standard_boltz_models_accept_seed_without_sampler_overrides(
    tmp_path, monkeypatch, model
):
    calls = []
    monkeypatch.setattr(
        runners.subprocess,
        'run',
        lambda cmd, **kwargs: calls.append(cmd) or SimpleNamespace(returncode=0),
    )
    runners.run_boltz(
        tmp_path / 'inputs',
        tmp_path / 'run',
        model=model,
        seed=7,
        python_exe=sys.executable,
    )
    command = calls[0]
    assert command[command.index('--seed') + 1] == '7'
    assert '--sampling_steps' not in command
    assert '--step_scale' not in command


def test_prediction_root_names_include_inference_method():
    assert runners.prediction_root_name('esmfold2') == 'esmfold2_prediction'
    assert runners.prediction_root_name('boltz1') == 'boltz1_prediction'
    assert runners.prediction_root_name('boltz2') == 'boltz2_prediction'
    assert runners.prediction_root_name('boltz1_ode') == 'boltz1_ode_prediction'
    assert runners.prediction_root_name('boltz2_ode') == 'boltz2_ode_prediction'
