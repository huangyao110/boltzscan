import json

import pytest

from boltzscan.predict.wash import wash_prediction_outputs


def _native_prediction(run, result_name, arm='model_1'):
    prediction = run / result_name / 'predictions' / arm
    prediction.mkdir(parents=True)
    artifact = prediction / 'result.cif'
    artifact.write_text('structure')
    return prediction.parent, artifact


def test_soft_wash_is_live_and_preserves_native_tree(tmp_path):
    run = tmp_path / 'run'
    source, artifact = _native_prediction(run, 'boltz_results_full')

    summary = wash_prediction_outputs(
        run,
        model='boltz2',
        mode='soft',
    )

    published = run / 'boltz2_prediction' / 'full'
    assert published.is_symlink()
    assert published.resolve() == source.resolve()
    assert artifact.is_file()
    assert (published / 'model_1' / 'result.cif').read_text() == 'structure'
    manifest = json.loads(summary.manifest.read_text())
    assert manifest['mode'] == 'soft'
    assert manifest['source_policy'] == 'preserve'


def test_hard_wash_materializes_view_without_removing_source(tmp_path):
    run = tmp_path / 'run'
    source, artifact = _native_prediction(run, 'boltz_results_full')
    wash_prediction_outputs(run, model='boltz2', mode='soft')

    summary = wash_prediction_outputs(
        run,
        model='boltz2',
        mode='hard',
    )

    published = run / 'boltz2_prediction' / 'full'
    published_artifact = published / 'model_1' / 'result.cif'
    assert published.is_dir()
    assert not published.is_symlink()
    assert artifact.is_file()
    assert published_artifact.samefile(artifact)
    assert summary.hard_linked == 1
    manifest = json.loads(summary.manifest.read_text())
    assert manifest['mode'] == 'hard'
    assert manifest['arms']['full']['source_preserved'] is True


def test_wash_discovers_native_esmfold_layout(tmp_path):
    native_run = tmp_path / 'native_run'
    native, _artifact = _native_prediction(native_run, 'esmfold_results_crop20')
    native_summary = wash_prediction_outputs(native_run, model='esmfold2')
    assert native_summary.arms == ('crop20',)
    assert (native_run / 'esmfold2_prediction' / 'crop20').resolve() == native.resolve()


def test_wash_rejects_a_model_name_that_conflicts_with_native_metadata(tmp_path):
    run = tmp_path / 'run'
    _source, _artifact = _native_prediction(run, 'boltz_results_full')
    (run / 'boltz_results_full' / 'inference_parameters.json').write_text(
        json.dumps({'engine': 'boltz', 'model': 'boltz2'})
    )

    with pytest.raises(ValueError, match='not requested model=boltz1'):
        wash_prediction_outputs(run, model='boltz1')
    assert not (run / 'boltz1_prediction' / 'full').exists()
