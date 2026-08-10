import json
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from boltzscan import tfdna
from boltzscan.pwmmap.dbd import dbd_crop_interval
from boltzscan.utils.boltz_input import (
    build_promoter_candidates_from_fimo,
    collapse_model_inputs,
    write_dbd_cropped_model_inputs,
    write_fimo_model_yamls,
    write_full_model_inputs,
)


def _promoter_candidates():
    return pd.DataFrame([
        {
            'candidate_id': 'candidate_p1', 'tf_name': 'TF1', 'tf_seq': 'M' * 324,
            'sequence_name': 'gene_promoter_1', 'motif_id': 'motif_a',
            'start': 10, 'stop': 13, 'strand': '+', 'pvalue': 1e-7, 'score': 20.0,
            'extracted_motif_seq': 'ATGC',
        },
        {
            # Same physical duplex on a different promoter: retain this edge,
            # but reuse the TF1 structural calculation through canonical dsDNA.
            'candidate_id': 'candidate_p2', 'tf_name': 'TF1', 'tf_seq': 'M' * 324,
            'sequence_name': 'gene_promoter_2', 'motif_id': 'motif_a',
            'start': 25, 'stop': 28, 'strand': '-', 'pvalue': 1e-6, 'score': 10.0,
            'extracted_motif_seq': 'GCAT',
        },
        {
            # Same dsDNA but a different TF is a different structural input.
            'candidate_id': 'candidate_p3', 'tf_name': 'TF2', 'tf_seq': 'M' * 324,
            'sequence_name': 'gene_promoter_1', 'motif_id': 'motif_b',
            'start': 40, 'stop': 43, 'strand': '+', 'pvalue': 1e-5, 'score': 8.0,
            'extracted_motif_seq': 'ATGC',
        },
    ])


def test_promoter_edges_are_retained_while_exact_tf_duplex_inputs_reuse(tmp_path):
    candidates, model_inputs, mapping = collapse_model_inputs(_promoter_candidates())

    assert len(candidates) == 3
    assert len(mapping) == 3
    assert len(model_inputs) == 2
    tf1_rows = mapping[mapping.tf_name == 'TF1']
    assert len(tf1_rows) == 2
    assert tf1_rows.model_id.nunique() == 1
    assert tf1_rows.sequence_name.nunique() == 2
    assert mapping.groupby('model_id').is_model_representative.sum().eq(1).all()
    assert mapping.loc[mapping.tf_name == 'TF1', 'model_id'].iloc[0] != (
        mapping.loc[mapping.tf_name == 'TF2', 'model_id'].iloc[0]
    )

    full_input_dir = write_full_model_inputs(
        model_inputs, tmp_path / 'inputs' / 'full'
    )
    assert len(list(full_input_dir.glob('*.yaml'))) == 2
    spec = yaml.safe_load(next(full_input_dir.glob('*.yaml')).read_text())
    assert [entry['protein']['id'] for entry in spec['sequences'] if 'protein' in entry] == [['A']]
    assert next(entry['protein']['msa'] for entry in spec['sequences'] if 'protein' in entry) == 'empty'
    assert [entry['dna']['id'] for entry in spec['sequences'] if 'dna' in entry] == [['B'], ['C']]


def test_fimo_builder_keeps_identical_dna_on_different_promoters(tmp_path):
    (tmp_path / 'tf.fa').write_text('>TF1\n' + 'M' * 324 + '\n')
    (tmp_path / 'promoters.fa').write_text(
        '>gene_promoter_1\nGGATGCAA\n>gene_promoter_2\nGGATGCAA\n'
    )
    (tmp_path / 'tf2pwms.json').write_text(json.dumps({'TF1': ['motif_a']}))
    fimo = pd.DataFrame([
        {
            'sequence_name': 'gene_promoter_1', 'motif_id': 'motif_a',
            'start': 3, 'stop': 6, 'strand': '+', 'pvalue': 1e-7, 'score': 20.0,
        },
        {
            'sequence_name': 'gene_promoter_2', 'motif_id': 'motif_a',
            'start': 3, 'stop': 6, 'strand': '+', 'pvalue': 1e-6, 'score': 10.0,
        },
    ])
    fimo.to_csv(tmp_path / 'fimo.csv', index=False)

    promoter_candidates = build_promoter_candidates_from_fimo(
        tmp_path / 'fimo.csv',
        tmp_path / 'tf.fa',
        tmp_path / 'tf2pwms.json',
        tmp_path / 'promoters.fa',
        dna_flank=0,
    )
    _, model_inputs, mapping = collapse_model_inputs(promoter_candidates)
    assert len(promoter_candidates) == len(mapping) == 2
    assert set(mapping.sequence_name) == {'gene_promoter_1', 'gene_promoter_2'}
    assert len(model_inputs) == 1


def test_dbd_crop_interval_uses_union_and_one_based_inclusive_coordinates():
    records = [
        ('TF1', 'PF03106', 171, 229),
        ('TF1', 'PF99999', 180, 200),
    ]
    assert dbd_crop_interval(records, 'TF1', protein_length=324, flank=20) == (151, 249)


def test_crop_writer_records_rin_style_dbd_plus_20_coordinates(tmp_path):
    candidates = _promoter_candidates().query("tf_name == 'TF1'")
    _, model_inputs, _ = collapse_model_inputs(candidates)
    fields = ['-'] * 21
    fields[0] = 'TF1'
    fields[4] = 'PF03106.1'
    fields[19] = '171'
    fields[20] = '229'
    domtbl = tmp_path / 'pfam.domtbl'
    domtbl.write_text(' '.join(fields) + '\n')

    summary = write_dbd_cropped_model_inputs(
        model_inputs,
        domtbl_path=domtbl,
        crop_input_dir=tmp_path / 'inputs' / 'crop20',
        crop_manifest_path=tmp_path / 'candidates' / 'crop_manifest.csv',
        flank=20,
    )
    manifest = pd.read_csv(summary.crop_manifest)
    assert set(manifest.crop_start_1based_inclusive) == {151}
    assert set(manifest.crop_stop_1based_inclusive) == {249}
    assert set(manifest.crop_length) == {99}
    crop_spec = yaml.safe_load(next(summary.crop_input_dir.glob('*.yaml')).read_text())
    protein = next(entry['protein']['sequence'] for entry in crop_spec['sequences'] if 'protein' in entry)
    assert len(protein) == 99


def test_fimo2yaml_writes_one_msa_backed_yaml_per_model(tmp_path):
    model_table = tmp_path / 'filtered_model_scan_level_res.csv'
    pd.DataFrame([{
        'model_id': 'model_1',
        'tf_name': 'TF.1',
        'tf_seq': 'MPEPT',
        'canonical_dsdna': 'ATGC',
    }]).to_csv(model_table, index=False)
    msa = tmp_path / 'msa' / 'TF_1' / '0.a3m'
    msa.parent.mkdir(parents=True)
    msa.write_text('>query_0\nMPEPT\n>hit\nMPE-T\n')
    fields = ['-'] * 21
    fields[0] = 'TF.1'
    fields[4] = 'PF03106.1'
    fields[19] = '2'
    fields[20] = '4'
    domtbl = tmp_path / 'pfam.domtbl'
    domtbl.write_text(' '.join(fields) + '\n')

    full_dir, crop_summary = write_fimo_model_yamls(
        model_table,
        tmp_path / 'msa',
        domtbl,
        tmp_path / 'inputs',
        crop=0,
    )

    full_yaml = full_dir / 'model_1.yaml'
    spec = yaml.safe_load(full_yaml.read_text())
    protein = next(entry['protein'] for entry in spec['sequences'] if 'protein' in entry)
    dna = [entry['dna']['sequence'] for entry in spec['sequences'] if 'dna' in entry]
    assert protein == {'id': ['A'], 'sequence': 'MPEPT', 'msa': str(msa.resolve())}
    assert dna == ['ATGC', 'GCAT']

    crop_yaml = crop_summary.crop_input_dir / 'model_1.yaml'
    crop_spec = yaml.safe_load(crop_yaml.read_text())
    crop_protein = next(
        entry['protein'] for entry in crop_spec['sequences'] if 'protein' in entry
    )
    crop_msa = next((crop_summary.crop_input_dir / 'msa').glob('*.a3m'))
    assert crop_protein == {
        'id': ['A'], 'sequence': 'PEP', 'msa': str(crop_msa.resolve()),
    }
    assert crop_msa.read_text().startswith('>query_0\nPEP\n')
    manifest = pd.read_csv(crop_summary.crop_manifest)
    assert manifest.loc[0, 'crop_msa_path'] == str(crop_msa.resolve())
    from boltzscan.predict.runners import preflight_prediction_inputs
    assert preflight_prediction_inputs(crop_summary.crop_input_dir, 'boltz2') == 'MSA'


def test_fimo2yaml_without_crop_writes_only_full_inputs(tmp_path):
    model_table = tmp_path / 'filtered_model_scan_level_res.csv'
    pd.DataFrame([{
        'model_id': 'model_1',
        'tf_name': 'TF1',
        'tf_seq': 'MPEPT',
        'canonical_dsdna': 'ATGC',
    }]).to_csv(model_table, index=False)
    msa = tmp_path / 'msa' / 'TF1' / '0.a3m'
    msa.parent.mkdir(parents=True)
    msa.write_text('>query_0\nMPEPT\n')

    full_dir, crop_summary = write_fimo_model_yamls(
        model_table,
        tmp_path / 'msa',
        None,
        tmp_path / 'inputs',
    )

    assert (full_dir / 'model_1.yaml').is_file()
    assert crop_summary is None
    assert {path.name for path in (tmp_path / 'inputs').iterdir()} == {'full'}


@pytest.mark.parametrize(('crop', 'prediction_arms', 'washed_arms'), [
    (None, ('full',), ('full',)),
    (20, ('crop20',), ('interface', 'crop20')),
])
def test_boltz_stage_runs_only_requested_arms(
    tmp_path, monkeypatch, crop, prediction_arms, washed_arms
):
    run = tmp_path / 'run'
    msa = tmp_path / 'query.a3m'
    msa.write_text('>query\nMPEPT\n')
    prepared_arms = ('full',) if crop is None else ('interface',)
    for arm in prepared_arms:
        input_dir = run / 'inputs' / arm
        input_dir.mkdir(parents=True)
        (input_dir / 'model.yaml').write_text(yaml.safe_dump({
            'version': 1,
            'sequences': [
                {'protein': {'id': ['A'], 'sequence': 'MPEPT', 'msa': str(msa)}},
                {'dna': {'id': ['B'], 'sequence': 'ACGT'}},
                {'dna': {'id': ['C'], 'sequence': 'ACGT'}},
            ],
        }, sort_keys=False))

    calls = []

    def fake_run_boltz(**kwargs):
        calls.append(kwargs)
        native = (
            kwargs['out_dir']
            / f"boltz_results_{kwargs['input_dir'].name}"
            / 'predictions'
        )
        native.mkdir(parents=True)
        (native / 'result.txt').write_text('native')
        return 0

    import boltzscan.predict.runners as runners
    import boltzscan.interface as interface_module

    monkeypatch.setattr(runners, 'run_boltz', fake_run_boltz)
    if crop is not None:
        def fake_find_interface(*_args, **_kwargs):
            source = run / 'inputs' / 'interface' / 'model.yaml'
            target = run / 'inputs' / 'crop20'
            target.mkdir(parents=True)
            (target / 'model.yaml').write_text(source.read_text())
            native = run / 'boltz_results_interface' / 'predictions'
            native.mkdir(parents=True)
            (native / 'result.txt').write_text('native')
            return SimpleNamespace(
                n_tfs=1,
                n_model_inputs=1,
                boundaries=run / 'inputs' / 'interface_boundaries.csv',
            )

        monkeypatch.setattr(interface_module, 'run_interface_stage', fake_find_interface)
    tfdna._run_predict_stage(
        out_dir=run,
        crop=crop,
        model='boltz2',
        seed=42,
    )

    assert {call['input_dir'].name for call in calls} == set(prediction_arms)
    assert all(call['out_dir'] == run for call in calls)
    assert all(call['model'] == 'boltz2' for call in calls)
    assert all(call['seed'] == 42 for call in calls)
    for arm in washed_arms:
        assert (run / f'boltz_results_{arm}' / 'predictions' / 'result.txt').is_file()
        assert (run / 'boltz2_prediction' / arm).is_symlink()
    manifest = json.loads((run / 'boltz2_prediction' / 'wash_manifest.json').read_text())
    assert set(manifest['arms']) == set(washed_arms)
    assert manifest['source_policy'] == 'preserve'


@pytest.mark.parametrize('model', ['boltz1', 'boltz2', 'boltz1_ode', 'boltz2_ode'])
def test_every_boltz_model_prepare_generates_msa_and_keeps_logs_at_run_root(
    tmp_path, monkeypatch, model
):
    run_dir = tmp_path / 'named_run'
    calls = {}

    monkeypatch.setattr(
        tfdna,
        '_resolve_pwm_inputs',
        lambda **_kwargs: (
            tmp_path / 'motifs', tmp_path / 'tf2pwms.json', tmp_path / 'pfam.domtbl',
        ),
    )

    from boltzscan.fimocistarget import runner as fimo_runner
    from boltzscan import msa as msa_module
    from boltzscan import interface as interface_module

    def fake_scan(**kwargs):
        calls['scan'] = kwargs
        model_level = run_dir / 'scan' / 'filtered_model_scan_level_res.csv'
        model_level.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([{
            'model_id': 'model_1',
            'tf_name': 'TF1',
            'tf_seq': 'MPEPT',
            'canonical_dsdna': 'ATGC',
        }]).to_csv(model_level, index=False)
        return SimpleNamespace(
            model_level=model_level,
            n_promoter_candidates=3,
            n_model_inputs=2,
        )

    def fake_msa(fasta, output, **_kwargs):
        calls['msa'] = (fasta, output)
        return SimpleNamespace(requested=1, completed=1, skipped=0)

    def fake_interface(model_inputs, output, msa_dir=None):
        calls['interface'] = (model_inputs, output, msa_dir)
        output.mkdir(parents=True, exist_ok=True)
        (output / 'model_1.yaml').write_text('test')
        return output

    monkeypatch.setattr(fimo_runner, 'run_fimo_scan', fake_scan)
    monkeypatch.setattr(msa_module, 'run_msa', fake_msa)
    monkeypatch.setattr(interface_module, 'write_interface_localizer_inputs', fake_interface)

    tfdna._run_prepare_stage(
        tf_fasta=tmp_path / 'tf.fasta',
        promoters=tmp_path / 'promoters.fasta',
        out_dir=run_dir,
        refs=tmp_path / 'refs',
        exclude_species=(),
        collapse_clusters=True,
        pvalue=1e-5,
        overlap_thresh=0.0,
        dna_flank=5,
        crop=20,
        model=model,
        reuse=False,
    )

    assert calls['msa'] == (tmp_path / 'tf.fasta', run_dir / 'msa')
    assert calls['interface'][1] == run_dir / 'inputs' / 'interface'
    assert calls['interface'][2] == run_dir / 'msa'
    assert not (run_dir / 'logs').exists()
    assert {
        path.name for path in run_dir.glob('*.log')
    } == {'map-pwm.log', 'fimo-scan.log', 'msa.log', 'fimo2yaml.log'}


def test_run_score_stage_uses_only_available_prediction_arms(tmp_path, monkeypatch):
    observed = {}

    def fake_score(**kwargs):
        observed.update(kwargs)
        return 'boltz2', ('crop20',)

    monkeypatch.setattr(tfdna, 'score_tf_dna_run', fake_score)
    result = tfdna._run_score_stage(
        out_dir=tmp_path / 'run',
        crop=20,
        processes=4,
        model='boltz2',
    )

    assert result == ('boltz2', ('crop20',))
    assert observed == {
        'out_dir': tmp_path / 'run',
        'processes': 4,
        'model': 'boltz2',
    }


def test_run_restores_model_crop_and_seed_from_run_config(tmp_path, monkeypatch):
    run_dir = tmp_path / 'named_run'
    tf_fasta = tmp_path / 'tf.fasta'
    promoters = tmp_path / 'promoters.fasta'
    tf_fasta.write_text('>TF1\nMPEPT\n')
    promoters.write_text('>gene1\nACGT\n')
    observed = {}

    monkeypatch.setattr(tfdna, '_run_prepare_stage', lambda **_kwargs: None)
    monkeypatch.setattr(
        tfdna,
        '_run_predict_stage',
        lambda **kwargs: observed.update(kwargs),
    )

    tfdna.run_tf_dna_workflow(
        tf_fasta=tf_fasta,
        promoters=promoters,
        out_dir=run_dir,
        model='boltz2',
        crop=20,
        seed=7,
        stage='prepare',
    )
    tfdna.run_tf_dna_workflow(
        tf_fasta=None,
        promoters=None,
        out_dir=run_dir,
        stage='predict',
    )

    assert json.loads((run_dir / 'run_config.json').read_text()) == {
        'schema_version': 1,
        'model': 'boltz2',
        'crop': 20,
        'seed': 7,
    }
    assert observed == {
        'out_dir': run_dir,
        'crop': 20,
        'model': 'boltz2',
        'seed': 7,
    }


@pytest.mark.parametrize('arm', ['full', 'crop20'])
def test_score_projection_restores_every_promoter_edge(tmp_path, monkeypatch, arm):
    run_dir = tmp_path / 'run'
    scan_dir = run_dir / 'scan'
    scan_dir.mkdir(parents=True)
    (run_dir / 'esmfold2_prediction' / arm).mkdir(parents=True)
    pd.DataFrame([{
        'model_id': 'model_1', 'tf_name': 'TF1', 'tf_seq': 'MPEPTIDE',
        'canonical_dsdna': 'ATGC', 'representative_candidate_id': 'candidate_1',
        'n_promoter_candidates': 2, 'n_distinct_promoters': 2,
    }]).to_csv(scan_dir / 'filtered_model_scan_level_res.csv', index=False)
    pd.DataFrame([
        {
            'candidate_id': 'candidate_1', 'model_id': 'model_1', 'tf_name': 'TF1',
            'sequence_name': 'promoter_gene_1', 'motif_id': 'motif_a', 'start': 1,
            'stop': 4, 'strand': '+', 'pvalue': 1e-6, 'score': 10.0,
        },
        {
            'candidate_id': 'candidate_2', 'model_id': 'model_1', 'tf_name': 'TF1',
            'sequence_name': 'promoter_gene_2', 'motif_id': 'motif_a', 'start': 5,
            'stop': 8, 'strand': '+', 'pvalue': 1e-5, 'score': 9.0,
        },
    ]).to_csv(scan_dir / 'filtered_promoter_scan_level_res.csv', index=False)

    def fake_score_ipsae_table(*, score_file, output, **_kwargs):
        table = pd.read_csv(score_file)
        table['ipsae'] = 0.77
        table['ipsae_iptm'] = 0.66
        table.to_csv(output, index=False, compression='gzip')
        return SimpleNamespace(scored_predictions=1, total_predictions=1, warnings=[])

    import boltzscan.utils.ipsae_score as ipsae_score

    monkeypatch.setattr(ipsae_score, 'score_ipsae_table', fake_score_ipsae_table)
    model, arms = tfdna.score_tf_dna_run(out_dir=run_dir, processes=1)

    assert (model, arms) == ('esmfold2', (arm,))
    result = pd.read_csv(run_dir / 'results' / f'{arm}_promoter_level_ipsae.csv.gz')
    assert len(result) == 2
    assert set(result.TG) == {'promoter_gene_1', 'promoter_gene_2'}
    assert set(result.ipsae) == {0.77}
