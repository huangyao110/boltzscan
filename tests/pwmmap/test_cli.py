from pathlib import Path
from types import SimpleNamespace

import pytest

from boltzscan.cli import (
    BOLTZSCAN_BANNER,
    _build_parser,
    _cmd_map_pwm,
    _terminal_banner,
)
from boltzscan.pwmmap import cluster, leaveout, mapper
from boltzscan.pwmmap.archive import (
    DEFAULT_REFERENCE_RELEASE_SHA256,
    DEFAULT_REFERENCE_RELEASE_URL,
)


def test_map_pwm_parses():
    p = _build_parser()
    ns = p.parse_args(["map-pwm", "-f", "cm.fasta", "-o", "data/pwms/cm_pwms",
                       "--domtbl", "x/pfam.domtbl", "--threshold-mode", "family"])
    assert ns.command == "map-pwm" and ns.proteins == "cm.fasta"
    assert ns.no_pwm_cluster is False


def test_map_pwm_unclustered_mode_is_explicit():
    p = _build_parser()
    ns = p.parse_args([
        'map-pwm', '-f', 'tf.fa', '-o', 'pwms', '--no-pwm-cluster',
    ])
    assert ns.no_pwm_cluster is True


@pytest.mark.parametrize('unclustered', [False, True])
def test_map_pwm_loso_is_orthogonal_to_clustering(tmp_path, monkeypatch, unclustered):
    output = tmp_path / ('raw' if unclustered else 'rep')
    argv = [
        'map-pwm', '-f', 'tf.fa', '-o', str(output),
        '--exclude-species', 'Solanum lycopersicum',
    ]
    if unclustered:
        argv.append('--no-pwm-cluster')
    args = _build_parser().parse_args(argv)
    calls = {}

    def fake_subset(source, destination, *, exclude_species):
        calls['subset'] = (Path(source), Path(destination), tuple(exclude_species))
        return SimpleNamespace(n_retained_dbd_rows=8, n_input_dbd_rows=10)

    def fake_cluster(*, refs_dir, cpu):
        calls['cluster'] = Path(refs_dir)
        return SimpleNamespace(n_clustered=6, n_clusters=3)

    def fake_map_species(**kwargs):
        calls['map'] = kwargs
        return SimpleNamespace(out_dir=output, n_mapped=1, n_species_tfs=1, n_motifs=2)

    monkeypatch.setattr(leaveout, 'create_reference_subset', fake_subset)
    monkeypatch.setattr(cluster, 'cluster_reference_motifs', fake_cluster)
    monkeypatch.setattr(mapper, 'map_species', fake_map_species)

    _cmd_map_pwm(args)

    loso_refs = output / 'refs_loso'
    assert calls['subset'] == (
        Path('data/pwms/_refs'), loso_refs, ('Solanum lycopersicum',),
    )
    assert calls['map']['refs_dir'] == loso_refs
    assert calls['map']['collapse_clusters'] is not unclustered
    assert ('cluster' in calls) is not unclustered


def test_map_pwm_help_documents_common_modes(capsys):
    parser = _build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(['map-pwm', '--help'])
    help_text = capsys.readouterr().out
    assert 'Examples:' in help_text
    assert '--exclude-species "Solanum lycopersicum"' in help_text
    assert '--no-pwm-cluster' in help_text
    assert '--threshold-mode global --threshold 0.80' in help_text
    assert 'OUT/tf2pwms.json' in help_text


def test_build_pwm_refs_parses():
    p = _build_parser()
    ns = p.parse_args(["build-pwm-refs", "--refs", "data/pwms/_refs", "-c", "16"])
    assert ns.command == "build-pwm-refs" and ns.cpu == 16
    assert ns.no_archive is False


def test_install_pwm_refs_uses_built_in_release_by_default():
    parser = _build_parser()
    ns = parser.parse_args(['install-pwm-refs'])
    assert ns.url == DEFAULT_REFERENCE_RELEASE_URL
    assert ns.sha256 == DEFAULT_REFERENCE_RELEASE_SHA256
    assert ns.refs == 'data/pwms/_refs'


def test_install_pwm_refs_accepts_release_override():
    parser = _build_parser()
    ns = parser.parse_args([
        'install-pwm-refs', '--url', 'https://example.org/pwm.tar.gz',
        '--sha256', 'a' * 64,
    ])
    assert ns.url == 'https://example.org/pwm.tar.gz'
    assert ns.sha256 == 'a' * 64


def test_cluster_motifs_is_internal_only():
    parser = _build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(['cluster-motifs'])


def test_cpg_islands_is_not_a_public_command():
    parser = _build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(['cpg-islands'])


def test_run_parses_concise_tf_dna_workflow():
    p = _build_parser()
    ns = p.parse_args([
        'run', '--tf', 'tf.fa', '--promoters', 'promoters.fa', '--run', 'run',
        '--exclude-species', 'Rosa chinensis', '--crop', '20', '--stage', 'prepare',
    ])
    assert ns.command == 'run'
    assert ns.crop == 20
    assert ns.stage == 'prepare'
    assert ns.model is None
    assert not hasattr(ns, 'engine')
    assert ns.seed is None
    assert ns.no_pwm_cluster is False
    assert ns.run == 'run'
    assert ns.pfam is None


def test_run_can_resume_from_only_the_run_directory():
    parser = _build_parser()
    ns = parser.parse_args([
        'run', '--run', 'named_run', '--stage', 'predict',
    ])

    assert ns.run == 'named_run'
    assert ns.tf is None
    assert ns.promoters is None
    assert ns.crop is None
    assert ns.model is None


def test_fimo2yaml_only_writes_full_inputs():
    parser = _build_parser()
    ns = parser.parse_args([
        'fimo2yaml', '--fimo', 'models.csv', '--msa-dir', 'msa', '--output', 'inputs',
    ])

    assert not hasattr(ns, 'crop')
    assert not hasattr(ns, 'domtbl')


def test_find_interface_restores_run_settings_by_default():
    ns = _build_parser().parse_args(['find-interface', '--run', 'named_run'])

    assert ns.run == 'named_run'
    assert ns.flank is None
    assert ns.cutoff == 5.0
    assert ns.model is None
    assert ns.seed is None


def test_run_help_exposes_only_production_parameters(capsys):
    parser = _build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(['run', '--help'])
    help_text = capsys.readouterr().out

    for public_option in ('--run', '--refs', '--pfam', '--stage', '--resume', '--model', '--crop'):
        assert public_option in help_text
    for internal_option in (
        '--out', '--motif-dir', '--tf2pwms', '--domtbl', '--cpu',
        '--processes', '--wash-mode', '--stages', '--reuse',
    ):
        assert internal_option not in help_text


@pytest.mark.parametrize('argv', [
    ['cistarg'],
    ['migrate-run'],
    ['run', '--out', 'run'],
    ['run', '--run', 'run', '--stages', 'prepare'],
    ['run', '--run', 'run', '--reuse'],
])
def test_removed_cli_compatibility_is_rejected(argv):
    with pytest.raises(SystemExit):
        _build_parser().parse_args(argv)


def test_root_help_prints_hash_boltzscan_banner(capsys):
    parser = _build_parser()
    parser.print_help()

    assert BOLTZSCAN_BANNER in capsys.readouterr().out


def test_hash_banner_uses_gradient_only_for_interactive_terminal(monkeypatch):
    class InteractiveStream:
        @staticmethod
        def isatty():
            return True

    monkeypatch.delenv('NO_COLOR', raising=False)
    monkeypatch.setenv('TERM', 'xterm-256color')
    colored = _terminal_banner(InteractiveStream())
    assert colored.count('\033[1;38;5;') == 7
    assert colored.count('\033[0m') == 7

    monkeypatch.setenv('NO_COLOR', '1')
    assert _terminal_banner(InteractiveStream()) == BOLTZSCAN_BANNER


def test_run_accepts_reference_store_outside_the_project():
    parser = _build_parser()
    ns = parser.parse_args([
        'run', '--run', 'RUN', '--refs', '/data/shared/boltzscan_pwm_refs',
        '--stage', 'predict', '--model', 'boltz2',
    ])

    assert ns.refs == '/data/shared/boltzscan_pwm_refs'
def test_fimo_scan_has_only_required_simple_inputs():
    p = _build_parser()
    ns = p.parse_args([
        'fimo-scan', '--promoters', 'promoters.fa', '--motif-dir', 'pwm/meme',
        '--tf2pwms', 'pwm/tf2pwms.json', '--tf', 'tf.fa', '--output', 'scan',
    ])
    assert ns.command == 'fimo-scan'
    assert ns.tf_fasta == 'tf.fa'
    assert ns.no_pwm_cluster is False


def test_promoter_is_the_only_public_bed_fasta_extraction_command():
    parser = _build_parser()
    ns = parser.parse_args([
        'promoter', '--gff', 'genes.gff3', '--genome', 'genome.fa',
        '--output', 'promoters', '--format', 'bed',
    ])
    assert ns.command == 'promoter'
    assert ns.format == 'bed'
    with pytest.raises(SystemExit):
        parser.parse_args(['extract-bed'])


def test_hit2fasta_accepts_promoter_scan_table_name():
    parser = _build_parser()
    ns = parser.parse_args([
        'hit2fasta', '--scan-table', 'scan/filtered_promoter_scan_level_res.csv',
        '--tf2pwms', 'pwm/tf2pwms.json', '--protein-fasta', 'tf.fa',
        '--output', 'candidate_tf_proteins.fasta',
    ])
    assert ns.hit_statistics.endswith('filtered_promoter_scan_level_res.csv')


def test_predict_selects_model_directly_and_exposes_only_seed():
    p = _build_parser()
    for model in ('esmfold2', 'boltz1', 'boltz2', 'boltz1_ode', 'boltz2_ode'):
        ns = p.parse_args([
            'predict', '--input-dir', 'inputs/full', '--output', 'run',
            '--model', model, '--seed', '42',
        ])
        assert ns.model == model
        assert ns.seed == 42
        assert not hasattr(ns, 'engine')
        assert not hasattr(ns, 'sampling_steps')
        assert not hasattr(ns, 'preprocessing_threads')
        assert not hasattr(ns, 'num_workers')
