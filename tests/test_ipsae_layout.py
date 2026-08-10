import pandas as pd

from boltzscan.cli import _build_parser
from boltzscan.utils import ipsae_score
from boltzscan.utils.ipsae_score import _resolve_predictions_dir


def test_ipsae_accepts_clean_per_arm_prediction_root(tmp_path):
    arm_dir = tmp_path / 'predictions' / 'full'
    model_dir = arm_dir / 'model_abc'
    model_dir.mkdir(parents=True)
    (model_dir / 'model_abc_model_0.cif').write_text('data_test')

    assert _resolve_predictions_dir(arm_dir) == arm_dir


def test_ipsae_cli_needs_only_run_directory():
    args = _build_parser().parse_args([
        'ipsae', '--run', 'RUN',
    ])

    assert args.run == 'RUN'
    assert args.model is None


def test_ipsae_standalone_output_stays_inside_prediction_arm(tmp_path, monkeypatch):
    arm_dir = tmp_path / 'predictions' / 'full'
    model_dir = arm_dir / 'model_abc'
    model_dir.mkdir(parents=True)
    (model_dir / 'model_abc.cif').write_text('data_test')
    metrics = {
        key: 0.5 for key in (
            'ipsae', 'ipsae_asym_min', 'boltz_iptm', 'ipsae_iptm', 'iptm_diff',
            'boltz_iptm_global', 'boltz_iptm_asym_min',
            'ipsae_iptm_asym_min', 'pDockQ',
        )
    }
    monkeypatch.setattr(
        ipsae_score,
        '_collect_ipsae_scores',
        lambda *_args, **_kwargs: ({'model_abc': metrics}, [], {'model_abc'}),
    )

    summary = ipsae_score.score_ipsae_table(arm_dir)

    assert summary.output == arm_dir / 'ipsae_scored.csv'
    result = pd.read_csv(summary.output)
    assert result['model_id'].tolist() == ['model_abc']
    assert 'TF' not in result and 'TG' not in result
