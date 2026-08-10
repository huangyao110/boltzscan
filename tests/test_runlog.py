from types import SimpleNamespace

import pytest

from boltzscan.cli import _command_log_path
from boltzscan.runlog import CommandLog


def test_command_log_records_options_output_and_completion(tmp_path, capsys):
    log_path = tmp_path / 'step.log'
    with CommandLog(
        log_path,
        ['boltzscan', 'map-pwm', '--no-pwm-cluster'],
        SimpleNamespace(command='map-pwm', no_pwm_cluster=True),
    ):
        print('mapped 3 TFs')

    text = log_path.read_text()
    assert '===BoltzScan===' in text
    assert 'Version: ' in text
    assert 'Command: boltzscan map-pwm --no-pwm-cluster' in text
    assert 'no_pwm_cluster = True' in text
    assert 'mapped 3 TFs' in text
    assert 'Status:  OK' in text
    assert 'mapped 3 TFs' in capsys.readouterr().out


def test_command_log_records_failure(tmp_path):
    log_path = tmp_path / 'failed.log'
    with pytest.raises(ValueError, match='bad input'):
        with CommandLog(
            log_path,
            ['boltzscan', 'fimo-scan'],
            SimpleNamespace(command='fimo-scan'),
        ):
            raise ValueError('bad input')

    text = log_path.read_text()
    assert 'Status:  FAILED' in text
    assert 'Error:   bad input' in text


def test_repeated_command_overwrites_previous_log(tmp_path):
    log_path = tmp_path / 'step.log'
    with CommandLog(log_path, ['boltzscan', 'first'], SimpleNamespace(command='first')):
        print('old run marker')
    with CommandLog(log_path, ['boltzscan', 'second'], SimpleNamespace(command='second')):
        print('new run marker')

    text = log_path.read_text()
    assert 'old run marker' not in text
    assert 'new run marker' in text
    assert text.count('BoltzScan') == 1


@pytest.mark.parametrize(
    ('args', 'expected'),
    [
        (SimpleNamespace(command='run', run='RUN'), 'RUN/run.log'),
        (SimpleNamespace(command='msa', output_dir='RUN/msa'), 'RUN/msa.log'),
        (SimpleNamespace(command='promoter', output='RUN/promoters'), 'RUN/promoter.log'),
        (SimpleNamespace(command='find-tf', output='RUN/tf'), 'RUN/find-tf.log'),
        (SimpleNamespace(command='map-pwm', output='RUN/pwm'), 'RUN/map-pwm.log'),
        (SimpleNamespace(command='fimo-scan', output='RUN/scan'), 'RUN/fimo-scan.log'),
        (SimpleNamespace(command='fimo2yaml', output='RUN/inputs'), 'RUN/fimo2yaml.log'),
        (SimpleNamespace(command='find-interface', run='RUN'), 'RUN/find-interface.log'),
        (SimpleNamespace(command='predict', output='RUN'), 'RUN/predict.log'),
        (SimpleNamespace(command='ipsae', run='RUN'), 'RUN/ipsae.log'),
        (
            SimpleNamespace(
                command='valid', output='RUN/results/validation.csv', res_dir='unused',
            ),
            'RUN/valid.log',
        ),
    ],
)
def test_workflow_logs_resolve_to_run_root(args, expected):
    assert str(_command_log_path(args)) == expected
