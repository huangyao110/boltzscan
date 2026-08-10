import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from boltzscan.fimocistarget import runner
from boltzscan.fimocistarget.runner import (
    MODEL_SCAN_FILENAME,
    PROMOTER_SCAN_FILENAME,
    RAW_FIMO_FILENAME,
    run_fimo_scan,
)


def _inputs(tmp_path, *, clustered=True, manifest=True):
    promoters = tmp_path / 'promoters.fa'
    promoters.write_text('>gene1\nGGATGCAA\n')
    tf_fasta = tmp_path / 'tf.fa'
    tf_fasta.write_text('>TF1\nMPEPTIDE\n')
    pwm_dir = tmp_path / 'pwm'
    (pwm_dir / 'meme').mkdir(parents=True)
    mapping = pwm_dir / 'tf2pwms.json'
    mapping.write_text(json.dumps({'TF1': ['rep_1']}))
    if manifest:
        (pwm_dir / 'pwm_mapping.json').write_text(json.dumps({
            'reference_dir': str(tmp_path / 'refs'),
            'pwm_clustered': clustered,
        }))
    return promoters, tf_fasta, pwm_dir, mapping


class _FakeScanner:
    def __init__(self, *, output_dir, **_kwargs):
        self.output_dir = output_dir

    def run(self):
        hits = pd.DataFrame([{
            'sequence_name': 'gene1',
            'start': 3,
            'stop': 6,
            'strand': '+',
            'score': 20.0,
            'pvalue': 1e-7,
            'qvalue': 1e-6,
            'motif_id': 'rep_1',
        }])
        hits.to_csv(f'{self.output_dir}/{RAW_FIMO_FILENAME}', index=False)
        return hits


def test_fimo_scan_writes_exactly_three_public_csvs(tmp_path, monkeypatch):
    promoters, tf_fasta, pwm_dir, mapping = _inputs(tmp_path)
    monkeypatch.setattr(runner, '_FimoScanner', _FakeScanner)

    summary = run_fimo_scan(
        target_fasta=promoters,
        motif_dir=pwm_dir / 'meme',
        tf2pwms_path=mapping,
        tf_fasta=tf_fasta,
        output_dir=tmp_path / 'scan',
        dna_flank=0,
    )

    assert summary.pwm_clustered is True
    assert {path.name for path in (tmp_path / 'scan').glob('*.csv')} == {
        RAW_FIMO_FILENAME,
        PROMOTER_SCAN_FILENAME,
        MODEL_SCAN_FILENAME,
    }
    promoter_level = pd.read_csv(summary.promoter_level)
    model_level = pd.read_csv(summary.model_level)
    assert len(promoter_level) == len(model_level) == 1
    assert promoter_level.loc[0, 'model_id'] == model_level.loc[0, 'model_id']


def test_fimo_scan_fails_closed_without_cluster_provenance(tmp_path):
    promoters, tf_fasta, pwm_dir, mapping = _inputs(tmp_path, manifest=False)

    with pytest.raises(ValueError, match='--no-pwm-cluster'):
        run_fimo_scan(
            target_fasta=promoters,
            motif_dir=pwm_dir / 'meme',
            tf2pwms_path=mapping,
            tf_fasta=tf_fasta,
            output_dir=tmp_path / 'scan',
        )


def test_fimo_scan_allows_unclustered_pwm_only_when_explicit(tmp_path, monkeypatch):
    promoters, tf_fasta, pwm_dir, mapping = _inputs(
        tmp_path, clustered=False, manifest=True
    )
    monkeypatch.setattr(runner, '_FimoScanner', _FakeScanner)

    summary = run_fimo_scan(
        target_fasta=promoters,
        motif_dir=pwm_dir / 'meme',
        tf2pwms_path=mapping,
        tf_fasta=tf_fasta,
        output_dir=tmp_path / 'scan',
        dna_flank=0,
        no_pwm_cluster=True,
    )

    assert summary.pwm_clustered is False


def _meme_file(path, motif_id='internal_id'):
    Path(path).write_text(
        'MEME version 5\n\n'
        'ALPHABET= ACGT\n\n'
        f'MOTIF {motif_id}\n'
        'letter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0\n'
        '0.90 0.03 0.03 0.04\n'
        '0.03 0.90 0.03 0.04\n'
        '0.03 0.03 0.90 0.04\n'
        '0.04 0.03 0.03 0.90\n'
    )


def test_external_fimo_adapter_normalizes_schema_and_motif_ids(tmp_path, monkeypatch):
    promoters = tmp_path / 'promoters.fa'
    promoters.write_text('>gene1\nACGTACGT\n')
    motif_dir = tmp_path / 'motifs'
    motif_dir.mkdir()
    _meme_file(motif_dir / 'rep_1.meme')
    output = tmp_path / 'scan'
    output.mkdir()
    commands = []

    monkeypatch.setattr(runner, 'resolve_executable', lambda name: '/managed/bin/fimo')

    def fake_run(command, **_kwargs):
        if command[-1] == '--version':
            return SimpleNamespace(returncode=0, stdout='5.5.7\n', stderr='')
        commands.append(command)
        background = Path(command[command.index('--bfile') + 1]).read_text()
        assert background == (
            'A 0.3\n'
            'C 0.2\n'
            'G 0.2\n'
            'T 0.3\n'
        )
        fimo_output = Path(command[command.index('--oc') + 1])
        fimo_output.mkdir()
        (fimo_output / 'fimo.tsv').write_text(
            'motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\t'
            'score\tp-value\tq-value\tmatched_sequence\n'
            'rep_1\t.\tgene1\t1\t4\t+\t9.5\t1e-6\t2e-6\tACGT\n'
            '# command line\n'
        )
        combined = Path(command[-2]).read_text()
        assert 'MOTIF rep_1\n' in combined
        assert 'MOTIF internal_id\n' not in combined
        return SimpleNamespace(returncode=0, stdout='', stderr='')

    monkeypatch.setattr(runner.subprocess, 'run', fake_run)
    scanner = runner._FimoScanner(
        fasta_file=promoters,
        motif_dir=motif_dir,
        output_dir=output,
        custom_bg={'A': 0.2, 'C': 0.1, 'G': 0.3, 'T': 0.4},
        pvalue_thresh=1e-4,
        motif_filter={'rep_1'},
    )
    hits = scanner.run()

    assert list(hits.columns) == [
        'sequence_name', 'start', 'stop', 'strand',
        'score', 'pvalue', 'qvalue', 'motif_id',
    ]
    assert hits.iloc[0].to_dict() == {
        'sequence_name': 'gene1', 'start': 1, 'stop': 4, 'strand': '+',
        'score': 9.5, 'pvalue': 1e-6, 'qvalue': 2e-6, 'motif_id': 'rep_1',
    }
    assert '--no-pgc' in commands[0]
    assert commands[0][commands[0].index('--motif-pseudo') + 1] == '0.1'
    assert commands[0][commands[0].index('--max-stored-scores') + 1] == '500000'


@pytest.mark.network
@pytest.mark.skipif(
    runner.resolve_executable('fimo') is None,
    reason='requires MEME Suite FIMO or `boltzscan doctor --fix`',
)
def test_real_fimo_roundtrip(tmp_path):
    promoters = tmp_path / 'promoters.fa'
    promoters.write_text('>gene1\nACGTACGT\n')
    motif_dir = tmp_path / 'motifs'
    motif_dir.mkdir()
    _meme_file(motif_dir / 'rep_1.meme')
    output = tmp_path / 'scan'
    output.mkdir()

    scanner = runner._FimoScanner(
        fasta_file=promoters,
        motif_dir=motif_dir,
        output_dir=output,
        custom_bg={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
        pvalue_thresh=0.1,
        motif_filter={'rep_1'},
    )
    hits = scanner.run()

    assert not hits.empty
    assert set(hits.motif_id) == {'rep_1'}
    assert set(hits.sequence_name) == {'gene1'}
    assert hits.start.min() >= 1
    assert hits.stop.max() <= 8
    assert (output / RAW_FIMO_FILENAME).is_file()
