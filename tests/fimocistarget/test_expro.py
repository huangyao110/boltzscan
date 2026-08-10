from types import SimpleNamespace

from boltzscan.fimocistarget.expro import extract_promoter_regions


def _inputs(tmp_path, output_format):
    genome = tmp_path / 'genome.fa'
    genome.write_text('>chr1\nAACCGGTTAACCGGTT\n')
    gff = tmp_path / 'genes.gff3'
    gff.write_text(
        '##gff-version 3\n'
        'chr1\ttest\tgene\t5\t12\t.\t+\t.\tID=gene1\n'
        'chr1\ttest\tCDS\t5\t8\t.\t+\t0\tParent=gene1\n'
        '###\n'
    )
    return SimpleNamespace(
        genome=str(genome),
        gff=str(gff),
        output=str(tmp_path / 'promoters'),
        format=output_format,
        upstream=2,
        downstream=2,
    )


def test_promoter_bed_format_writes_only_promoter_bed(tmp_path):
    extract_promoter_regions(_inputs(tmp_path, 'bed'))

    assert (tmp_path / 'promoters.bed').read_text() == (
        'chr1\t2\t6\tgene1_promoter\t1000\t+\n'
    )
    assert not (tmp_path / 'promoters.fasta').exists()


def test_promoter_fasta_format_writes_only_promoter_fasta(tmp_path):
    args = _inputs(tmp_path, 'fasta')
    args.output = str(tmp_path / 'promoters.fasta')
    extract_promoter_regions(args)

    fasta = (tmp_path / 'promoters.fasta').read_text()
    assert fasta.startswith('>gene1 promoter_cds_anchor chr1:3-6(+) anchor:5\n')
    assert fasta.endswith('CCgg\n')
    assert not (tmp_path / 'promoters.bed').exists()
