from Bio import SeqIO

from boltzscan.utils.find_tf import find_transcription_factors


def _domtbl_line(protein_id):
    return (
        f"{protein_id} - 54 SRF-TF PF00319.24 48 1e-20 90.0 0.0 1 1 "
        "1e-20 1e-20 85.0 0.0 1 48 3 50 3 50 0.99 MADS-box\n"
    )


def test_reused_domtbl_is_filtered_to_input_fasta(tmp_path):
    proteins = tmp_path / "proteins.fasta"
    proteins.write_text(">wanted\nMKRIST\n")
    domtbl = tmp_path / "species.domtbl"
    domtbl.write_text(_domtbl_line("wanted") + _domtbl_line("other"))

    summary = find_transcription_factors(
        proteins,
        tmp_path / "out",
        domtbl=domtbl,
    )

    assert summary.n_total == 1
    assert summary.n_tf == 1
    assert summary.family_counts == {"M-type_MADS": 1}
    annotation_rows = summary.annotation_tsv.read_text().splitlines()
    assert len(annotation_rows) == 2
    assert annotation_rows[1].startswith("wanted\tM-type_MADS\t")
    assert [record.id for record in SeqIO.parse(summary.tf_fasta, "fasta")] == ["wanted"]
