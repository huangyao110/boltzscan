from boltzscan.pwmmap import dbd

def test_parse_domtbl_keeps_only_dbd_pfams(tiny_domtbl):
    rows = dbd.parse_domtbl_dbds(tiny_domtbl)
    assert rows == [("geneA", "PF00010", 10, 65)]  # PF99999 dropped

def test_extract_dbds_slices_envelope(tiny_domtbl, tiny_proteins):
    recs = dbd.extract_dbds(tiny_proteins, domtbl=tiny_domtbl)
    assert len(recs) == 1
    r = recs[0]
    assert r.tf_id == "geneA" and r.pfam_acc == "PF00010" and r.family == "bHLH"
    assert len(r.seq) == 65 - 10 + 1   # 1-based inclusive envelope

def test_write_dbd_fasta_id_format(tiny_domtbl, tiny_proteins, tmp_path):
    recs = dbd.extract_dbds(tiny_proteins, domtbl=tiny_domtbl)
    out = tmp_path / "dbd.fasta"
    dbd.write_dbd_fasta(recs, out)
    text = out.read_text()
    assert text.startswith(">geneA__PF00010__0\n")
