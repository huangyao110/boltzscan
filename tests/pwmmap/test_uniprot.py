from boltzscan.pwmmap.sources import uniprot


def test_parse_fasta_body():
    fa = ">sp|P00001|X\nMKTA\nYYYY\n"
    assert uniprot._seq_from_fasta(fa) == "MKTAYYYY"


def test_resolve_uses_uniprot_id_directly(monkeypatch, tmp_path):
    from boltzscan.pwmmap.models import RefTf
    monkeypatch.setattr(uniprot, "fetch_fasta", lambda acc: ">x\nMKK\n" if acc == "P9" else None)
    refs = [RefTf("jaspar:MA1", "jaspar", "sp", "bZIP", "MA1", ["MA1"], ["P9"])]
    seqs = uniprot.resolve_sequences(refs, cache_path=tmp_path/"cache.json")
    assert seqs["jaspar:MA1"] == "MKK"
