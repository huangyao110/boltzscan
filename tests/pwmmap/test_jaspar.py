from boltzscan.pwmmap.sources import jaspar
from boltzscan.pwmmap import pwmio


def test_detail_to_ref_and_pwm(tmp_path, monkeypatch):
    detail = {"matrix_id": "MA0570.1", "name": "ABF1", "family": ["bZIP"],
              "species": [{"name": "Arabidopsis thaliana"}], "uniprot_ids": ["P00001"],
              "pfm": {"A": [5,1], "C": [1,1], "G": [1,1], "T": [1,5]}}
    monkeypatch.setattr(jaspar, "iter_jaspar_matrices", lambda **k: iter([detail]))
    refs = jaspar.jaspar_refs_and_pwms(tmp_path/"raw", tmp_path/"txt", tmp_path/"meme")
    assert len(refs) == 1
    r = refs[0]
    assert r.source == "jaspar" and r.motif_ids == ["MA0570.1"]
    assert r.uniprot_ids == ["P00001"]
    assert (tmp_path/"txt"/"MA0570.1.txt").exists()
    assert (tmp_path/"meme"/"MA0570.1.meme").exists()
