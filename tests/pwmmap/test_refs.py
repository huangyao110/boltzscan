from boltzscan.pwmmap import refs as R
from boltzscan.pwmmap.models import RefTf


def test_build_reference_db_offline(tmp_path, monkeypatch):
    # Two refs sharing a (real) Pfam DBD so hmmsearch finds them.
    refs = [RefTf("cisbp:G1","cisbp","Arabidopsis","bHLH","G1",["M001_3.00"]),
            RefTf("jaspar:MA1.1","jaspar","Arabidopsis","bHLH","MA1.1",["MA1.1"],["P9"])]
    # a real bHLH DBD sequence (~60 aa) padded; both refs get the same for the test
    bhlh = ("MKRAHHNALERRRRDHIKDSFSSLRDSVPSLQGEKASRAQILDKATEYIQYMRRKNHTHQQDIDDLKRQNALLEQQVRAL")
    monkeypatch.setattr(R, "_gather_cisbp", lambda d, refresh: ([refs[0]], None))
    monkeypatch.setattr(R, "_gather_jaspar", lambda d, t, m, refresh: [refs[1]])
    monkeypatch.setattr(R, "_resolve_seqs",
                        lambda rr, c: {"cisbp:G1": bhlh, "jaspar:MA1.1": bhlh})
    store = R.build_reference_db(tmp_path/"_refs", include_cisbp=True, include_jaspar=True)
    idx = R.load_ref_index(tmp_path/"_refs")
    assert store.ref_dbd_fasta.exists()
    # both refs produced at least one DBD record tagged with a Pfam acc
    assert any(v["pfam_acc"].startswith("PF") for v in idx.values())
