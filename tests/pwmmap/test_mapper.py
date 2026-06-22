import json
from boltzscan.pwmmap import mapper, align
from boltzscan.pwmmap.dbd import DbdRecord


def _make_store(tmp_path):
    root = tmp_path/"_refs"; (root/"motif_store"/"txt").mkdir(parents=True)
    (root/"motif_store"/"meme").mkdir(parents=True)
    (root/"ref_dbd.fasta").write_text(">cisbp:G1__PF00010__0\nMKRAHH\n")
    (root/"ref_index.tsv").write_text(
        "ref_id\tsource\tspecies\tfamily\tpfam_acc\tdbd_seq_id\tmotif_ids\n"
        "cisbp:G1\tcisbp\tAth\tbHLH\tPF00010\tcisbp:G1__PF00010__0\tM001_3.00\n")
    for ext, d in (("txt","txt"),("meme","meme")):
        (root/"motif_store"/d/f"M001_3.00.{ext}").write_text("x")
    return root


def test_map_species_transfers_motif_above_threshold(tmp_path, monkeypatch):
    root = _make_store(tmp_path)
    sp = tmp_path/"cm.fasta"; sp.write_text(">CM1\nMKRAHHWXYZ\n")
    # species DBD extraction -> one bHLH DBD on CM1
    monkeypatch.setattr(mapper.dbd, "extract_dbds",
        lambda *a, **k: [DbdRecord("CM1","PF00010","bHLH",1,6,"MKRAHH")])
    # blast -> CM1 matches ref G1 DBD at 0.95 (above bHLH cutoff 0.69)
    monkeypatch.setattr(mapper, "_blast",
        lambda *a, **k: [align.Hit("CM1","PF00010","cisbp:G1__PF00010__0","PF00010",0.95)])
    s = mapper.map_species(sp, tmp_path/"cm_pwms", refs_dir=root, threshold_mode="family")
    j = json.loads((tmp_path/"cm_pwms"/"tf2pwms.json").read_text())
    assert j["CM1"] == ["M001_3.00"]
    assert (tmp_path/"cm_pwms"/"txt"/"M001_3.00.txt").exists()
    assert s.n_mapped == 1


def test_below_threshold_not_transferred(tmp_path, monkeypatch):
    root = _make_store(tmp_path)
    sp = tmp_path/"cm.fasta"; sp.write_text(">CM1\nMKRAHHWXYZ\n")
    monkeypatch.setattr(mapper.dbd, "extract_dbds",
        lambda *a, **k: [DbdRecord("CM1","PF00010","bHLH",1,6,"MKRAHH")])
    monkeypatch.setattr(mapper, "_blast",
        lambda *a, **k: [align.Hit("CM1","PF00010","cisbp:G1__PF00010__0","PF00010",0.40)])
    s = mapper.map_species(sp, tmp_path/"cm_pwms", refs_dir=root, threshold_mode="family")
    assert s.n_mapped == 0
    assert json.loads((tmp_path/"cm_pwms"/"tf2pwms.json").read_text()) == {}
