from boltzscan.pwmmap.sources import cisbp

def _fake_tf_info(tmp_path):
    header = ("TF_ID\tFamily_ID\tTSource_ID\tMotif_ID\tMSource_ID\tDBID\tTF_Name\t"
              "TF_Species\tTF_Status\tFamily_Name\tDBDs\tDBD_Count\n")
    rows = [
        "T1\tF1\tS\tM001_3.00\tMS\tAT1G01250\tERF1\tArabidopsis_thaliana\tD\tAP2\tAP2\t1",
        "T2\tF1\tS\tM002_3.00\tMS\tAT1G01260\tERF2\tArabidopsis_thaliana\tD\tAP2\tAP2\t1",
        "T3\tF1\tS\t.\tMS\tAT9G99999\tX\tArabidopsis_thaliana\tD\tAP2\tAP2\t1",     # no motif
        "T4\tF1\tS\tM003_3.00\tMS\tSORGHUM01\tY\tSorghum_bicolor\tI\tAP2\tAP2\t1", # inferred
        "T1\tF1\tS\tM010_3.00\tMS\tAT1G01250\tERF1\tArabidopsis_thaliana\tD\tAP2\tAP2\t1", # 2nd motif for T1
    ]
    p = tmp_path / "TF_Information.txt"
    p.write_text(header + "\n".join(rows) + "\n")
    return p

def test_parse_keeps_determined_with_motif_and_groups(tmp_path):
    refs = {r.dbid: r for r in cisbp.parse_cisbp_refs(_fake_tf_info(tmp_path))}
    assert set(refs) == {"AT1G01250", "AT1G01260"}     # T3 (no motif) & T4 (inferred) dropped
    assert sorted(refs["AT1G01250"].motif_ids) == ["M001_3.00", "M010_3.00"]
    assert refs["AT1G01250"].source == "cisbp"
    assert refs["AT1G01260"].family == "AP2"
