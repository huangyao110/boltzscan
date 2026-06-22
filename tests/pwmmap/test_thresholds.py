from boltzscan.pwmmap import thresholds as T

def test_pfam_maps_to_family():
    assert T.family_for_pfam("PF00010") == "bHLH"
    assert T.family_for_pfam("PF00046") == "Homeobox"
    assert T.family_for_pfam("PF99999") is None

def test_family_cutoff_used_in_family_mode():
    # bHLH has a calibrated cutoff distinct from the default
    assert T.cutoff_for("PF00010", mode="family") == T.FAMILY_CUTOFF["bHLH"]

def test_unknown_family_falls_back_to_default():
    assert T.cutoff_for("PF99999", mode="family") == T.DEFAULT_CUTOFF

def test_global_mode_ignores_family():
    assert T.cutoff_for("PF00010", mode="global", global_thr=0.55) == 0.55
