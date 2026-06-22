from boltzscan.fimocistarget import cpg


def test_no_island_in_at_rich():
    assert cpg.find_cpg_islands("AT" * 300, window=200) == []


def test_detects_cpg_island():
    seq = "AT" * 200 + "CG" * 200 + "AT" * 200      # CG core at 400..800
    isl = cpg.find_cpg_islands(seq, window=200, min_gc=0.5, min_oe=0.6)
    assert len(isl) >= 1
    big = max(isl, key=lambda i: i["length"])
    # merged island covers the CG core; boundary half-windows can extend it
    # ~100bp into AT flanks, diluting GC but staying >= min_gc.
    assert big["gc"] >= 0.5 and big["oe"] >= 1.0 and big["n_cpg"] > 100
    assert big["length"] >= 200 and big["start"] <= 400 and big["end"] >= 800


def test_short_seq_returns_empty():
    assert cpg.find_cpg_islands("CG" * 50, window=200) == []
