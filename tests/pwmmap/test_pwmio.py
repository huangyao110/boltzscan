from boltzscan.pwmmap import pwmio


def test_pfm_to_txt_normalizes_counts(tmp_path):
    pfm = {"A": [10, 0], "C": [0, 0], "G": [0, 10], "T": [0, 0]}
    p = pwmio.pfm_to_txt(pfm, tmp_path / "MA0001.1.txt")
    lines = p.read_text().strip().splitlines()
    assert lines[0].split() == ["Pos", "A", "C", "G", "T"]
    # position 1 is all-A -> A prob 1.0; position 2 all-G -> G prob 1.0
    assert lines[1].split()[1] == "1.0"
    assert lines[2].split()[3] == "1.0"


def test_write_txt_and_meme_emits_both(tmp_path):
    pfm = {"A": [5, 1], "C": [1, 1], "G": [1, 1], "T": [1, 5]}
    txt, meme = pwmio.write_txt_and_meme(pfm, "MA0001.1", tmp_path/"txt", tmp_path/"meme")
    assert txt.exists() and meme.exists()
    assert "MOTIF MA0001.1" in meme.read_text()


def test_copy_cisbp_skips_meme_format(tmp_path):
    src = tmp_path / "rose_x_pwm_src.txt"
    src.write_text("MEME version 5.5.8\n\nALPHABET= ACGU\n")
    result = pwmio.copy_cisbp_pwm(src, "rose_x_pwm", tmp_path / "txt", tmp_path / "meme")
    assert result == (None, None)
    assert not (tmp_path / "meme" / "rose_x_pwm.meme").exists()


def test_copy_cisbp_normal(tmp_path):
    src = tmp_path / "valid_pwm_src.txt"
    src.write_text("Pos\tA\tC\tG\tT\n1\t1.0\t0.0\t0.0\t0.0\n")
    dst_txt, meme_path = pwmio.copy_cisbp_pwm(src, "valid_pwm", tmp_path / "txt", tmp_path / "meme")
    assert dst_txt is not None and dst_txt.exists()
    assert meme_path is not None and meme_path.exists()
    assert "MOTIF" in meme_path.read_text()
