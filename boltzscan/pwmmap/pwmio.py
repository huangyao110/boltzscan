"""Convert PWMs between JASPAR PFM (counts), cisBP-style txt (probabilities),
and MEME, reusing the existing txt->MEME converter."""
import shutil
from pathlib import Path

from boltzscan.utils.io_utils import txt_to_meme


def pfm_to_txt(pfm, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["A", "C", "G", "T"]
    width = len(pfm["A"])
    with open(path, "w") as fh:
        fh.write("Pos\tA\tC\tG\tT\n")
        for i in range(width):
            counts = [float(pfm[c][i]) for c in cols]
            tot = sum(counts) or 1.0
            probs = [x / tot for x in counts]
            fh.write(f"{i+1}\t" + "\t".join(repr(round(p, 6)) for p in probs) + "\n")
    return path


def write_txt_and_meme(pfm, motif_id, txt_dir, meme_dir):
    txt_dir, meme_dir = Path(txt_dir), Path(meme_dir)
    txt = pfm_to_txt(pfm, txt_dir / f"{motif_id}.txt")
    meme_dir.mkdir(parents=True, exist_ok=True)
    txt_to_meme(input_path=str(txt), output_path=str(meme_dir / f"{motif_id}.meme"),
                force=True)
    return txt, meme_dir / f"{motif_id}.meme"


def copy_cisbp_pwm(src_txt, motif_id, txt_dir, meme_dir):
    # Read source before doing anything so we can detect degenerate files.
    src_txt = Path(src_txt)
    raw = src_txt.read_text(errors="replace") if src_txt.exists() else ""
    lines = [l for l in raw.splitlines() if l.strip()]

    # Detect degenerate: empty, MEME-format, or no row starting with digit or "Pos".
    if not lines:
        print(f"[pwmio] skip degenerate cisBP PWM {motif_id}")
        return (None, None)
    if lines[0].startswith("MEME version"):
        print(f"[pwmio] skip degenerate cisBP PWM {motif_id}")
        return (None, None)
    has_data = any(l[0].isdigit() or l.startswith("Pos") for l in lines)
    if not has_data:
        print(f"[pwmio] skip degenerate cisBP PWM {motif_id}")
        return (None, None)

    txt_dir, meme_dir = Path(txt_dir), Path(meme_dir)
    txt_dir.mkdir(parents=True, exist_ok=True)
    meme_dir.mkdir(parents=True, exist_ok=True)
    dst_txt = txt_dir / f"{motif_id}.txt"
    shutil.copyfile(src_txt, dst_txt)
    meme_path = meme_dir / f"{motif_id}.meme"
    txt_to_meme(input_path=str(dst_txt), output_path=str(meme_path), force=True)
    return dst_txt, meme_path
