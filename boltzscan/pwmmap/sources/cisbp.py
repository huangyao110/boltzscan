"""Download cisBP entire dataset (TF_Information + PWMs) and parse high-quality refs."""
import os
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

from boltzscan.pwmmap.models import RefTf

BASE = "https://cisbp.ccbr.utoronto.ca/data/3_10/DataFiles/Bulk_downloads/EntireDataset"
TF_INFO_URL = f"{BASE}/TF_Information.zip"
PWMS_URL = f"{BASE}/PWMs.zip"


def _download_zip(url, dest_dir):
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp:
            for chunk in r.iter_content(chunk_size=1 << 20):
                tmp.write(chunk)
            tmp_path = tmp.name
    try:
        with zipfile.ZipFile(tmp_path) as z:
            z.extractall(dest_dir)
    finally:
        os.unlink(tmp_path)
    return dest_dir


def _extract_local_zip(zip_path):
    """Extract a local .zip in place (cisBP wraps its .txt files in inner zips)."""
    with zipfile.ZipFile(zip_path) as z:
        z.extractall(Path(zip_path).parent)


def download_cisbp(dest, refresh=False):
    """Return (TF_Information.txt path, dir containing M*.txt PWMs).

    The cisBP entire-dataset archives are nested: TF_Information.zip contains
    inner *.txt.zip files, and PWMs.zip unpacks to a pwms/ subdir of M*.txt.
    """
    dest = Path(dest)
    tf_dir = dest / "tf_information"
    pwms_root = dest / "pwms"

    tf_info = next(dest.glob("**/TF_Information.txt"), None)
    if refresh or tf_info is None:
        _download_zip(TF_INFO_URL, tf_dir)                       # outer -> inner *.txt.zip
        inner = next(tf_dir.glob("**/TF_Information.txt.zip"), None)
        if inner is not None:
            _extract_local_zip(inner)                            # inner -> TF_Information.txt
        tf_info = next(tf_dir.glob("**/TF_Information.txt"))

    pwm_txt = next(pwms_root.glob("**/M*.txt"), None) if pwms_root.exists() else None
    if refresh or pwm_txt is None:
        _download_zip(PWMS_URL, pwms_root)
        for z in list(pwms_root.glob("**/*.zip")):              # tolerate nested zips
            _extract_local_zip(z)
        pwm_txt = next(pwms_root.glob("**/M*.txt"))
    pwms_dir = pwm_txt.parent                                    # the actual dir of M*.txt

    return tf_info, pwms_dir


def parse_cisbp_refs(tf_info_path):
    df = pd.read_csv(tf_info_path, sep="\t", dtype=str).fillna("")
    df = df[(df["TF_Status"] == "D") & (df["Motif_ID"] != "") & (df["Motif_ID"] != ".")]
    by_dbid = defaultdict(lambda: {"motifs": set(), "family": "", "species": ""})
    for _, row in df.iterrows():
        d = by_dbid[row["DBID"]]
        d["motifs"].add(row["Motif_ID"])
        d["family"] = row["Family_Name"]
        d["species"] = row["TF_Species"]
    return [
        # uniprot_ids intentionally empty: cisBP TF_Information.txt has no UniProt column;
        # sequences are resolved downstream via the UniProt ID-mapping step (Task 5).
        RefTf(ref_id=f"cisbp:{dbid}", source="cisbp", species=v["species"],
              family=v["family"], dbid=dbid, motif_ids=sorted(v["motifs"]))
        for dbid, v in by_dbid.items()
    ]
