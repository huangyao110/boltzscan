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


def download_cisbp(dest, refresh=False):
    dest = Path(dest)
    tf_info = next(dest.glob("**/TF_Information.txt"), None)
    pwms_dir = next((p for p in dest.glob("**/pwms_all_motifs") if p.is_dir()), None) \
        or next((p for p in dest.glob("**/pwms*") if p.is_dir()), None)
    if refresh or tf_info is None:
        _download_zip(TF_INFO_URL, dest / "tf_information")
        tf_info = next((dest / "tf_information").glob("**/TF_Information.txt"))
    if refresh or pwms_dir is None:
        _download_zip(PWMS_URL, dest / "pwms")
        pwms_dir = dest / "pwms"
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
