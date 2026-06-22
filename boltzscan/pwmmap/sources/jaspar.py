"""Download JASPAR CORE matrices via REST and convert PFMs to txt+meme."""
import json
from pathlib import Path

import requests

from boltzscan.pwmmap.models import RefTf
from boltzscan.pwmmap import pwmio

API = "https://jaspar.elixir.no/api/v1/matrix/"


def iter_jaspar_matrices(collection="CORE", page_size=500):
    url = f"{API}?collection={collection}&page_size={page_size}"
    while url:
        r = requests.get(url, headers={"Accept": "application/json"}, timeout=120)
        r.raise_for_status()
        page = r.json()
        for stub in page["results"]:
            d = requests.get(stub["url"], headers={"Accept": "application/json"},
                             timeout=120)
            d.raise_for_status()
            yield d.json()
        url = page.get("next")


def jaspar_refs_and_pwms(dest, txt_dir, meme_dir, refresh=False, **iter_kwargs):
    dest = Path(dest); dest.mkdir(parents=True, exist_ok=True)
    refs = []
    for det in iter_jaspar_matrices(**iter_kwargs):
        mid = det["matrix_id"]
        pfm = det.get("pfm")
        if not pfm:
            continue
        (dest / f"{mid}.json").write_text(json.dumps(det))
        pwmio.write_txt_and_meme(pfm, mid, txt_dir, meme_dir)
        species = "; ".join(s.get("name", "") for s in det.get("species", [])) or "?"
        fam = (det.get("family") or ["?"])[0]
        refs.append(RefTf(ref_id=f"jaspar:{mid}", source="jaspar", species=species,
                          family=fam, dbid=mid, motif_ids=[mid],
                          uniprot_ids=det.get("uniprot_ids") or []))
    return refs
