"""Download JASPAR CORE matrices via REST and convert PFMs to txt+meme."""
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

from boltzscan.pwmmap.models import RefTf
from boltzscan.pwmmap import pwmio

API = "https://jaspar.elixir.no/api/v1/matrix/"


def _get_json(url, tries=3):
    """GET JSON with retry on transient errors; None on persistent failure."""
    for attempt in range(tries):
        try:
            r = requests.get(url, headers={"Accept": "application/json"}, timeout=120)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503):
                time.sleep(1 + attempt)
                continue
            return None
        except requests.RequestException:
            time.sleep(1 + attempt)
    return None


def _list_stubs(collection="CORE", page_size=500):
    """Paginate the matrix list (fast) -> list of {matrix_id, url} stubs."""
    url = f"{API}?collection={collection}&page_size={page_size}"
    stubs = []
    while url:
        page = _get_json(url)
        if not page:
            break
        stubs.extend(page.get("results", []))
        url = page.get("next")
    return stubs


def iter_jaspar_matrices(collection="CORE", page_size=500, max_workers=16):
    """Yield JASPAR matrix detail dicts, fetching the per-matrix details
    concurrently (the list is small; the details are the slow part). A detail
    that fails after retries is skipped rather than crashing the whole build."""
    stubs = _list_stubs(collection, page_size)
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        for det in ex.map(lambda s: _get_json(s["url"]), stubs):
            if det is not None:
                yield det


def jaspar_refs_and_pwms(dest, txt_dir, meme_dir, refresh=False, **iter_kwargs):
    dest = Path(dest); dest.mkdir(parents=True, exist_ok=True)
    refs = []
    for det in iter_jaspar_matrices(**iter_kwargs):
        mid = det.get("matrix_id")
        pfm = det.get("pfm")
        if not mid or not pfm:
            continue
        (dest / f"{mid}.json").write_text(json.dumps(det))
        try:
            pwmio.write_txt_and_meme(pfm, mid, txt_dir, meme_dir)
        except Exception:                                # malformed pfm -> skip
            continue
        species = "; ".join(s.get("name", "") for s in det.get("species", [])) or "?"
        fam = (det.get("family") or ["?"])[0]
        refs.append(RefTf(ref_id=f"jaspar:{mid}", source="jaspar", species=species,
                          family=fam, dbid=mid, motif_ids=[mid],
                          uniprot_ids=det.get("uniprot_ids") or []))
    return refs
