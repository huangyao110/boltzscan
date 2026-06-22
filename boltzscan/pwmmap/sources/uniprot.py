"""Resolve reference TF protein sequences from UniProt."""
import json
import time
from pathlib import Path

import requests

FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"
IDMAP = "https://rest.uniprot.org/idmapping"


def _seq_from_fasta(text):
    if not text or not text.startswith(">"):
        return None
    return "".join(text.splitlines()[1:])


def fetch_fasta(acc):
    """Return raw FASTA text for a UniProt accession, or None on failure."""
    try:
        r = requests.get(FASTA.format(acc), timeout=60)
        if r.status_code == 200 and r.text.startswith(">"):
            return r.text
    except requests.RequestException:
        return None
    return None


def map_ids_to_uniprot(ids, from_db="Ensembl_Genomes", to_db="UniProtKB"):
    """Batch map source IDs -> UniProt accessions; returns {src_id: uniprot_acc}."""
    if not ids:
        return {}
    sub = requests.post(f"{IDMAP}/run",
                        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
                        timeout=120)
    sub.raise_for_status()
    job = sub.json()["jobId"]
    for _ in range(60):
        st = requests.get(f"{IDMAP}/status/{job}", timeout=120).json()
        if st.get("results") is not None or st.get("jobStatus") == "FINISHED":
            break
        time.sleep(3)
    res = requests.get(f"{IDMAP}/results/{job}?size=500", timeout=120).json()
    return {x["from"]: x["to"]["primaryAccession"] if isinstance(x["to"], dict)
            else x["to"] for x in res.get("results", [])}


def resolve_sequences(refs, cache_path):
    cache_path = Path(cache_path)
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    out = {}
    pending_cisbp = []
    for ref in refs:
        if ref.ref_id in cache:
            if cache[ref.ref_id]:
                out[ref.ref_id] = cache[ref.ref_id]
            continue
        acc = ref.uniprot_ids[0] if ref.uniprot_ids else None
        if acc:
            seq = _seq_from_fasta(fetch_fasta(acc))
            cache[ref.ref_id] = seq or ""
            if seq:
                out[ref.ref_id] = seq
        else:
            pending_cisbp.append(ref)
    # cisBP refs without a direct UniProt id: batch-map their DBID
    mapped = map_ids_to_uniprot([r.dbid for r in pending_cisbp]) if pending_cisbp else {}
    for r in pending_cisbp:
        acc = mapped.get(r.dbid)
        seq = _seq_from_fasta(fetch_fasta(acc)) if acc else None
        cache[r.ref_id] = seq or ""
        if seq:
            out[r.ref_id] = seq
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache))
    return out
