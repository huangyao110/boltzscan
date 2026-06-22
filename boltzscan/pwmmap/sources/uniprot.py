"""Resolve reference TF protein sequences from UniProt.

JASPAR refs carry `uniprot_ids` -> fetch the accession's FASTA directly.
cisBP refs carry only a source `DBID` (heterogeneous namespaces: AGI locus,
Ensembl gene, FlyBase, SGD, ...) -> resolve via the UniProt *search* API
keyed by the DBID + organism name, which works across namespaces (idmapping
needs a different `from_db` per namespace and misses several, including
Arabidopsis AGI codes). Resolution is parallelized with a thread pool and
cached per ref_id so a re-run is cheap.
"""
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"
SEARCH = "https://rest.uniprot.org/uniprotkb/search"
IDMAP = "https://rest.uniprot.org/idmapping"


def _seq_from_fasta(text):
    if not text or not text.startswith(">"):
        return None
    seq = []
    for line in text.splitlines()[1:]:
        if line.startswith(">"):     # stop at a second record (defensive)
            break
        seq.append(line.strip())
    return "".join(seq)


def _get_with_retry(url, params=None, tries=3):
    for attempt in range(tries):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 200:
                return r.text
            if r.status_code in (429, 500, 502, 503):
                time.sleep(1 + attempt)
                continue
            return None
        except requests.RequestException:
            time.sleep(1 + attempt)
    return None


def fetch_fasta(acc):
    """Return raw FASTA text for a UniProt accession, or None on failure."""
    text = _get_with_retry(FASTA.format(acc))
    return text if (text and text.startswith(">")) else None


def search_fasta(dbid, species=None):
    """Resolve a source DBID to a UniProt sequence via the search API.

    Queries the DBID, optionally scoped to the organism, and returns the top
    hit's raw FASTA. Works across source namespaces (AGI, Ensembl, FlyBase, ...).
    """
    if not dbid:
        return None
    query = dbid
    if species:
        query = f'({dbid}) AND (organism_name:"{species.replace("_", " ")}")'
    text = _get_with_retry(SEARCH, params={"query": query, "format": "fasta", "size": 1})
    return text if (text and text.startswith(">")) else None


def map_ids_to_uniprot(ids, from_db="Ensembl_Genomes", to_db="UniProtKB"):
    """Batch map source IDs -> UniProt accessions (kept for compatibility; the
    search-based resolver below is preferred for heterogeneous cisBP DBIDs)."""
    if not ids:
        return {}
    sub = requests.post(f"{IDMAP}/run",
                        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
                        timeout=120)
    sub.raise_for_status()
    job = sub.json()["jobId"]
    for _ in range(60):
        st = requests.get(f"{IDMAP}/status/{job}", timeout=120).json()
        if st.get("jobStatus") == "FINISHED":
            break
        time.sleep(3)
    res = requests.get(f"{IDMAP}/results/{job}?size=500", timeout=120).json()
    return {x["from"]: x["to"]["primaryAccession"] if isinstance(x["to"], dict)
            else x["to"] for x in res.get("results", [])}


def _resolve_one(ref):
    if ref.uniprot_ids:
        return _seq_from_fasta(fetch_fasta(ref.uniprot_ids[0]))
    return _seq_from_fasta(search_fasta(ref.dbid, ref.species))


def resolve_sequences(refs, cache_path, max_workers=12):
    """Return {ref_id: protein_seq} for refs resolvable to a UniProt sequence.

    Cached per ref_id in cache_path (value "" marks a previously-unresolved ref
    so re-runs don't re-query it). Resolution runs concurrently.
    """
    cache_path = Path(cache_path)
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    out = {}
    todo = []
    for ref in refs:
        if ref.ref_id in cache:
            if cache[ref.ref_id]:
                out[ref.ref_id] = cache[ref.ref_id]
            continue
        todo.append(ref)

    if todo:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            results = list(ex.map(_resolve_one, todo))
        for ref, seq in zip(todo, results):
            cache[ref.ref_id] = seq or ""
            if seq:
                out[ref.ref_id] = seq
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache))

    return out
