#!/usr/bin/env python
"""Build an auditable general protein--DNA Pfam-domain whitelist.

The primary evidence is the official Pfam2GO mapping to GO:0003677 (DNA
binding) or one of its ``is_a`` descendants.  The production plant-TF list is
then retained, and a small reviewed extension covers clear DBDs whose Pfam2GO
records currently lack a DNA-binding term.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re

import pandas as pd
import requests

from boltzscan.utils.find_tf import DBD_ACCS


GO_ROOT = "GO:0003677"
GO_OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
PFAM2GO_URL = "https://current.geneontology.org/ontology/external2go/pfam2go"

# Reviewed additions are deliberately exact accessions rather than a broad
# name-pattern rule.  Each one names a canonical DNA-binding fold/domain.
CURATED_EXTENSIONS = {
    "PF00505": "HMG-box DNA-binding domain",
    "PF09011": "second HMG-box DNA-binding domain",
    "PF03299": "transcription factor AP-2 DNA-binding family",
    "PF10497": "4CXXC zinc-finger DNA-binding domain",
    "PF10491": "Nrf1 DNA-binding family",
    "PF04545": "sigma-70 region 4 DNA-binding domain",
    "PF04814": "HNF-1 N-terminal DNA-binding family",
    "PF02467": "WhiB transcriptional-regulator family",
    "PF12844": "helix-turn-helix domain 19",
    "PF13560": "helix-turn-helix domain 31",
    "PF01381": "helix-turn-helix domain 3",
    "PF13443": "helix-turn-helix domain 26",
    "PF13744": "helix-turn-helix domain 37",
    "PF13412": "helix-turn-helix domain 24",
    "PF15731": "MqsA antitoxin DNA-binding domain",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url, path):
    path = Path(path)
    if path.is_file():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    staged = path.with_suffix(path.suffix + ".part")
    with requests.get(url, stream=True, timeout=(20, 180)) as response:
        response.raise_for_status()
        with staged.open("wb") as handle:
            for block in response.iter_content(1024 * 1024):
                if block:
                    handle.write(block)
    staged.replace(path)
    return path


def go_terms(obo_path):
    parents = defaultdict(set)
    names = {}
    current = None
    with Path(obo_path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                current = None
            elif line.startswith("id: GO:"):
                current = line.split(": ", 1)[1]
            elif current and line.startswith("name: "):
                names[current] = line[6:]
            elif current and line.startswith("is_a: GO:"):
                parents[current].add(line.split()[1])

    def is_dna_binding(term):
        pending, seen = [term], set()
        while pending:
            node = pending.pop()
            if node == GO_ROOT:
                return True
            if node in seen:
                continue
            seen.add(node)
            pending.extend(parents.get(node, ()))
        return False

    return names, {term for term in names if is_dna_binding(term)}


def pfam_mapping(path, dna_binding_terms):
    pattern = re.compile(
        r"^Pfam:(PF\d+)\s+(\S+)\s+>.*;\s+(GO:\d+)\s*$"
    )
    names, terms = {}, defaultdict(set)
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            match = pattern.match(raw.strip())
            if not match:
                continue
            accession, name, go_term = match.groups()
            names[accession] = name
            if go_term in dna_binding_terms:
                terms[accession].add(go_term)
    return names, terms


def observed_counts(domtbl):
    counts = defaultdict(int)
    names = {}
    if domtbl is None:
        return counts, names
    with Path(domtbl).open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            accession = fields[4].split(".")[0]
            counts[accession] += 1
            names[accession] = fields[3]
    return counts, names


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--observed-domtbl")
    args = parser.parse_args()

    output = Path(args.output)
    cache = Path(args.cache_dir)
    obo = download(GO_OBO_URL, cache / "go-basic.obo")
    pfam2go = download(PFAM2GO_URL, cache / "pfam2go")
    go_names, dna_terms = go_terms(obo)
    pfam_names, mapped_terms = pfam_mapping(pfam2go, dna_terms)
    counts, observed_names = observed_counts(args.observed_domtbl)

    accessions = set(mapped_terms) | set(DBD_ACCS) | set(CURATED_EXTENSIONS)
    rows = []
    for accession in sorted(accessions):
        sources = []
        if accession in mapped_terms:
            sources.append("GO:0003677-descendant")
        if accession in DBD_ACCS:
            sources.append("BoltzScan-plant-DBD")
        if accession in CURATED_EXTENSIONS:
            sources.append("curated-extension")
        terms = sorted(mapped_terms.get(accession, ()))
        rows.append({
            "pfam_acc": accession,
            "pfam_name": pfam_names.get(accession, observed_names.get(accession, "")),
            "evidence": ";".join(sources),
            "go_terms": ";".join(terms),
            "go_names": ";".join(go_names[term] for term in terms),
            "curation_note": CURATED_EXTENSIONS.get(accession, ""),
            "observed_domain_hits": counts.get(accession, 0),
        })
    frame = pd.DataFrame(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, sep="\t", index=False)
    metadata = {
        "generated": datetime.now().astimezone().isoformat(timespec="seconds"),
        "go_root": GO_ROOT,
        "go_root_name": go_names[GO_ROOT],
        "go_descendants_including_root": len(dna_terms),
        "pfam_accessions": len(frame),
        "sources": {
            "go_obo": {"url": GO_OBO_URL, "sha256": sha256(obo)},
            "pfam2go": {"url": PFAM2GO_URL, "sha256": sha256(pfam2go)},
        },
        "selection": "GO is_a descendants + production plant DBDs + reviewed exact extensions",
    }
    output.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote {len(frame)} Pfam accessions -> {output}")


if __name__ == "__main__":
    main()
