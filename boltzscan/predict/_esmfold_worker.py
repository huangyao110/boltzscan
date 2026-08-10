"""ESMFold2 protein-DNA folding worker — runs INSIDE the esmfold2 conda env.

Invoked out-of-process by boltzscan.predict.runners.run_esmfold. Reads a
directory of Boltz-style input YAMLs
(protein chain + msa + dna strands), folds each complex with ESMFold2 + ESMC,
and writes a Boltz-compatible output tree so the same downstream ranking/ipSAE
works for either engine:

  <out_dir>/predictions/<name>/<name>_model_0.cif
  <out_dir>/predictions/<name>/pae_<name>_model_0.npz
  <out_dir>/predictions/<name>/plddt_<name>_model_0.npz
  <out_dir>/predictions/<name>/confidence_<name>_model_0.json

Usage (esmfold2 env python):
  python _esmfold_worker.py --input-dir <yaml_dir> --out-dir <out> \
      [--seed N]
"""
import argparse
import json
import os
import sys
import traceback
from pathlib import Path

import yaml

DEFAULT_ESMFOLD2 = os.environ.get(
    "BOLTZSCAN_ESMFOLD2_WEIGHTS", "biohub/ESMFold2")
DEFAULT_ESMC = os.environ.get(
    "BOLTZSCAN_ESMC_WEIGHTS", "biohub/ESMC-6B")


def _build_input(yaml_path):
    from esm.utils.msa.msa import MSA
    from esm.utils.structure.input_builder import (
        ProteinInput, DNAInput, RNAInput, StructurePredictionInput,
    )
    spec = yaml.safe_load(Path(yaml_path).read_text())
    seqs = []
    for entry in spec.get("sequences", []):
        if "protein" in entry:
            p = entry["protein"]
            msa = None
            mp = p.get("msa")
            if mp and mp != "empty" and Path(mp).exists():
                msa = MSA.from_a3m(mp)
            seqs.append(ProteinInput(id=p["id"], sequence=p["sequence"], msa=msa))
        elif "dna" in entry:
            d = entry["dna"]
            seqs.append(DNAInput(id=d["id"], sequence=d["sequence"]))
        elif "rna" in entry:
            r = entry["rna"]
            seqs.append(RNAInput(id=r["id"], sequence=r["sequence"]))
    return StructurePredictionInput(sequences=seqs)


def _artifact_paths(outd, name):
    stem = f"{name}_model_0"
    return (
        Path(outd) / f"{stem}.cif",
        Path(outd) / f"pae_{stem}.npz",
        Path(outd) / f"plddt_{stem}.npz",
        Path(outd) / f"confidence_{stem}.json",
    )


def _as_numpy(value):
    if hasattr(value, "detach"):
        value = value.detach()
    if hasattr(value, "float"):
        value = value.float()
    if hasattr(value, "cpu"):
        value = value.cpu()
    return value.numpy() if hasattr(value, "numpy") else value


def _as_float(value):
    if value is None:
        return None
    if hasattr(value, "detach"):
        value = value.detach()
    if hasattr(value, "cpu"):
        value = value.cpu()
    return float(value)


def _write_artifacts(res, name, outd):
    """Write the minimal Boltz/ipSAE artifact set for one ESMFold2 result."""
    import io
    import numpy as np
    import biotite.structure.io.pdbx as pdbx

    cif_path, pae_path, plddt_path, confidence_path = _artifact_paths(outd, name)
    cif_in = pdbx.CIFFile.read(io.StringIO(res.complex.to_mmcif()))
    structure = pdbx.get_structure(cif_in, model=1)
    structure.set_annotation("occupancy", np.ones(structure.array_length()))
    if "b_factor" not in structure.get_annotation_categories():
        structure.set_annotation("b_factor", np.zeros(structure.array_length()))
    cif_out = pdbx.CIFFile()
    pdbx.set_structure(cif_out, structure)
    cif_out.write(str(cif_path))

    pae = getattr(res, "pae", None)
    plddt = getattr(res, "plddt", None)
    if pae is None or plddt is None:
        raise ValueError("ESMFold2 result is missing PAE or pLDDT")
    np.savez_compressed(pae_path, pae=_as_numpy(pae))
    plddt = np.asarray(_as_numpy(plddt))
    if plddt.size and plddt.max() > 1.5:
        plddt = plddt / 100.0
    np.savez_compressed(plddt_path, plddt=plddt)

    conf = {
        "iptm": _as_float(getattr(res, "iptm", None)),
        "ptm": _as_float(getattr(res, "ptm", None)),
        "engine": "esmfold2",
    }
    pair_iptm = getattr(res, "pair_chains_iptm", None)
    if pair_iptm is not None:
        matrix = np.asarray(_as_numpy(pair_iptm))
        conf["pair_chains_iptm"] = {
            str(i): {str(j): float(matrix[i, j]) for j in range(matrix.shape[1])}
            for i in range(matrix.shape[0])
        }
        lookup = getattr(getattr(res.complex, "metadata", None), "chain_lookup", None)
        if lookup:
            conf["chain_order"] = [lookup[i] for i in sorted(lookup)]
    confidence_path.write_text(json.dumps(conf))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--predictions-dir", default=None,
                    help="Optional direct directory for per-model prediction folders")
    ap.add_argument("--weights", default=DEFAULT_ESMFOLD2)
    ap.add_argument("--esmc-weights", default=DEFAULT_ESMC)
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    # Set this before importing torch or creating a CUDA context.  The vendored
    # ESMFold2 seed context also enables deterministic algorithms around input
    # preparation and folding; cuBLAS needs its workspace policy earlier.
    if args.seed is not None:
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

    # A CCD pickle produced by a newer RDKit patch release can emit one
    # compatibility warning per molecule while it is loaded.  The molecules
    # remain readable, but thousands of duplicate lines obscure the concise
    # prediction summary.  Suppressing RDKit warnings does not alter inference.
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")
    except ImportError:
        pass

    import torch
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    from esm.models.esmfold2.processor import ESMFold2InputBuilder

    yamls = sorted(Path(args.input_dir).glob("*.yaml"))
    if not yamls:
        sys.exit(f"no .yaml inputs under {args.input_dir}")
    print(f"[esmfold] {len(yamls)} inputs; loading ESMFold2 ...", file=sys.stderr)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = ESMFold2Model.from_pretrained(args.weights, load_esmc=False).to(device).eval()
    model.load_esmc(args.esmc_weights)
    builder = ESMFold2InputBuilder()

    pred_root = Path(args.predictions_dir) if args.predictions_dir else Path(args.out_dir) / "predictions"
    pred_root.mkdir(parents=True, exist_ok=True)
    n_ok = n_fail = 0
    for y in yamls:
        name = y.stem
        outd = pred_root / name
        if all(path.exists() for path in _artifact_paths(outd, name)):
            n_ok += 1
            continue
        try:
            inp = _build_input(y)
            fold_kwargs = {"complex_id": name}
            if args.seed is not None:
                fold_kwargs["seed"] = args.seed
            res = builder.fold(model, inp, **fold_kwargs)
            if isinstance(res, list):
                res = res[0]
            outd.mkdir(parents=True, exist_ok=True)
            _write_artifacts(res, name, outd)
            n_ok += 1
            print(f"[esmfold] {name} iptm={_as_float(getattr(res, 'iptm', None))}", file=sys.stderr)
        except Exception:
            n_fail += 1
            print(f"[esmfold] FAILED {name}:\n{traceback.format_exc()}", file=sys.stderr)
    print(f"[esmfold] done: {n_ok} ok, {n_fail} failed -> {pred_root}", file=sys.stderr)
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
