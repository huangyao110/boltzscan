"""ESMFold2 protein-DNA folding worker — runs INSIDE the esmfold2 conda env.

Invoked out-of-process by boltzscan.predict.runners.run_esmfold (mirrors the
esm_embed ESMC worker pattern). Reads a directory of Boltz-style input YAMLs
(protein chain + msa + dna strands), folds each complex with ESMFold2 + ESMC,
and writes a Boltz-compatible output tree so the same downstream ranking/ipSAE
works for either engine:

  <out_dir>/predictions/<name>/<name>_model_0.cif
  <out_dir>/predictions/<name>/confidence_<name>_model_0.json   # {iptm, ptm, ...}

Usage (esmfold2 env python):
  python _esmfold_worker.py --input-dir <yaml_dir> --out-dir <out> \
      [--num-loops 3 --sampling-steps 50 --diffusion-samples 1 --seed 42]
"""
import argparse
import json
import os
import sys
import traceback
from pathlib import Path

import yaml

DEFAULT_ESMFOLD2 = os.environ.get(
    "BOLTZSCAN_ESMFOLD2_WEIGHTS", "/media/zlab/新加卷/esm_weights/ESMFold2")
DEFAULT_ESMC = os.environ.get(
    "BOLTZSCAN_ESMC_WEIGHTS", "/media/zlab/新加卷/esm_weights/ESMC-6B")


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--weights", default=DEFAULT_ESMFOLD2)
    ap.add_argument("--esmc-weights", default=DEFAULT_ESMC)
    ap.add_argument("--num-loops", type=int, default=3)
    ap.add_argument("--sampling-steps", type=int, default=50)
    ap.add_argument("--diffusion-samples", type=int, default=1)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

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

    pred_root = Path(args.out_dir) / "predictions"
    pred_root.mkdir(parents=True, exist_ok=True)
    n_ok = n_fail = 0
    for y in yamls:
        name = y.stem
        outd = pred_root / name
        if (outd / f"{name}_model_0.cif").exists():
            n_ok += 1
            continue
        try:
            inp = _build_input(y)
            res = builder.fold(model, inp, num_loops=args.num_loops,
                               num_sampling_steps=args.sampling_steps,
                               num_diffusion_samples=args.diffusion_samples,
                               seed=args.seed, complex_id=name)
            if isinstance(res, list):
                res = res[0]
            outd.mkdir(parents=True, exist_ok=True)
            (outd / f"{name}_model_0.cif").write_text(res.complex.to_mmcif())
            conf = {
                "iptm": float(res.iptm) if res.iptm is not None else None,
                "ptm": float(res.ptm) if res.ptm is not None else None,
                "engine": "esmfold2",
            }
            (outd / f"confidence_{name}_model_0.json").write_text(json.dumps(conf))
            n_ok += 1
            print(f"[esmfold] {name} iptm={conf['iptm']}", file=sys.stderr)
        except Exception:
            n_fail += 1
            print(f"[esmfold] FAILED {name}:\n{traceback.format_exc()}", file=sys.stderr)
    print(f"[esmfold] done: {n_ok} ok, {n_fail} failed -> {pred_root}", file=sys.stderr)


if __name__ == "__main__":
    main()
