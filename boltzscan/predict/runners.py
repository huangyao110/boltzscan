"""Run protein-DNA structure prediction with a choice of engine.

Both engines live in their own conda env and are driven out-of-process from the
same Boltz-style YAML input dir:

- boltz   -> the vendored Boltz (`boltz_ode`/`boltz2`) in the `boltz` env. Fixes
             the old `bscan.py boltzscan` bug of calling boltz with sys.executable
             (the boltzscan env has no torch); we use the boltz env python +
             PYTHONPATH=vendored src so `import boltz` finds the boltz_ode source.
- esmfold -> ESMFold2 (+ESMC) in the `esmfold2` env via _esmfold_worker.py.

Both write a Boltz-compatible output tree (predictions/<name>/...cif +
confidence_*.json) so the same downstream ranking/ipSAE applies either way.
Env paths are overridable via the BOLTZSCAN_* environment variables.
"""
import os
import subprocess
import sys
from pathlib import Path

_PKG = Path(__file__).resolve().parents[1]          # boltzscan/
_BOLTZ_SRC = _PKG / "boltz" / "src"
_BOLTZ_MAIN = _BOLTZ_SRC / "boltz" / "main.py"
_ESMFOLD_WORKER = Path(__file__).resolve().parent / "_esmfold_worker.py"

DEFAULT_BOLTZ_PY = os.environ.get(
    "BOLTZSCAN_BOLTZ_PYTHON", "/home/zlab/miniconda3/envs/boltz/bin/python")
DEFAULT_ESMFOLD_PY = os.environ.get(
    "BOLTZSCAN_ESMFOLD_PYTHON", "/home/zlab/miniconda3/envs/esmfold2/bin/python")


def run_boltz(input_dir, out_dir, model="boltz_ode", sampling_steps=2,
              seed=42, step_scale=1.0, preprocessing_threads=1, num_workers=2,
              override=False, python_exe=None):
    py = python_exe or DEFAULT_BOLTZ_PY
    env = dict(os.environ)
    env["PYTHONPATH"] = str(_BOLTZ_SRC) + os.pathsep + env.get("PYTHONPATH", "")
    cmd = [py, str(_BOLTZ_MAIN), "predict", str(input_dir),
           "--out_dir", str(out_dir), "--model", model,
           "--sampling_steps", str(sampling_steps), "--seed", str(seed),
           "--write_full_pae", "--step_scale", str(step_scale),
           "--preprocessing-threads", str(preprocessing_threads),
           "--num_workers", str(num_workers)]
    if override:
        cmd.append("--override")
    print("[predict:boltz] " + " ".join(cmd), file=sys.stderr)
    return subprocess.run(cmd, env=env).returncode


def run_esmfold(input_dir, out_dir, num_loops=3, sampling_steps=50,
                diffusion_samples=1, seed=42, python_exe=None):
    py = python_exe or DEFAULT_ESMFOLD_PY
    cmd = [py, str(_ESMFOLD_WORKER), "--input-dir", str(input_dir),
           "--out-dir", str(out_dir), "--num-loops", str(num_loops),
           "--sampling-steps", str(sampling_steps),
           "--diffusion-samples", str(diffusion_samples), "--seed", str(seed)]
    print("[predict:esmfold] " + " ".join(cmd), file=sys.stderr)
    return subprocess.run(cmd).returncode
