"""Run protein-DNA structure prediction by selecting one public model.

Inference uses this repository's ``boltz/`` and ``esm/`` source trees. By
default they run with the active BoltzScan Python; environment-variable
overrides remain available when dependency isolation is needed.

Each engine keeps its native output layout. The separate ``wash`` command can
publish a method-named soft or hard view without altering these source trees.
``inference_parameters.json`` is written beside the native prediction folder.
Env paths are overridable via the BOLTZSCAN_* environment variables.
"""
import csv
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

_PKG = Path(__file__).resolve().parents[1]          # boltzscan/
_BOLTZ_SRC = _PKG / "boltz" / "src"
_ESM_SRC = _PKG / "esm"
_BOLTZ_WORKER = Path(__file__).resolve().parent / "_boltz_worker.py"
_ESMFOLD_WORKER = Path(__file__).resolve().parent / "_esmfold_worker.py"


PREDICTION_MODELS = (
    "esmfold2", "boltz1", "boltz2", "boltz1_ode", "boltz2_ode",
)
_BOLTZ_MODELS = frozenset(PREDICTION_MODELS[1:])
_ODE_DEFAULT_SAMPLING_STEPS = 2
_ODE_DEFAULT_STEP_SCALE = 1.0
_BOLTZ_NUM_WORKERS = 2
_SEQUENCE_ALPHABETS = {
    "protein": frozenset("ARNDCEQGHILKMFPSTWYVXJBZOU"),
    "dna": frozenset("ACGTN"),
}
_RUNTIME_MODULES = {
    "boltz": ("torch", "pytorch_lightning", "rdkit"),
    "esmfold2": ("torch", "transformers", "biotite"),
}


def _runtime_python(explicit, env_var, runtime):
    """Use an explicit interpreter, otherwise the active Python with its extra."""
    configured = explicit or os.environ.get(env_var)
    if configured:
        return configured
    missing = [name for name in _RUNTIME_MODULES[runtime]
               if importlib.util.find_spec(name) is None]
    if missing:
        raise ModuleNotFoundError(
            f"{runtime} inference dependencies are missing from {sys.executable}: "
            f"{', '.join(missing)}. Install the `{runtime}` extra in this environment "
            f"or set {env_var} to another Python interpreter."
        )
    return sys.executable


def engine_for_model(model):
    """Resolve the implementation engine from the single public model name."""
    if model == "esmfold2":
        return "esmfold"
    if model in _BOLTZ_MODELS:
        return "boltz"
    raise ValueError(f"Unsupported prediction model: {model}")


def prediction_root_name(model):
    engine_for_model(model)
    return f"{model}_prediction"


def available_cpu_count():
    """Return CPUs available to this process, respecting Linux CPU affinity."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def native_prediction_dir(input_dir, out_dir, model):
    """Return the engine-native prediction directory for one input arm."""
    arm = Path(input_dir).stem
    prefix = "esmfold_results_" if engine_for_model(model) == "esmfold" else "boltz_results_"
    return Path(out_dir) / f"{prefix}{arm}" / "predictions"


def _msa_query_sequence(path):
    """Return the first MSA sequence, stripped as Boltz interprets A3M rows."""
    path = Path(path)
    try:
        if path.suffix == ".a3m":
            with path.open() as handle:
                sequences = (
                    line.strip() for line in handle
                    if line.strip() and not line.startswith((">", "#"))
                )
                rows = list(sequences)
        elif path.suffix == ".csv":
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                if set(reader.fieldnames or ()) != {"key", "sequence"}:
                    raise ValueError(f"Invalid MSA CSV columns in {path}; expected key,sequence")
                rows = [str(row["sequence"]).strip() for row in reader if row["sequence"]]
        else:
            raise ValueError(f"Unsupported MSA file {path}; expected .a3m or .csv")
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f"Cannot read MSA {path}: {exc}") from None
    if not rows:
        raise ValueError(f"MSA contains no sequences: {path}")
    for row in rows:
        invalid = sorted({c for c in row if c != "-" and c.upper() not in _SEQUENCE_ALPHABETS["protein"]})
        if invalid:
            raise ValueError(f"Illegal MSA character(s) {invalid} in {path}")
    return "".join(c for c in rows[0] if c != "-" and not c.islower()).upper()


def preflight_prediction_inputs(input_dir, model):
    """Validate TF-DNA YAMLs and report their MSA mode before inference."""
    engine = engine_for_model(model)
    input_dir = Path(input_dir)
    yamls = sorted(input_dir.glob("*.yaml")) if input_dir.is_dir() else []
    if not yamls:
        raise ValueError(f"No .yaml prediction inputs found under {input_dir}")

    proteins = proteins_with_msa = 0
    msa_queries = {}
    for yaml_path in yamls:
        try:
            spec = yaml.safe_load(yaml_path.read_text())
        except (OSError, UnicodeError, yaml.YAMLError) as exc:
            raise ValueError(f"Cannot read prediction YAML {yaml_path}: {exc}") from None
        if not isinstance(spec, dict) or spec.get("version", 1) != 1:
            raise ValueError(f"Invalid version/schema in prediction YAML: {yaml_path}")
        entries = spec.get("sequences")
        if not isinstance(entries, list) or not entries:
            raise ValueError(f"Missing non-empty sequences list in {yaml_path}")

        seen_ids = set()
        entity_counts = {"protein": 0, "dna": 0}
        dna_sequences = []
        for entry in entries:
            if not isinstance(entry, dict) or len(entry) != 1:
                raise ValueError(f"Each sequence entry must contain one entity in {yaml_path}")
            entity = next(iter(entry))
            if entity not in _SEQUENCE_ALPHABETS or not isinstance(entry[entity], dict):
                raise ValueError(f"BoltzScan predict supports only protein and DNA in {yaml_path}")
            payload = entry[entity]
            ids = payload.get("id")
            ids = [ids] if isinstance(ids, str) else ids
            if not isinstance(ids, list) or not ids or not all(isinstance(i, str) and i for i in ids):
                raise ValueError(f"Invalid {entity} chain id in {yaml_path}")
            if len(set(ids)) != len(ids) or seen_ids.intersection(ids):
                raise ValueError(f"Duplicate chain id in {yaml_path}: {ids}")
            seen_ids.update(ids)

            sequence = payload.get("sequence")
            if not isinstance(sequence, str) or not sequence:
                raise ValueError(f"Missing {entity} sequence for chain {ids} in {yaml_path}")
            invalid = sorted(set(sequence).difference(_SEQUENCE_ALPHABETS[entity]))
            if invalid:
                raise ValueError(
                    f"Illegal {entity} character(s) {invalid} for chain {ids} in {yaml_path}"
                )
            entity_counts[entity] += len(ids)

            if entity != "protein":
                dna_sequences.extend([sequence] * len(ids))
                continue
            proteins += len(ids)
            msa = payload.get("msa")
            has_msa = isinstance(msa, str) and msa not in {"", "empty"}
            if not has_msa:
                if engine == "boltz":
                    raise ValueError(
                        f"{model} requires a local MSA for protein chain {ids} in {yaml_path}"
                    )
                continue
            msa_path = Path(msa)
            if not msa_path.is_file():
                raise ValueError(f"MSA file not found for protein chain {ids}: {msa_path}")
            if engine == "esmfold" and msa_path.suffix != ".a3m":
                raise ValueError(f"ESMFold2 MSA must be .a3m: {msa_path}")
            if msa_path not in msa_queries:
                msa_queries[msa_path] = _msa_query_sequence(msa_path)
            query = msa_queries[msa_path]
            if query != sequence:
                raise ValueError(
                    f"Protein/MSA query mismatch in {yaml_path} chain {ids}: "
                    f"protein length {len(sequence)}, MSA query length {len(query)} ({msa_path})"
                )
            proteins_with_msa += len(ids)

        if entity_counts != {"protein": 1, "dna": 2}:
            raise ValueError(
                f"Prediction YAML must contain 1 protein and 2 DNA chains: {yaml_path}"
            )
        reverse_complement = dna_sequences[0].translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
        if dna_sequences[1] != reverse_complement:
            raise ValueError(f"DNA chains are not reverse complements in {yaml_path}")

    msa_mode = (
        "MSA" if proteins_with_msa == proteins
        else "no-MSA" if proteins_with_msa == 0
        else "mixed"
    )
    print(
        f"Prediction preflight: model={model}; inputs={len(yamls)}; "
        f"proteins={proteins}; MSA mode={msa_mode} ({proteins_with_msa}/{proteins}); "
        "sequences valid"
    )
    return msa_mode


def _write_inference_parameters(metadata_dir, *, engine, model, input_dir,
                                output_dir, parameters):
    metadata_dir = Path(metadata_dir)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "engine": engine,
        "model": model,
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "msa_source": "input_yaml",
        "parameters": parameters,
    }
    (metadata_dir / "inference_parameters.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )


def _runtime_environment(python_exe, source_dir, source_package, model):
    """Validate an inference runtime and put its local source first."""
    python_path = Path(python_exe)
    if not python_path.is_file():
        raise FileNotFoundError(
            f"{model} Python interpreter not found: {python_path}. "
            f"Set BOLTZSCAN_{'BOLTZ' if model.startswith('boltz') else 'ESMFOLD'}_PYTHON."
        )
    if not (source_dir / source_package).is_dir():
        raise FileNotFoundError(f"Local {model} source tree not found: {source_dir}")
    env = dict(os.environ)
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        str(source_dir) if not existing else str(source_dir) + os.pathsep + existing
    )
    return env


def run_boltz(input_dir, out_dir, model="boltz1_ode", seed=None, python_exe=None):
    if model not in _BOLTZ_MODELS:
        raise ValueError(f"Unsupported Boltz model: {model}")

    if model in {"boltz1_ode", "boltz2_ode"}:
        command_sampling_steps = _ODE_DEFAULT_SAMPLING_STEPS
        command_step_scale = _ODE_DEFAULT_STEP_SCALE
        configuration_source = f"boltzscan_{model}"
    else:
        command_sampling_steps = None
        command_step_scale = None
        configuration_source = "boltz_native_default"

    preprocessing_threads = available_cpu_count()
    py = _runtime_python(
        python_exe, "BOLTZSCAN_BOLTZ_PYTHON", "boltz"
    )
    env = _runtime_environment(py, _BOLTZ_SRC, "boltz", model)
    native_result_dir = Path(out_dir) / f"boltz_results_{Path(input_dir).stem}"
    _write_inference_parameters(
        native_result_dir,
        engine="boltz",
        model=model,
        input_dir=input_dir,
        output_dir=native_result_dir / "predictions",
        parameters={
            "configuration_source": configuration_source,
            "sampling_steps": command_sampling_steps,
            "sampling_steps_source": configuration_source,
            "seed": seed,
            "seed_source": "explicit" if seed is not None else "native_default",
            "step_scale": command_step_scale,
            "step_scale_source": configuration_source,
            "write_full_pae": True,
            "preprocessing_threads": preprocessing_threads,
            "preprocessing_threads_source": "available_cpu_affinity",
            "num_workers": _BOLTZ_NUM_WORKERS,
            "num_workers_source": "boltzscan_internal_default",
            "override": False,
        },
    )
    input_dir = Path(input_dir)
    with tempfile.TemporaryDirectory(prefix="boltzscan-inputs-") as tmp:
        boltz_input = input_dir
        if input_dir.is_dir() and any(
            child.suffix.lower() not in {".yml", ".yaml"}
            for child in input_dir.iterdir()
        ):
            boltz_input = Path(tmp) / input_dir.name
            boltz_input.mkdir()
            for yaml_path in sorted(input_dir.glob("*.yaml")):
                (boltz_input / yaml_path.name).symlink_to(yaml_path.resolve())

        cmd = [py, str(_BOLTZ_WORKER), "predict", str(boltz_input),
               "--out_dir", str(out_dir), "--model", model,
               "--write_full_pae",
               "--preprocessing-threads", str(preprocessing_threads),
               "--num_workers", str(_BOLTZ_NUM_WORKERS)]
        if seed is not None:
            cmd.extend(["--seed", str(seed)])
        if command_sampling_steps is not None:
            cmd.extend(["--sampling_steps", str(command_sampling_steps)])
        if command_step_scale is not None:
            cmd.extend(["--step_scale", str(command_step_scale)])
        print("[predict:boltz] " + " ".join(cmd), file=sys.stderr)
        return subprocess.run(cmd, env=env).returncode


def run_esmfold(input_dir, out_dir, seed=None, python_exe=None):
    py = _runtime_python(
        python_exe, "BOLTZSCAN_ESMFOLD_PYTHON", "esmfold2"
    )
    env = _runtime_environment(py, _ESM_SRC, "esm", "esmfold2")
    esmfold2_weights = env.get("BOLTZSCAN_ESMFOLD2_WEIGHTS", "biohub/ESMFold2")
    esmc_weights = env.get("BOLTZSCAN_ESMC_WEIGHTS", "biohub/ESMC-6B")
    ccd_path = env.get("ESMCFOLD_CCD_PATH")
    if ccd_path is None:
        local_ccd = Path(esmfold2_weights).expanduser() / "ccd.pkl"
        if local_ccd.is_file():
            ccd_path = str(local_ccd)
            env["ESMCFOLD_CCD_PATH"] = ccd_path
    _write_inference_parameters(
        out_dir,
        engine="esmfold",
        model="esmfold2",
        input_dir=input_dir,
        output_dir=Path(out_dir) / "predictions",
        parameters={
            "configuration_source": "esmfold2_native_default",
            "seed": seed,
            "seed_source": "explicit" if seed is not None else "native_default",
            "deterministic": seed is not None,
            "deterministic_source": (
                "explicit_seed" if seed is not None else "native_default"
            ),
            "esmfold2_weights": esmfold2_weights,
            "esmc_weights": esmc_weights,
            "ccd_path": ccd_path or "biohub/ESMFold2:ccd.pkl",
        },
    )
    cmd = [py, str(_ESMFOLD_WORKER), "--input-dir", str(input_dir),
           "--out-dir", str(out_dir),
           "--weights", esmfold2_weights,
           "--esmc-weights", esmc_weights]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    print("[predict:esmfold] " + " ".join(cmd), file=sys.stderr)
    return subprocess.run(cmd, env=env).returncode


def run_prediction(input_dir, out_dir, model="esmfold2", seed=None):
    """Run one input directory through the selected public model."""
    preflight_prediction_inputs(input_dir, model)
    if engine_for_model(model) == "esmfold":
        native_root = native_prediction_dir(input_dir, out_dir, model).parent
        return run_esmfold(input_dir=input_dir, out_dir=native_root, seed=seed)
    return run_boltz(input_dir=input_dir, out_dir=out_dir, model=model, seed=seed)
