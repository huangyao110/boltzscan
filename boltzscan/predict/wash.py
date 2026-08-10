"""Publish native prediction outputs in the shared BoltzScan layout."""

import errno
import filecmp
import json
import os
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from boltzscan.predict.runners import engine_for_model, prediction_root_name


@dataclass(frozen=True)
class WashSummary:
    method_root: Path
    manifest: Path
    mode: str
    arms: tuple[str, ...]
    files: int
    hard_linked: int
    copied: int


def _add_source(sources, arm, source):
    source = Path(source)
    previous = sources.get(arm)
    if previous is not None and previous.resolve() != source.resolve():
        raise ValueError(
            f"Multiple native prediction directories found for arm {arm}: "
            f"{previous} and {source}"
        )
    sources[arm] = source


def discover_native_predictions(run_dir, model):
    """Return ``arm -> native predictions directory`` without changing it."""
    run_dir = Path(run_dir)
    engine = engine_for_model(model)

    prefix = "boltz_results_" if engine == "boltz" else "esmfold_results_"
    sources = {}
    for result_dir in sorted(run_dir.glob(f"{prefix}*")):
        source = result_dir / "predictions"
        if result_dir.is_dir() and source.is_dir():
            _add_source(sources, result_dir.name[len(prefix):], source)

    return sources


def _tree_counts(root):
    files = directories = 0
    for _current, dirnames, filenames in os.walk(root):
        directories += len(dirnames)
        files += len(filenames)
    return files, directories


def _soft_publish(source, destination):
    if destination.is_symlink():
        if destination.resolve() != source.resolve():
            raise FileExistsError(
                f"Wash destination points to a different source: {destination}"
            )
        return "reused"
    if destination.exists():
        raise FileExistsError(
            f"Wash destination is already a real path; use --mode hard to update it: "
            f"{destination}"
        )
    destination.symlink_to(
        os.path.relpath(source.resolve(), destination.parent.resolve()),
        target_is_directory=True,
    )
    return "created"


def _link_or_copy(source, destination):
    try:
        os.link(source, destination)
        return "hard_linked"
    except OSError as exc:
        fallback_errors = {
            errno.EXDEV,
            errno.EPERM,
            errno.EACCES,
            getattr(errno, "ENOTSUP", errno.EPERM),
            getattr(errno, "EOPNOTSUPP", errno.EPERM),
        }
        if exc.errno not in fallback_errors:
            raise
    shutil.copy2(source, destination)
    return "copied"


def _hard_publish(source, destination):
    if destination.is_symlink():
        if destination.resolve() != source.resolve():
            raise FileExistsError(
                f"Refusing to replace an unrelated wash symlink: {destination}"
            )
        destination.unlink()
    elif destination.exists() and not destination.is_dir():
        raise FileExistsError(f"Wash destination is not a directory: {destination}")
    destination.mkdir(parents=True, exist_ok=True)

    stats = {"hard_linked": 0, "copied": 0, "reused": 0}
    for current, dirnames, filenames in os.walk(source):
        relative = Path(current).relative_to(source)
        target_dir = destination / relative
        target_dir.mkdir(parents=True, exist_ok=True)
        for dirname in dirnames:
            (target_dir / dirname).mkdir(exist_ok=True)
        for filename in filenames:
            source_file = Path(current) / filename
            target_file = target_dir / filename
            if target_file.exists():
                if source_file.samefile(target_file):
                    stats["reused"] += 1
                    continue
                if target_file.is_file() and filecmp.cmp(
                    source_file, target_file, shallow=False
                ):
                    stats["reused"] += 1
                    continue
                raise FileExistsError(
                    f"Wash destination contains different data: {target_file}"
                )
            result = _link_or_copy(source_file, target_file)
            stats[result] += 1
    return stats


def _load_inference_parameters(source):
    candidates = (
        source.parent / "inference_parameters.json",
        source / "inference_parameters.json",
        source.parent.parent / "inference_parameters.json",
    )
    for candidate in candidates:
        if candidate.is_file():
            try:
                return str(candidate.resolve()), json.loads(candidate.read_text())
            except (OSError, json.JSONDecodeError):
                return str(candidate.resolve()), None
    return None, None


def wash_prediction_outputs(
    run_dir,
    *,
    model="esmfold2",
    mode="soft",
    arms=None,
):
    """Create a tracked normalized view while preserving all native outputs.

    ``soft`` creates relative directory symlinks and stays live while inference
    runs. ``hard`` creates real directories using hard links when possible and
    copies only when the filesystem cannot hard-link; rerun it to refresh a
    snapshot after more native predictions appear.
    """
    if mode not in {"soft", "hard"}:
        raise ValueError(f"Unsupported wash mode: {mode}")
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    engine = engine_for_model(model)
    sources = discover_native_predictions(run_dir, model)
    if arms is not None:
        requested = tuple(dict.fromkeys(arms))
        missing = [arm for arm in requested if arm not in sources]
        if missing:
            raise FileNotFoundError(
                f"Native prediction directories not found for arms: {', '.join(missing)}"
            )
        sources = {arm: sources[arm] for arm in requested}
    if not sources:
        raise FileNotFoundError(
            f"No native {engine} prediction directories found under {run_dir}"
        )

    method_model = model
    source_metadata = {}
    for arm, source in sources.items():
        parameters_file, parameters = _load_inference_parameters(source)
        recorded_model = parameters.get("model") if isinstance(parameters, dict) else None
        if recorded_model is not None and recorded_model != method_model:
            raise ValueError(
                f"Native {arm} metadata records model={recorded_model}, "
                f"not requested model={method_model}"
            )
        source_metadata[arm] = (parameters_file, parameters)

    method_root = run_dir / prediction_root_name(model)
    method_root.mkdir(parents=True, exist_ok=True)
    arm_manifest = {}
    total_files = total_hard_linked = total_copied = 0

    for arm, source in sorted(sources.items()):
        destination = method_root / arm
        source_files, source_directories = _tree_counts(source)
        total_files += source_files
        parameters_file, parameters = source_metadata[arm]
        if mode == "soft":
            status = _soft_publish(source, destination)
            stats = {"status": status, "hard_linked": 0, "copied": 0, "reused": 0}
        else:
            stats = _hard_publish(source, destination)
            stats["status"] = "materialized"
            total_hard_linked += stats["hard_linked"]
            total_copied += stats["copied"]

        arm_manifest[arm] = {
            "source": str(source.resolve()),
            "destination": str(destination.absolute()),
            "source_preserved": True,
            "source_files": source_files,
            "source_directories": source_directories,
            "publish": stats,
            "inference_parameters_file": parameters_file,
            "inference_parameters": parameters,
        }

    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(run_dir.resolve()),
        "engine": engine,
        "model": method_model,
        "mode": mode,
        "source_policy": "preserve",
        "hard_mode_policy": "hard-link files; copy only when hard links are unavailable",
        "arms": arm_manifest,
    }
    manifest_path = method_root / "wash_manifest.json"
    temporary_manifest = method_root / ".wash_manifest.json.tmp"
    temporary_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    temporary_manifest.replace(manifest_path)
    return WashSummary(
        method_root=method_root,
        manifest=manifest_path,
        mode=mode,
        arms=tuple(sorted(sources)),
        files=total_files,
        hard_linked=total_hard_linked,
        copied=total_copied,
    )
