"""Build and resolve the small Pfam library required by BoltzScan.

The complete Pfam-A library is a maintainer input.  Runtime TF discovery and
PWM transfer need only the profiles referenced by the plant-TF classification
rules, so release artifacts carry this deterministic subset instead.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tempfile


RUNTIME_PFAM_DIR = "pfam"
RUNTIME_PFAM_NAME = "plant_tf_pfam.hmm"
RUNTIME_PFAM_MANIFEST = "plant_tf_pfam_manifest.json"
DEFAULT_FULL_PFAM = Path(__file__).resolve().parents[2] / "data" / "pfam" / "Pfam-A.hmm"


@dataclass(frozen=True)
class RuntimePfamSummary:
    hmm: Path
    manifest: Path
    n_profiles: int
    size_bytes: int
    sha256: str


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def required_pfam_accessions() -> frozenset[str]:
    """Return every Pfam accession used by TF and TE-exclusion rules."""
    from boltzscan.utils.find_tf import DBD_ACCS, RULES, TE_FORBIDDEN

    accessions = set(DBD_ACCS) | set(TE_FORBIDDEN)

    def collect(value):
        if isinstance(value, str) and value.startswith("PF"):
            accessions.add(value.split(".")[0])
        elif isinstance(value, (list, tuple, set, frozenset)):
            for item in value:
                collect(item)

    for _family, all_of, none_of in RULES:
        collect(all_of)
        collect(none_of)
    return frozenset(accessions)


def runtime_pfam_paths(refs_dir: str | Path) -> tuple[Path, Path]:
    root = Path(refs_dir) / RUNTIME_PFAM_DIR
    return root / RUNTIME_PFAM_NAME, root / RUNTIME_PFAM_MANIFEST


def packaged_runtime_pfam_paths() -> tuple[Path, Path]:
    root = Path(__file__).resolve().parents[1] / "data" / "pfam"
    return root / RUNTIME_PFAM_NAME, root / RUNTIME_PFAM_MANIFEST


def _open_hmm(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def _profile_accessions(path: str | Path) -> tuple[str, ...]:
    accessions = []
    with _open_hmm(Path(path)) as handle:
        for line in handle:
            if line.startswith("ACC   "):
                accessions.append(line.split()[1].split(".")[0])
    return tuple(accessions)


def validate_runtime_pfam(
    hmm: str | Path,
    manifest: str | Path | None = None,
) -> RuntimePfamSummary:
    """Validate profile membership and, when present, its release manifest."""
    hmm = Path(hmm)
    if not hmm.is_file():
        raise FileNotFoundError(f"BoltzScan Pfam HMM not found: {hmm}")
    accessions = _profile_accessions(hmm)
    observed = set(accessions)
    required = set(required_pfam_accessions())
    missing = sorted(required - observed)
    extra = sorted(observed - required)
    duplicates = sorted(acc for acc, count in Counter(accessions).items() if count > 1)
    if missing or extra or duplicates or len(accessions) != len(required):
        details = []
        if missing:
            details.append(f"missing={','.join(missing[:10])}")
        if extra:
            details.append(f"extra={','.join(extra[:10])}")
        if duplicates:
            details.append(f"duplicates={','.join(duplicates[:10])}")
        raise ValueError(
            f"Invalid BoltzScan Pfam subset {hmm}: " + "; ".join(details)
        )

    digest = sha256_file(hmm)
    manifest_path = Path(manifest) if manifest is not None else hmm.with_name(
        RUNTIME_PFAM_MANIFEST
    )
    if not manifest_path.is_file():
        raise FileNotFoundError(f"BoltzScan Pfam manifest not found: {manifest_path}")
    try:
        metadata = json.loads(manifest_path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid BoltzScan Pfam manifest: {manifest_path}") from exc
    subset = metadata.get("subset", {})
    if metadata.get("schema_version") != 1:
        raise ValueError(f"Unsupported BoltzScan Pfam manifest: {manifest_path}")
    if subset.get("sha256") != digest:
        raise ValueError(
            f"BoltzScan Pfam SHA256 mismatch: manifest={subset.get('sha256')}, file={digest}"
        )
    if subset.get("profile_count") != len(accessions):
        raise ValueError(f"BoltzScan Pfam profile count mismatch: {manifest_path}")
    if sorted(subset.get("accessions", [])) != sorted(required):
        raise ValueError(f"BoltzScan Pfam accession manifest mismatch: {manifest_path}")
    return RuntimePfamSummary(
        hmm=hmm,
        manifest=manifest_path,
        n_profiles=len(accessions),
        size_bytes=hmm.stat().st_size,
        sha256=digest,
    )


def extract_runtime_pfam(
    source_pfam: str | Path,
    destination_dir: str | Path,
) -> RuntimePfamSummary:
    """Extract the exact plant-TF rule profiles from a full Pfam-A library."""
    source_pfam = Path(source_pfam)
    if not source_pfam.is_file():
        raise FileNotFoundError(f"Full Pfam-A HMM not found: {source_pfam}")
    destination_dir = Path(destination_dir)
    destination_dir.mkdir(parents=True, exist_ok=True)
    output = destination_dir / RUNTIME_PFAM_NAME
    manifest = destination_dir / RUNTIME_PFAM_MANIFEST
    required = set(required_pfam_accessions())
    selected = set()
    total_profiles = 0
    first_profile_date = None

    with tempfile.NamedTemporaryFile(
        mode="w",
        dir=destination_dir,
        prefix=f".{RUNTIME_PFAM_NAME}.",
        suffix=".tmp",
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)
        try:
            block = []
            accession = None
            profile_date = None
            with _open_hmm(source_pfam) as source:
                for line in source:
                    block.append(line)
                    if line.startswith("ACC   "):
                        accession = line.split()[1].split(".")[0]
                    elif line.startswith("DATE  "):
                        profile_date = line.removeprefix("DATE  ").strip()
                    if line.rstrip() != "//":
                        continue
                    total_profiles += 1
                    if accession in required:
                        if accession in selected:
                            raise ValueError(
                                f"Duplicate required Pfam accession in {source_pfam}: {accession}"
                            )
                        temporary.writelines(block)
                        selected.add(accession)
                        if first_profile_date is None:
                            first_profile_date = profile_date
                    block = []
                    accession = None
                    profile_date = None
            if block:
                raise ValueError(f"Truncated Pfam HMM block in {source_pfam}")
        except Exception:
            temporary_path.unlink(missing_ok=True)
            raise

    missing = sorted(required - selected)
    if missing:
        temporary_path.unlink(missing_ok=True)
        raise ValueError(
            f"Full Pfam library is missing {len(missing)} required profiles: "
            + ", ".join(missing[:10])
        )
    temporary_path.replace(output)
    output.chmod(0o644)
    digest = sha256_file(output)
    metadata = {
        "schema_version": 1,
        "artifact": "boltzscan-plant-tf-pfam",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "filename": source_pfam.name,
            "size_bytes": source_pfam.stat().st_size,
            "sha256": sha256_file(source_pfam),
            "profile_count": total_profiles,
            "first_selected_profile_date": first_profile_date,
        },
        "subset": {
            "filename": RUNTIME_PFAM_NAME,
            "profile_count": len(selected),
            "size_bytes": output.stat().st_size,
            "sha256": digest,
            "accessions": sorted(selected),
        },
        "selection": (
            "All Pfam accessions referenced by BoltzScan plant-TF family rules, "
            "DBD mapping, and transposon-exclusion rules."
        ),
        "license": "Pfam data are distributed under CC0.",
        "upstream": "https://ftp.ebi.ac.uk/pub/databases/Pfam/",
    }
    manifest.write_text(json.dumps(metadata, indent=2) + "\n")
    manifest.chmod(0o644)
    return validate_runtime_pfam(output, manifest)


def install_packaged_runtime_pfam(refs_dir: str | Path) -> RuntimePfamSummary:
    """Copy the wheel's pinned subset into a reference store when necessary."""
    source_hmm, source_manifest = packaged_runtime_pfam_paths()
    packaged = validate_runtime_pfam(source_hmm, source_manifest)
    target_hmm, target_manifest = runtime_pfam_paths(refs_dir)
    if target_hmm.exists() or target_manifest.exists():
        return validate_runtime_pfam(target_hmm, target_manifest)
    target_hmm.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(packaged.hmm, target_hmm)
    shutil.copy2(packaged.manifest, target_manifest)
    return validate_runtime_pfam(target_hmm, target_manifest)


def resolve_runtime_pfam(
    refs_dir: str | Path = "data/pwms/_refs",
    explicit: str | Path | None = None,
) -> Path:
    """Resolve explicit HMM, reference-owned subset, then packaged subset."""
    if explicit is not None:
        path = Path(explicit).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"Pfam HMM not found: {path}")
        return path.resolve()

    refs_hmm, refs_manifest = runtime_pfam_paths(refs_dir)
    if refs_hmm.exists() or refs_manifest.exists():
        return validate_runtime_pfam(refs_hmm, refs_manifest).hmm.resolve()
    packaged_hmm, packaged_manifest = packaged_runtime_pfam_paths()
    try:
        return validate_runtime_pfam(packaged_hmm, packaged_manifest).hmm.resolve()
    except (FileNotFoundError, ValueError) as exc:
        raise FileNotFoundError(
            "BoltzScan's plant-TF Pfam subset is unavailable. Run "
            f"`boltzscan install-pwm-refs --refs {refs_dir} --replace` or pass --pfam."
        ) from exc
