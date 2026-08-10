"""Resolve and install BoltzScan's external command-line tools.

The managed toolchain is deliberately separate from the Python environment
that contains BoltzScan.  This avoids changing a user's active Conda
environment and means callers never need to activate the managed environment;
BoltzScan resolves the executables by their absolute paths.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
from urllib.request import Request, urlopen

from boltzscan import __version__


TOOLCHAIN_VERSION = "2026.08"
TOOLCHAIN_PACKAGES = ("blast=2.16.0", "hmmer=3.4", "meme=5.5.7")
TOOLCHAIN_EXECUTABLES = ("blastp", "makeblastdb", "hmmsearch", "fimo", "tomtom")

# Static, official micromamba release artifacts mirrored from conda-forge.
# Only platforms on which the complete Bioconda toolchain is supported are
# bootstrapped. Other systems can still use an existing conda/mamba executable.
MICROMAMBA_VERSION = "2.8.1-0"
_MICROMAMBA_RELEASES = {
    ("Linux", "x86_64"): (
        "micromamba-linux-64",
        "9689782d863c05a1bf5d2d371ba527104e7a4eb4310c1637d8653b751aed9c82",
    ),
    ("Linux", "aarch64"): (
        "micromamba-linux-aarch64",
        "e5ba23b5945aa49dfd11022e592a510d2686a8feee810e00140b73c9fdf0ba2a",
    ),
}


@dataclass(frozen=True)
class ToolchainInstall:
    root: Path
    prefix: Path
    manager: Path
    manifest: Path
    executables: dict[str, Path]


def default_data_root() -> Path:
    """Return the user-owned BoltzScan data directory without extra packages."""
    override = os.environ.get("BOLTZSCAN_HOME")
    if override:
        return Path(override).expanduser()
    if os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
        return base / "BoltzScan"
    if platform.system() == "Darwin":
        return Path.home() / "Library" / "Application Support" / "BoltzScan"
    base = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))
    return base / "boltzscan"


def managed_tool_root(tool_dir=None) -> Path:
    """Return the versioned managed-tool directory."""
    value = tool_dir or os.environ.get("BOLTZSCAN_TOOL_DIR")
    if value:
        return Path(value).expanduser()
    return default_data_root() / "toolchains" / TOOLCHAIN_VERSION


def managed_prefix(tool_dir=None) -> Path:
    return managed_tool_root(tool_dir) / "env"


def managed_binary(name: str, tool_dir=None) -> Path:
    prefix = managed_prefix(tool_dir)
    if os.name == "nt":
        return prefix / "Scripts" / f"{name}.exe"
    return prefix / "bin" / name


def _resolve_explicit(value) -> str | None:
    if value is None:
        return None
    candidate = Path(value).expanduser()
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate.resolve())
    return shutil.which(str(value))


def resolve_executable(name: str, *, explicit=None, tool_dir=None) -> str | None:
    """Resolve explicit path, managed tool, then the process ``PATH``."""
    if explicit is not None:
        return _resolve_explicit(explicit)
    managed = managed_binary(name, tool_dir)
    if managed.is_file() and os.access(managed, os.X_OK):
        return str(managed.resolve())
    return shutil.which(name)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_verified(url: str, destination: Path, expected_sha256: str, report) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": f"BoltzScan/{__version__}"})
    with tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)
        try:
            report(f"Downloading {url}")
            with urlopen(request, timeout=600) as response:
                shutil.copyfileobj(response, temporary)
            temporary.flush()
            os.fsync(temporary.fileno())
            observed = _sha256(temporary_path)
            if observed.lower() != expected_sha256.lower():
                raise ValueError(
                    "micromamba SHA256 mismatch: "
                    f"expected {expected_sha256}, observed {observed}"
                )
            temporary_path.chmod(0o755)
            temporary_path.replace(destination)
        except Exception:
            temporary_path.unlink(missing_ok=True)
            raise


def bootstrap_micromamba(tool_dir=None, *, report=print) -> Path:
    """Install one checksum-pinned micromamba binary under the managed root."""
    root = managed_tool_root(tool_dir)
    executable = root / "bootstrap" / "micromamba"
    if executable.is_file() and os.access(executable, os.X_OK):
        return executable

    release = _MICROMAMBA_RELEASES.get((platform.system(), platform.machine()))
    if release is None:
        raise RuntimeError(
            "No automatic micromamba bootstrap is available for "
            f"{platform.system()} {platform.machine()}. Install conda, mamba, or "
            "micromamba once, then rerun `boltzscan doctor --fix`."
        )
    asset, digest = release
    url = (
        "https://github.com/mamba-org/micromamba-releases/releases/download/"
        f"{MICROMAMBA_VERSION}/{asset}"
    )
    _download_verified(url, executable, digest, report)
    return executable


def find_package_manager(tool_dir=None, *, bootstrap=False, report=print) -> Path | None:
    """Find a Conda-compatible manager without modifying its current environment."""
    bootstrapped = managed_tool_root(tool_dir) / "bootstrap" / "micromamba"
    if bootstrapped.is_file() and os.access(bootstrapped, os.X_OK):
        return bootstrapped
    for name in ("micromamba", "mamba", "conda"):
        found = shutil.which(name)
        if found:
            return Path(found).resolve()
    if bootstrap:
        return bootstrap_micromamba(tool_dir, report=report)
    return None


def _manager_command(manager: Path, action: str, prefix: Path, packages) -> list[str]:
    return [
        str(manager),
        action,
        "--yes",
        "--prefix",
        str(prefix),
        "--override-channels",
        "--channel",
        "conda-forge",
        "--channel",
        "bioconda",
        "--strict-channel-priority",
        *packages,
    ]


def _probe_version(name: str, executable: Path) -> str:
    flags = {
        "blastp": ("-version",),
        "makeblastdb": ("-version",),
        "hmmsearch": ("-h",),
        "fimo": ("--version",),
        "tomtom": ("-version",),
    }.get(name, ("--version",))
    try:
        result = subprocess.run(
            [str(executable), *flags],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    output = result.stdout or result.stderr
    return next((line.strip() for line in output.splitlines() if line.strip()), "unknown")


def install_toolchain(tool_dir=None, *, report=print) -> ToolchainInstall:
    """Create or update a private, pinned external-tool environment."""
    root = managed_tool_root(tool_dir)
    prefix = managed_prefix(tool_dir)
    root.mkdir(parents=True, exist_ok=True)
    manager = find_package_manager(tool_dir, bootstrap=True, report=report)
    action = "install" if (prefix / "conda-meta" / "history").is_file() else "create"
    command = _manager_command(manager, action, prefix, TOOLCHAIN_PACKAGES)
    report(f"Installing BoltzScan toolchain in {prefix} with {manager.name}")
    environment = os.environ.copy()
    environment.setdefault("MAMBA_ROOT_PREFIX", str(root / "mamba-root"))
    try:
        subprocess.run(command, check=True, env=environment)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"{manager.name} could not install the BoltzScan toolchain "
            f"(exit code {exc.returncode})"
        ) from exc

    executables = {
        name: managed_binary(name, tool_dir)
        for name in TOOLCHAIN_EXECUTABLES
    }
    missing = [
        name for name, path in executables.items()
        if not path.is_file() or not os.access(path, os.X_OK)
    ]
    if missing:
        raise RuntimeError(
            f"Toolchain installation completed but these executables are missing: "
            f"{', '.join(missing)}"
        )

    manifest = root / "install_manifest.json"
    manifest.write_text(json.dumps({
        "schema_version": 1,
        "toolchain_version": TOOLCHAIN_VERSION,
        "boltzscan_version": __version__,
        "installed_at_utc": datetime.now(timezone.utc).isoformat(),
        "manager": str(manager),
        "prefix": str(prefix),
        "packages": list(TOOLCHAIN_PACKAGES),
        "executables": {
            name: {
                "path": str(path),
                "version": _probe_version(name, path),
            }
            for name, path in executables.items()
        },
    }, indent=2) + "\n")
    return ToolchainInstall(root, prefix, manager, manifest, executables)


def executable_version(name: str, executable) -> str:
    return _probe_version(name, Path(executable))
