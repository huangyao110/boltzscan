"""Read-only environment checks and explicit repair for BoltzScan."""

from __future__ import annotations

from dataclasses import dataclass
import importlib.util
from pathlib import Path
import platform
import sys

from boltzscan.pwmmap.archive import validate_runtime_store
from boltzscan.pwmmap.dbd import DEFAULT_PFAM
from boltzscan.toolchain import (
    PROFILE_EXECUTABLES,
    executable_version,
    install_toolchain,
    managed_tool_root,
    resolve_executable,
)


@dataclass(frozen=True)
class DoctorCheck:
    component: str
    status: str
    detail: str
    required: bool = True
    fixable: bool = False


@dataclass(frozen=True)
class DoctorSummary:
    checks: tuple[DoctorCheck, ...]
    tool_root: Path

    @property
    def failed_required(self):
        return tuple(check for check in self.checks if check.required and check.status != "OK")


def _python_checks():
    version = platform.python_version()
    supported = (3, 10) <= sys.version_info[:2] < (3, 13)
    yield DoctorCheck(
        "Python",
        "OK" if supported else "MISSING",
        f"{version} ({sys.executable})" + ("" if supported else "; requires >=3.10,<3.13"),
    )
    imports = {
        "Bio": ("Biopython", "pip install boltzscan"),
        "numpy": ("NumPy", "pip install boltzscan"),
        "pandas": ("pandas", "pip install boltzscan"),
        "yaml": ("PyYAML", "pip install boltzscan"),
        "requests": ("requests", "pip install boltzscan"),
        "tqdm": ("tqdm", "pip install boltzscan"),
    }
    for module, (label, install_command) in imports.items():
        available = importlib.util.find_spec(module) is not None
        hint = "installed" if available else f"missing; run `{install_command}`"
        yield DoctorCheck(label, "OK" if available else "MISSING", hint)


def _vendored_source_checks():
    package_root = Path(__file__).resolve().parent
    sources = {
        "Boltz 2.2.1 source": package_root / "boltz" / "src" / "boltz" / "main.py",
        "ESMFold2 source": package_root / "esm" / "esm" / "__init__.py",
    }
    for label, source in sources.items():
        available = source.is_file()
        detail = (
            str(source)
            if available
            else f"missing from installation: {source}; reinstall the BoltzScan wheel"
        )
        yield DoctorCheck(label, "OK" if available else "MISSING", detail)


def _tool_checks(profile, tool_dir):
    for name in PROFILE_EXECUTABLES[profile]:
        executable = resolve_executable(name, tool_dir=tool_dir)
        if executable is None:
            yield DoctorCheck(
                name,
                "MISSING",
                "not found; run `boltzscan doctor --fix"
                + (" --profile refs-builder" if profile == "refs-builder" else "")
                + "`",
                fixable=True,
            )
            continue
        yield DoctorCheck(name, "OK", f"{executable_version(name, executable)} [{executable}]")


def _data_checks(refs, pfam):
    refs = Path(refs).expanduser()
    try:
        counts = validate_runtime_store(refs)
    except (FileNotFoundError, OSError, ValueError) as exc:
        yield DoctorCheck(
            "PWM references",
            "MISSING",
            f"{exc}; run `boltzscan install-pwm-refs --refs {refs}`",
            fixable=False,
        )
    else:
        yield DoctorCheck(
            "PWM references",
            "OK",
            f"{refs.resolve()} ({counts['n_dbd_rows']} DBD rows, "
            f"{counts['n_meme_motifs']} motifs)",
        )

    pfam = Path(pfam).expanduser()
    if pfam.is_file():
        yield DoctorCheck("Pfam-A.hmm", "OK", f"{pfam.resolve()} ({pfam.stat().st_size} bytes)")
    else:
        yield DoctorCheck(
            "Pfam-A.hmm",
            "MISSING",
            f"not found at {pfam}; pass --pfam or install the Pfam-A HMM database",
            fixable=False,
        )


def inspect_environment(*, profile="runtime", refs="data/pwms/_refs", pfam=None, tool_dir=None):
    if profile not in PROFILE_EXECUTABLES:
        raise ValueError(f"Unknown doctor profile: {profile}")
    checks = tuple([
        *_python_checks(),
        *_vendored_source_checks(),
        *_tool_checks(profile, tool_dir),
        *_data_checks(refs, pfam or DEFAULT_PFAM),
    ])
    return DoctorSummary(checks, managed_tool_root(tool_dir))


def print_summary(summary: DoctorSummary):
    print("BoltzScan doctor")
    print(f"Managed tools: {summary.tool_root}")
    print()
    width = max(len(check.component) for check in summary.checks)
    for check in summary.checks:
        suffix = "" if check.required else " (optional)"
        print(f"{check.component:<{width}}  {check.status:<7}  {check.detail}{suffix}")
    print()
    if summary.failed_required:
        print(f"Result: {len(summary.failed_required)} required check(s) failed")
    else:
        print("Result: environment is ready")


def run_doctor(*, fix=False, profile="runtime", refs="data/pwms/_refs", pfam=None,
               tool_dir=None, report=print):
    if fix:
        install_toolchain(profile=profile, tool_dir=tool_dir, report=report)
    summary = inspect_environment(
        profile=profile,
        refs=refs,
        pfam=pfam,
        tool_dir=tool_dir,
    )
    print_summary(summary)
    return summary
