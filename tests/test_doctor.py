import json
from pathlib import Path

from boltzscan.cli import _build_parser, _command_log_path
from boltzscan.doctor import _vendored_source_checks
from boltzscan import toolchain


def _executable(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\nexit 0\n")
    path.chmod(0o755)
    return path


def test_doctor_cli_is_read_only_by_default(tmp_path):
    args = _build_parser().parse_args([
        "doctor", "--tool-dir", str(tmp_path), "--refs", "shared/refs",
    ])

    assert args.command == "doctor"
    assert args.fix is False
    assert not hasattr(args, "profile")
    assert args.refs == "shared/refs"
    assert _command_log_path(args) == tmp_path / "doctor.log"


def test_managed_executable_precedes_path(tmp_path, monkeypatch):
    managed = _executable(toolchain.managed_binary("blastp", tmp_path))
    monkeypatch.setattr(toolchain.shutil, "which", lambda name: f"/system/{name}")

    assert toolchain.resolve_executable("blastp", tool_dir=tmp_path) == str(managed.resolve())


def test_explicit_executable_precedes_managed(tmp_path):
    _executable(toolchain.managed_binary("blastp", tmp_path))
    explicit = _executable(tmp_path / "custom" / "blastp")

    assert toolchain.resolve_executable(
        "blastp", explicit=explicit, tool_dir=tmp_path,
    ) == str(explicit.resolve())


def test_install_toolchain_uses_isolated_prefix_and_writes_manifest(tmp_path, monkeypatch):
    manager = _executable(tmp_path / "conda")
    calls = []

    monkeypatch.setattr(
        toolchain,
        "find_package_manager",
        lambda *args, **kwargs: manager,
    )
    monkeypatch.setattr(toolchain, "_probe_version", lambda name, path: f"{name} test")

    def fake_run(command, *, check, env):
        calls.append((command, check, env))
        for name in toolchain.TOOLCHAIN_EXECUTABLES:
            _executable(toolchain.managed_binary(name, tmp_path / "tools"))

    monkeypatch.setattr(toolchain.subprocess, "run", fake_run)

    result = toolchain.install_toolchain(
        tmp_path / "tools", report=lambda message: None,
    )

    command, check, environment = calls[0]
    assert check is True
    assert command[1] == "create"
    assert command[command.index("--prefix") + 1] == str(result.prefix)
    assert "blast=2.16.0" in command
    assert "hmmer=3.4" in command
    assert "meme=5.5.7" in command
    assert environment["MAMBA_ROOT_PREFIX"].startswith(str(tmp_path / "tools"))

    manifest = json.loads(result.manifest.read_text())
    assert set(manifest["executables"]) == {
        "blastp", "makeblastdb", "hmmsearch", "fimo", "tomtom",
    }


def test_toolchain_includes_fimo_and_tomtom():
    assert "meme=5.5.7" in toolchain.TOOLCHAIN_PACKAGES
    assert {"fimo", "tomtom"}.issubset(toolchain.TOOLCHAIN_EXECUTABLES)


def test_vendored_inference_sources_are_available():
    checks = tuple(_vendored_source_checks())

    assert {check.component for check in checks} == {
        "Boltz 2.2.1 source", "ESMFold2 source",
    }
    assert all(check.status == "OK" for check in checks)

    manifest = json.loads(
        (Path(toolchain.__file__).parent / "VENDORED_INFERENCE.json").read_text()
    )
    assert manifest["boltz"]["version"] == "2.2.1"
    assert manifest["boltz"]["revision"] == "cb04aeccdd480fd4db707f0bbafde538397fa2ac"
    assert manifest["esm"]["revision"] == "b6b0e8815dc3174a9fe9afd36b222aea57d536c5"
