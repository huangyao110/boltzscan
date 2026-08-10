"""Package and install the small runtime PWM reference store."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile
from urllib.parse import parse_qs, quote, urlparse
from urllib.request import Request, urlopen

from boltzscan import __version__
from boltzscan.pwmmap.pfam import (
    RUNTIME_PFAM_MANIFEST,
    RUNTIME_PFAM_NAME,
    install_packaged_runtime_pfam,
    runtime_pfam_paths,
    validate_runtime_pfam,
)


ARCHIVE_ROOT = '_refs'
DEFAULT_REFERENCE_RELEASE = '2026-08-08'
DEFAULT_REFERENCE_RELEASE_URL = (
    'https://github.com/huangyao110/boltzscan/releases/download/'
    'pwm-refs-2026-08-08/boltzscan_pwm_refs_20260808.tar.gz'
)
DEFAULT_REFERENCE_RELEASE_CHECKSUM_URL = (
    'https://github.com/huangyao110/boltzscan/releases/download/'
    'pwm-refs-2026-08-08/boltzscan_pwm_refs_20260808.tar.gz.sha256'
)
DEFAULT_REFERENCE_RELEASE_SHA256 = (
    'be71be84ff7467d2c5d2e00c3d4aa108e821a056da4dc333d446db442417eae5'
)
DOWNLOAD_TIMEOUT_SECONDS = 30
REQUIRED_FILES = (
    'build_manifest.json',
    'ref_dbd.fasta',
    'ref_index.tsv',
    'ref_proteins.fasta',
    'motif_clusters.tsv',
    f'pfam/{RUNTIME_PFAM_NAME}',
    f'pfam/{RUNTIME_PFAM_MANIFEST}',
)


@dataclass(frozen=True)
class ReferenceArchiveSummary:
    archive: Path
    checksum_file: Path
    sha256: str
    size_bytes: int
    n_dbd_rows: int
    n_txt_motifs: int
    n_meme_motifs: int
    n_pfam_profiles: int


@dataclass(frozen=True)
class ReferenceInstallSummary:
    refs_dir: Path
    archive_sha256: str
    backup_dir: Path | None
    n_dbd_rows: int
    n_txt_motifs: int
    n_meme_motifs: int
    n_pfam_profiles: int


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _line_count(path, *, header=True):
    with Path(path).open() as handle:
        count = sum(1 for _ in handle)
    return max(0, count - int(header))


def validate_runtime_store(refs_dir):
    refs_dir = Path(refs_dir)
    missing = [name for name in REQUIRED_FILES if not (refs_dir / name).is_file()]
    txt_dir = refs_dir / 'motif_store' / 'txt'
    meme_dir = refs_dir / 'motif_store' / 'meme'
    if not txt_dir.is_dir():
        missing.append('motif_store/txt/')
    if not meme_dir.is_dir():
        missing.append('motif_store/meme/')
    if missing:
        raise FileNotFoundError(
            f"Incomplete PWM reference store at {refs_dir}: missing {', '.join(missing)}"
        )

    n_txt = sum(1 for _ in txt_dir.glob('*.txt'))
    n_meme = sum(1 for _ in meme_dir.glob('*.meme'))
    if not n_txt or not n_meme:
        raise ValueError(f'PWM reference motif store is empty: {refs_dir / "motif_store"}')
    pfam = validate_runtime_pfam(*runtime_pfam_paths(refs_dir))
    return {
        'n_dbd_rows': _line_count(refs_dir / 'ref_index.tsv'),
        'n_txt_motifs': n_txt,
        'n_meme_motifs': n_meme,
        'n_cluster_rows': _line_count(refs_dir / 'motif_clusters.tsv'),
        'n_pfam_profiles': pfam.n_profiles,
        'pfam_sha256': pfam.sha256,
    }


def _release_manifest(refs_dir, counts):
    refs_dir = Path(refs_dir)
    build_manifest = json.loads((refs_dir / 'build_manifest.json').read_text())
    return {
        'schema_version': 1,
        'artifact': 'boltzscan-pwm-reference-runtime',
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'boltzscan_version': __version__,
        'contents': (
            'Runtime DBD index, representative motif map, and scan-ready PWM files. '
            'The compact plant-TF Pfam HMM required at runtime is included. Raw '
            'download caches, logs, Pfam domtblout, and Tomtom work files are excluded.'
        ),
        'counts': counts,
        'build': build_manifest,
        'sources': {
            'cisbp': 'https://cisbp.ccbr.utoronto.ca/data/3_10/DataFiles/Bulk_downloads/EntireDataset',
            'jaspar': 'https://jaspar.elixir.no/api/v1/matrix/',
            'uniprot': 'https://rest.uniprot.org/',
            'pfam': 'https://ftp.ebi.ac.uk/pub/databases/Pfam/',
        },
        'source_notice': (
            'JASPAR states CC BY 4.0. Confirm CIS-BP redistribution terms before '
            'publishing this combined artifact publicly. Pfam data are CC0.'
        ),
        'core_sha256': {
            name: _sha256(refs_dir / name)
            for name in REQUIRED_FILES
        },
    }


def _normalise_tar_info(info):
    info.uid = 0
    info.gid = 0
    info.uname = ''
    info.gname = ''
    return info


def pack_reference_store(refs_dir, archive):
    """Write a runtime-only ``tar.gz`` and adjacent SHA256 file."""
    refs_dir = Path(refs_dir)
    archive = Path(archive)
    counts = validate_runtime_store(refs_dir)
    release_manifest = refs_dir / 'release_manifest.json'
    release_manifest.write_text(
        json.dumps(_release_manifest(refs_dir, counts), indent=2) + '\n'
    )

    archive.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=archive.parent,
        prefix=f'.{archive.name}.',
        suffix='.tmp',
        delete=False,
    ) as handle:
        temporary_archive = Path(handle.name)
    try:
        with tarfile.open(temporary_archive, 'w:gz', compresslevel=6) as tar:
            for name in (*REQUIRED_FILES, 'release_manifest.json'):
                tar.add(
                    refs_dir / name,
                    arcname=f'{ARCHIVE_ROOT}/{name}',
                    filter=_normalise_tar_info,
                )
            tar.add(
                refs_dir / 'motif_store',
                arcname=f'{ARCHIVE_ROOT}/motif_store',
                filter=_normalise_tar_info,
            )
        temporary_archive.replace(archive)
    except Exception:
        temporary_archive.unlink(missing_ok=True)
        raise

    archive.chmod(0o644)
    digest = _sha256(archive)
    checksum_file = Path(str(archive) + '.sha256')
    checksum_file.write_text(f'{digest}  {archive.name}\n')
    checksum_file.chmod(0o644)
    return ReferenceArchiveSummary(
        archive=archive,
        checksum_file=checksum_file,
        sha256=digest,
        size_bytes=archive.stat().st_size,
        **{key: counts[key] for key in ('n_dbd_rows', 'n_txt_motifs', 'n_meme_motifs')},
        n_pfam_profiles=counts['n_pfam_profiles'],
    )


def _download_url(url):
    """Convert a Google Drive share page to its unauthenticated download URL."""
    parsed = urlparse(str(url))
    if parsed.hostname not in {'drive.google.com', 'www.drive.google.com'}:
        return str(url)

    parts = [part for part in parsed.path.split('/') if part]
    file_id = None
    if len(parts) >= 3 and parts[:2] == ['file', 'd']:
        file_id = parts[2]
    elif parts and parts[0] in {'open', 'uc'}:
        file_id = parse_qs(parsed.query).get('id', [None])[0]
    if not file_id:
        return str(url)
    return (
        'https://drive.usercontent.google.com/download'
        f'?id={quote(file_id, safe="")}&export=download&confirm=t'
    )


def _format_size(size_bytes):
    size_bytes = int(size_bytes)
    if size_bytes >= 1024 * 1024:
        return f'{size_bytes / 1024 / 1024:.1f} MiB'
    if size_bytes >= 1024:
        return f'{size_bytes / 1024:.1f} KiB'
    return f'{size_bytes} bytes'


def _download(url, output, *, progress=None):
    report = progress or (lambda message: None)
    parsed = urlparse(str(url))
    if parsed.scheme == '' and Path(url).is_file():
        report(f'Using local PWM reference archive: {url}')
        shutil.copyfile(url, output)
        report(f'Archive copied: {_format_size(Path(output).stat().st_size)}')
        return
    download_url = _download_url(url)
    source = (
        'Google Drive'
        if urlparse(download_url).hostname == 'drive.usercontent.google.com'
        else urlparse(download_url).hostname or 'download server'
    )
    report(
        f'Connecting to {source} '
        f'(timeout: {DOWNLOAD_TIMEOUT_SECONDS} seconds)...'
    )
    request = Request(
        download_url,
        headers={'User-Agent': f'BoltzScan/{__version__}'},
    )
    try:
        with urlopen(
            request,
            timeout=DOWNLOAD_TIMEOUT_SECONDS,
        ) as response, Path(output).open('wb') as handle:
            try:
                total_bytes = int(response.headers.get('Content-Length', 0))
            except (TypeError, ValueError):
                total_bytes = 0
            if total_bytes:
                report(f'Download started: {_format_size(total_bytes)}')
            else:
                report('Download started: server did not report the file size')

            downloaded = 0
            next_percent = 10
            next_bytes = 10 * 1024 * 1024
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                downloaded += len(chunk)
                if total_bytes:
                    percent = min(100, downloaded * 100 // total_bytes)
                    while percent >= next_percent:
                        report(
                            f'Download: {next_percent:3d}% '
                            f'({_format_size(downloaded)} / {_format_size(total_bytes)})'
                        )
                        next_percent += 10
                elif downloaded >= next_bytes:
                    report(f'Downloaded: {_format_size(downloaded)}')
                    next_bytes += 10 * 1024 * 1024
            report(f'Download complete: {_format_size(downloaded)}')
    except OSError as exc:
        raise RuntimeError(
            f'Cannot download the PWM reference release from {source}: {exc}. '
            'The remote host may be blocked or unavailable. Download the tar.gz in a '
            'browser, then pass its local path with `install-pwm-refs --url FILE`.'
        ) from exc


def _safe_extract_runtime_archive(archive, destination):
    destination = Path(destination)
    with tarfile.open(archive, 'r:gz') as tar:
        members = tar.getmembers()
        for member in members:
            parts = PurePosixPath(member.name).parts
            if not parts or parts[0] != ARCHIVE_ROOT or '..' in parts:
                raise ValueError(f'Unsafe or unexpected archive member: {member.name}')
            if not (member.isdir() or member.isfile()):
                raise ValueError(f'Archive links/devices are not allowed: {member.name}')
            target = destination.joinpath(*parts)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            source = tar.extractfile(member)
            if source is None:
                raise ValueError(f'Cannot read archive member: {member.name}')
            with source, target.open('wb') as output:
                shutil.copyfileobj(source, output)


def install_reference_store(url, refs_dir, sha256, *, replace=False, progress=None):
    """Download, verify, safely extract, and atomically install a release."""
    report = progress or (lambda message: None)
    expected = str(sha256).strip().split()[0].lower()
    if len(expected) != 64 or any(char not in '0123456789abcdef' for char in expected):
        raise ValueError('--sha256 must be a 64-character hexadecimal digest')

    refs_dir = Path(refs_dir)
    refs_dir.parent.mkdir(parents=True, exist_ok=True)
    if refs_dir.exists() and any(refs_dir.iterdir()) and not replace:
        raise FileExistsError(
            f'Reference directory is not empty: {refs_dir}. Use --replace to move it '
            'to a timestamped backup before installing.'
        )
    if refs_dir.exists() and any(refs_dir.iterdir()):
        report('Existing reference store remains in place until the new archive is verified.')

    with tempfile.TemporaryDirectory(dir=refs_dir.parent) as temporary:
        temporary = Path(temporary)
        downloaded = temporary / 'pwm_refs.tar.gz'
        _download(url, downloaded, progress=report)
        report('Verifying SHA256...')
        observed = _sha256(downloaded)
        if observed != expected:
            raise ValueError(
                f'PWM reference archive SHA256 mismatch: expected {expected}, got {observed}'
            )
        report('SHA256 verified.')
        extracted = temporary / 'extracted'
        report('Extracting PWM reference archive...')
        _safe_extract_runtime_archive(downloaded, extracted)
        staged_refs = extracted / ARCHIVE_ROOT
        if not (staged_refs / 'release_manifest.json').is_file():
            raise ValueError('PWM reference archive is missing release_manifest.json')
        staged_hmm, staged_manifest = runtime_pfam_paths(staged_refs)
        if not staged_hmm.exists() and not staged_manifest.exists():
            # Compatibility bridge for the pre-Pfam built-in release. New
            # archives always carry these files and are validated below.
            report('Adding the packaged compact plant-TF Pfam library...')
            install_packaged_runtime_pfam(staged_refs)
        elif not staged_hmm.is_file() or not staged_manifest.is_file():
            raise ValueError(
                'PWM reference archive contains an incomplete compact Pfam library'
            )
        report('Validating extracted reference store...')
        counts = validate_runtime_store(staged_refs)
        report('Reference store validated; installing atomically...')

        backup = None
        if refs_dir.exists():
            if any(refs_dir.iterdir()):
                stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                backup = refs_dir.with_name(f'{refs_dir.name}.backup_{stamp}')
                if backup.exists():
                    raise FileExistsError(f'Refusing to overwrite backup: {backup}')
                refs_dir.rename(backup)
            else:
                refs_dir.rmdir()
        try:
            staged_refs.rename(refs_dir)
        except Exception:
            if backup is not None and not refs_dir.exists():
                backup.rename(refs_dir)
            raise

    return ReferenceInstallSummary(
        refs_dir=refs_dir,
        archive_sha256=observed,
        backup_dir=backup,
        **{key: counts[key] for key in ('n_dbd_rows', 'n_txt_motifs', 'n_meme_motifs')},
        n_pfam_profiles=counts['n_pfam_profiles'],
    )
