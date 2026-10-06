"""Obtain and check the raw inputs of the case study.

Every input the analysis reads is deposited as one archive whose layout mirrors
``data/raw``.  :func:`fetch_inputs` downloads it and unpacks it into
:data:`analysis.paths.RAW`; :func:`verify_inputs` checks each required file
against the SHA-256 checksums in ``data/raw/MANIFEST.sha256``.
"""
from __future__ import annotations

import hashlib
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path

from . import paths

# The deposit is not yet public; both are set once the DOI is assigned.
DEPOSIT_DOI: str | None = None
DEPOSIT_URL: str | None = None      # direct link to the deposit's .zip archive

MANIFEST = paths.RAW / "MANIFEST.sha256"

# What the analysis reads from data/raw, relative to it.  The GDC clinical
# tables are a directory read whole, so every file in it is required.
REQUIRED_FILES = (
    "delivered/sarcoma_plus_NciLandscape.tpm.N_33967.260709T0925.tsv",
    "delivered/sarcoma_plus_NciLandscape.pca_coords.log2_tpmPlusOne.N_33967.260705T1319.tsv",
    "gencode_v36_gene_map.tsv",
    "brohl_supplement/NIHMS1759318-supplement-2.xlsx",   # Brohl Table S1
    "brohl_supplement/NIHMS1759318-supplement-3.xlsx",   # Brohl Table S2 (TCR)
)
REQUIRED_DIRS = ("gdc_clinical",)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def required_inputs() -> list[str]:
    """Paths, relative to ``data/raw``, of every input the analysis reads."""
    files = list(REQUIRED_FILES)
    for d in REQUIRED_DIRS:
        files += sorted(p.relative_to(paths.RAW).as_posix()
                        for p in (paths.RAW / d).glob("*") if p.is_file())
    return files


def fetch_inputs(force: bool = False) -> Path:
    """Download the deposited inputs and unpack them into ``data/raw``.

    Skips the download when every input is already present and matches the
    manifest, unless ``force`` is set.  Returns ``paths.RAW``.
    """
    if DEPOSIT_URL is None:
        raise RuntimeError(
            "The input deposit has no public URL yet (analysis.fetch.DEPOSIT_URL "
            "is unset).  Until it does, place the raw inputs under "
            f"{paths.RAW} by hand and check them with verify_inputs().")
    if not force and MANIFEST.exists():
        try:
            verify_inputs()
            print(f"inputs already present and verified in {paths.RAW}")
            return paths.RAW
        except (FileNotFoundError, ValueError):
            pass

    paths.RAW.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp) / "deposit.zip"
        print(f"downloading {DEPOSIT_URL}")
        with urllib.request.urlopen(DEPOSIT_URL) as r, open(archive, "wb") as fh:
            shutil.copyfileobj(r, fh)
        with zipfile.ZipFile(archive) as z:
            root = paths.RAW.resolve()
            for name in z.namelist():
                if not (root / name).resolve().is_relative_to(root):
                    raise ValueError(f"archive member outside data/raw: {name}")
            z.extractall(paths.RAW)
    print(f"unpacked into {paths.RAW}")
    verify_inputs()
    return paths.RAW


def _read_manifest() -> dict[str, str]:
    if not MANIFEST.exists():
        raise FileNotFoundError(f"no manifest at {MANIFEST}")
    out = {}
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            out[name.lstrip("*")] = digest
    return out


def verify_inputs() -> None:
    """Check every required input exists and matches ``MANIFEST.sha256``.

    Raises ``FileNotFoundError`` naming the missing files, or ``ValueError``
    naming the files whose checksum differs from the manifest.
    """
    manifest = _read_manifest()
    wanted = sorted(set(required_inputs()) | set(manifest))
    missing = [f for f in wanted if not (paths.RAW / f).is_file()]
    if missing:
        raise FileNotFoundError(
            f"missing from {paths.RAW}:\n  " + "\n  ".join(missing)
            + "\nRun analysis.fetch.fetch_inputs() to obtain them.")
    unlisted = [f for f in wanted if f not in manifest]
    mismatched = [f for f in wanted
                  if f in manifest and _sha256(paths.RAW / f) != manifest[f]]
    if unlisted or mismatched:
        msg = []
        if mismatched:
            msg.append("checksum differs from MANIFEST.sha256:\n  "
                       + "\n  ".join(mismatched))
        if unlisted:
            msg.append("required but not in MANIFEST.sha256:\n  "
                       + "\n  ".join(unlisted))
        raise ValueError("\n".join(msg))
    print(f"{len(wanted)} input files present and verified in {paths.RAW}")


def write_manifest() -> Path:
    """Regenerate ``data/raw/MANIFEST.sha256`` from the required inputs.

    For the authors, after a deliberate change to the deposited inputs.  The
    format is that of ``sha256sum``, so ``sha256sum -c MANIFEST.sha256`` run
    inside ``data/raw`` checks it too.
    """
    files = sorted(required_inputs())
    missing = [f for f in files if not (paths.RAW / f).is_file()]
    if missing:
        raise FileNotFoundError("cannot write the manifest, missing:\n  "
                                + "\n  ".join(missing))
    lines = [f"{_sha256(paths.RAW / f)}  {f}" for f in files]
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {MANIFEST} ({len(files)} files)")
    return MANIFEST


def run_all() -> None:
    """Make sure the inputs are in place: fetch them if possible, then verify."""
    if DEPOSIT_URL is None:
        verify_inputs()
    else:
        fetch_inputs()
