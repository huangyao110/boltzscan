"""Keep the known RDKit CCD pickle-version warning out of experiment logs."""

try:
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.warning")
except ImportError:
    pass
