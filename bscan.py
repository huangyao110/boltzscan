"""Compat shim: keeps `python bscan.py <cmd>` working after the package refactor."""
from boltzscan.cli import main

if __name__ == '__main__':
    main()
