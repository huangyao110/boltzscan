# Repository Guidelines

## Project Structure & Module Organization

BoltzScan is an installable Python package for GRN inference workflows. Core source lives in `boltzscan/`: `cli.py` defines the console entry point, `fimocistarget/` handles promoter extraction and FIMO/cisTarget-style motif scoring, `cropd2p/` prepares protein-DNA Boltz inputs, and `utils/` contains shared path, IPSAE, and FIMO-to-Boltz helpers. `bscan.py` is a compatibility shim for `python bscan.py <cmd>`. Notebooks such as `boltzscan.ipynb` and `tair_tf2pwm.ipynb` are exploratory workflows. `png/` stores README assets. Large local inputs and run outputs belong under ignored paths such as `data/`, `tasks/`, and `boltz_results_*/`.

## Build, Test, and Development Commands

- `pip install -e .`: install the package and expose the `boltzscan` CLI.
- `pip install -e .[cistarg]`: install optional FIMO/cisTarget dependencies.
- `boltzscan --help` or `python -m boltzscan --help`: inspect available subcommands.
- `python bscan.py promoter -gff <genes.gff3> -g <genome.fa> -o promoters.fasta`: extract promoter FASTA.
- `python bscan.py cistarg -f <promoters.fa> -m <meme_dir> -t <tf2pwms.json> -o <out>`: build motif score outputs.
- `python -m build`: build distributions when the `build` package is installed.

Use Python 3.11 for local work; `pyproject.toml` allows Python `>=3.10,<3.13`.

## Coding Style & Naming Conventions

Follow the existing Python style: 4-space indentation, `snake_case` functions and variables, concise module-level helpers, and `Path` for filesystem paths where practical. Keep CLI wiring in `boltzscan/cli.py`; place substantive workflow logic in the relevant package module. Do not normalize existing mixed English/Chinese comments or messages unless touching the surrounding logic for a reason.

## Testing Guidelines

There is no repository-level test suite or lint config. Before submitting changes, run the specific CLI path you changed with small fixture inputs when possible, and at minimum run `python -m compileall boltzscan` to catch syntax/import issues. If adding tests, use `tests/` with `test_<module>.py` naming and document any new `pytest` command in this file or the PR.

## Commit & Pull Request Guidelines

Recent history uses short imperative or summary-style subjects, for example `Refactor into installable Python package`. Keep commit subjects concise and focused. PRs should describe the workflow affected, list commands or notebooks used for validation, link related issues when available, and include screenshots only for documentation or visual-output changes.

## Security & Configuration Tips

Do not commit large datasets, generated Boltz runs, archives, logs, or machine-specific outputs. Treat vendored external trees such as `boltzscan/boltz/` and `boltzscan/cropd2p/disopred/` as separate projects with their own dependencies and conventions.
