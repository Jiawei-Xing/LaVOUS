# Repository Guidelines

## Project Structure & Module Organization

This repository contains the LAVOUS/SingleCellStochastics Python package for lineage-aware BM/OU models on single-cell RNA-seq counts. Core package code lives in `singlecellstochastics/`; important modules include `preprocess.py`, `likelihood.py`, `elbo.py`, `optimize.py`, `plasticity.py`, `ou_diff.py`, `calibrate.py`, and `reconstruct.py`. Example inputs and expected outputs are under `examples/`. Developer notes are in `docs/source_map.md`. Larger exploratory or benchmark material is kept in `expression_simulation/`, `lineage_simulation/`, `real_data/`, and `test/`.

## Build, Test, and Development Commands

Run commands from the repository root.

```bash
pip install -e .
python -m compileall singlecellstochastics
python - <<'PY'
import singlecellstochastics
print(singlecellstochastics.__version__)
PY
```

`pip install -e .` installs the package and console scripts locally. `compileall` catches syntax errors across the package. The import check verifies that the editable install resolves correctly. CLI entry points include `lavous-heritability`, `lavous-diff`, `lavous-calibrate`, `lavous-reconstruct`, and `lavous-simulate`.

## Coding Style & Naming Conventions

Use Python with 4-space indentation and keep imports grouped as standard library, third-party, then local package imports. Follow the existing module style: lowercase module names, snake_case functions and variables, and concise docstrings for public workflows or numerical routines. `black` is listed as a project dependency; use it for formatting touched Python files when practical. Avoid committing generated caches such as `__pycache__/`, `.ipynb_checkpoints/`, and local `wandb/` runs.

## Testing Guidelines

There is no formal pytest suite configured yet. For code changes, at minimum run `python -m compileall singlecellstochastics` and an import check. For workflow changes, run the relevant CLI against `examples/input_data/` and write temporary output outside tracked example results, for example under `/tmp/lavous-check/`. Name future automated tests `test_*.py` and place them in a dedicated test directory rather than mixing them with large benchmark artifacts.

## Commit & Pull Request Guidelines

Recent history uses short, imperative, lower-case commit messages, such as `fix optimization nan` or `update empirical calibrations`; keep new commits similarly focused. Pull requests should describe the workflow affected, list validation commands run, and call out changes to output schemas, CLI options, dependencies, or statistical assumptions. Include links to related issues and small example outputs when behavior changes.
