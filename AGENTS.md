# Repository Guidelines

## Project Structure & Module Organization

This repository contains a Python package for annotating VCF files. Core source code lives in `vcf_annotator/`: `cli.py` handles command-line parsing, `processor.py` coordinates record processing, `transcripts.py` and `chromosome.py` manage reference data, and annotator implementations live under `vcf_annotator/annotators/`. Tests are in `tests/`, with reusable fixtures and small VCF/reference files in `tests/data/`. Behaviour-driven splice-distance scenarios are stored in `tests/features/`.

## Build, Test, and Development Commands

- `uv sync --extra test`: create or update the local environment with runtime and test dependencies from `pyproject.toml` and `uv.lock`.
- `uv run pytest`: run the full test suite configured under `tests/`.
- `uv run pytest tests/test_cli.py`: run CLI integration tests only.
- `uv run python -m vcf_annotator.cli --input tests/data/input.vcf --output /tmp/annotated.vcf`: exercise the CLI locally.
- `docker build -t vcf-annotator .`: build the container image when validating Docker changes.

## Coding Style & Naming Conventions

Use Python 3.10+ syntax and keep modules typed where practical. Follow PEP 8 with 4-space indentation, `snake_case` functions and variables, `PascalCase` classes, and uppercase constants such as `DATA_DIR`. Keep annotator-specific logic inside `vcf_annotator/annotators/` and prefer small helper functions over duplicating parsing or formatting logic. Preserve existing `from __future__ import annotations` usage in new package modules.

## Testing Guidelines

The project uses `pytest` and `pytest-bdd`. Name test files `test_*.py`, matching the pytest configuration in `pyproject.toml`. Add focused unit tests for helpers and annotators, and CLI-level tests when changing argument parsing, VCF output, or TSV output. Keep test fixtures small and deterministic in `tests/data/`; document any fixture regeneration steps in `README.md`.

## Commit & Pull Request Guidelines

Recent history uses concise conventional-style prefixes such as `docs:`, `test:`, `build:`, and `chore:`; continue that pattern where it fits. Keep commits scoped to one logical change. Pull requests should describe the user-visible effect, list test commands run, and call out any changes to fixture data, CLI arguments, Docker behaviour, or annotation output formats.

## Security & Configuration Tips

Do not commit private genomic data or large raw reference downloads. Use the sampled fixtures in `tests/data/` for tests, and keep generated outputs in temporary paths such as `/tmp` unless they are intentional fixtures.
