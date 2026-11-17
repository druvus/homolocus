# homolocus

Find homologous genomic regions in bacterial genomes from the proteins of a reference region.
The authoritative design is `SPECIFICATION.md`. Read the relevant section before implementing anything.
The package lives in `src/homolocus/`. Open source, MIT licensed.

## Commands

- Tests: `pytest -m "not diamond"` (fast). Full: `pytest` (needs the `diamond` executable).
- Lint/format: `ruff check --fix` and `ruff format`
- Types: `mypy src`

## Workflow

- Test first. Write a failing test for the behaviour, then implement it. Never mark work done with failing tests.
- One milestone per session. Plan first, wait for approval, then implement.
- Small commits with a message that names the spec section implemented (e.g. "Add GeneMatch aggregation (spec 5.7)").

## Spec discipline

- Follow tags in SPECIFICATION.md: [DECIDED] is binding, [PROPOSED] is implemented as written but kept easy to change, [OPEN] is never resolved silently. For an [OPEN] item, suggest a default and ask me.
- Do not add features listed as out of scope for v0.1.
- If code and spec disagree, stop and tell me instead of choosing.

## Conventions

- Internal coordinates are 0-based half-open. User-facing coordinates (CLI, TSV, plots) are 1-based inclusive. Convert only in `coordinates.py` and the I/O layer.
- Gene ID format: `{contig_id}|{start_1based}-{end_1based}|{+|-}`. Never key on `locus_tag`.
- Detection, matching and scoring code is pure: no file I/O, no subprocess, no plotting. DIAMOND is isolated in `search/diamond.py`. Plotting only consumes result models.
- Keep raw evidence: individual HSPs are never discarded when aggregated.
- Python 3.10+, full type hints, dataclasses for models, small functions with a single responsibility, docstrings on public API.
- No broad `except`. Warnings go through the `warnings`/`logging` modules, never `print`.
- Tests use small synthetic genomes and proteins and are deterministic. Tests that need DIAMOND carry `@pytest.mark.diamond`.

## Do not

- Do not edit `SPECIFICATION.md` without asking.
- Do not add dependencies beyond those in the spec without asking.
