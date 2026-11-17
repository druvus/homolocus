# homolocus

Find homologous genomic regions in bacterial genomes from the proteins of a reference region.
The authoritative design is `SPECIFICATION.md`. Read the relevant section before implementing anything.
The package lives in `src/homolocus/`. Open source, MIT licensed.

## Commands

- Tests: `pytest -m "not diamond"` (fast). Full: `pytest` (needs the `diamond` executable).
- Lint/format: `ruff check --fix` and `ruff format`
- Types: `mypy src`
- CI: `.github/workflows/tests.yml` runs lint, types and the full tests (Python 3.10–3.13, DIAMOND pinned by version and SHA-256). Upgrading DIAMOND means changing both values there.
- Quality gates (spec 17). The tools are optional external checkouts next to this repo (`../crapper`, `../dryer`, `../mutator`):
  - CRAP below 10 per function: `../crapper/crapper src`
  - No duplication at 0.8 or above: `../dryer/dryer --threshold 0.8 src`
  - No surviving mutants in changed code: `../mutator/mutator --test-command "pytest -q -m 'not diamond'" <changed file>`

## Workflow

- Every code change follows the `ub-mode` skill (`.claude/skills/ub-mode/SKILL.md`): First Make it Work, Then Make it Right, Follow Practice. The overrides below take precedence.
- Test first. Write a failing test for the behaviour, see it fail for the expected reason, then implement it. Never mark work done with failing tests.
- One milestone per session. Plan first, wait for approval, then implement.
- Small commits. Each message names the spec section implemented (e.g. "Add GeneMatch aggregation (spec 5.7)") and says why, what changed, and how it was proved (tests and quality gates run).
- Keep a log of each action and the reason for it. Report the log at the end of the session in chat; do not commit it.

## ub-mode overrides

Where ub-mode and this file disagree, this file wins.

- Do not push to `main`. Push to the session's feature branch and open a pull request. The user merges.
- "Redesign from first principles" and "Attack the premise" never override a [DECIDED] spec item or the rule against editing `SPECIFICATION.md` without asking. Propose the change instead.
- Worktree cleanup never deletes uncommitted work.

## Spec discipline

- Follow tags in SPECIFICATION.md: [DECIDED] is binding, [PROPOSED] is implemented as written but kept easy to change, [OPEN] is never resolved silently. For an [OPEN] item, suggest a default and ask me.
- Do not add features listed as out of scope for v0.1.
- If code and spec disagree, stop and tell me instead of choosing.

## Conventions

- Internal coordinates are 0-based half-open. User-facing coordinates (CLI, TSV, plots) are 1-based inclusive. Convert only in `coordinates.py` and the I/O layer.
- Gene ID format: `{contig_id}|{start_1based}-{end_1based}|{+|-}`. Never key on `locus_tag`.
- Detection, matching and scoring code is pure: no file I/O, no subprocess, no plotting. DIAMOND is isolated in `search/diamond.py`. Plotting only consumes result models.
- Parse and validate outside data once, where it enters. Inside the package, trust it.
- Models make illegal states impossible to build (strand, intervals, coverage and identity ranges), enforced at construction.
- Keep raw evidence: individual HSPs are never discarded when aggregated.
- Use the glossary terms (spec 20) in names, tests and docs.
- Python 3.10+, full type hints, dataclasses for models, small functions with a single responsibility, docstrings on public API.
- Do not silence mypy: no blanket `# type: ignore`, no `cast` to hide an error.
- No broad `except`. Warnings go through the `warnings`/`logging` modules, never `print`.
- Tests check behaviour through public functions, use small synthetic genomes and proteins, and are deterministic. Tests that need DIAMOND carry `@pytest.mark.diamond`.

## Do not

- Do not edit `SPECIFICATION.md` without asking.
- Do not add dependencies beyond those in the spec without asking.
