"""README.md stays true to the package: flags, defaults, output files and examples."""

import contextlib
import dataclasses
import io
import re
import shlex
from pathlib import Path

import pytest

from genome_files import fasta, genbank
from homolocus import DiamondSearch, FamilyGrouping, PlotOptions, RegionDetection, ScoringWeights
from homolocus.cli import build_parser
from sample_world import results as sample_results

README = (Path(__file__).parent.parent / "README.md").read_text(encoding="utf-8")
USER_PART = README.split("## Development")[0]
CONFIGS = (
    (DiamondSearch, ""),
    (FamilyGrouping, ""),
    (RegionDetection, ""),
    (ScoringWeights, "weight-"),
    (PlotOptions, "plot-"),
)
TABLE_ROW = re.compile(r"^\| `(--[a-z-]+)` \| `([^`]*)` \|", re.MULTILINE)


def code_blocks(language: str) -> list[str]:
    return re.findall(rf"```{language}\n(.*?)```", README, re.DOTALL)


def find_options() -> set[str]:
    """Every option ``homolocus --help`` and ``homolocus find --help`` list."""
    help_text = io.StringIO()
    for argv in (["--help"], ["find", "--help"]):
        with contextlib.redirect_stdout(help_text), contextlib.suppress(SystemExit):
            build_parser().parse_args(argv)
    return set(re.findall(r"(?<![\w-])--[a-z][a-z-]*", help_text.getvalue()))


def config_defaults() -> dict[str, str]:
    return {
        f"--{prefix}{field.name.replace('_', '-')}": str(field.default)
        for config, prefix in CONFIGS
        for field in dataclasses.fields(config)
    }


def tree_paths(block: str) -> set[str]:
    """Files of a ``tree``-style block, as paths below its root directory."""
    paths, stack = set(), []
    for line in block.splitlines()[1:]:
        match = re.match(r"^((?:│   |    )*)[├└]── (\S+)", line)
        if not match:
            continue
        depth, name = len(match.group(1)) // 4, match.group(2)
        stack[depth:] = [name.rstrip("/")]
        if not name.endswith("/") and name != "...":
            paths.add("/".join(stack))
    return paths


def test_every_flag_in_the_user_docs_exists() -> None:
    named = set(re.findall(r"(?<![\w-])--[a-z][a-z-]*", USER_PART))
    assert named, "the README names no flags"
    assert named <= find_options()


def test_option_tables_cover_every_setting_with_its_default() -> None:
    rows = dict(TABLE_ROW.findall(README))
    defaults = config_defaults()
    assert set(defaults) <= set(rows)
    for flag, default in defaults.items():
        shown = rows[flag]
        assert shown == ("unset" if default == "None" else default), flag


def test_output_tree_matches_what_write_produces(tmp_path) -> None:
    (block,) = [b for b in code_blocks("text") if b.startswith("results/")]
    documented = tree_paths(block)
    sample_results().write(tmp_path)
    written = {str(p.relative_to(tmp_path)) for p in tmp_path.rglob("*") if p.is_file()}
    assert documented <= written
    first_only = {p for p in written if not re.search(r"candidate_00[2-9]", p)}
    assert first_only <= documented


def test_command_line_examples_parse() -> None:
    commands = [
        shlex.split(line.replace("\\\n", " "))
        for block in code_blocks("bash")
        for line in re.split(r"(?<!\\)\n", block)
        if line.startswith("homolocus find")
    ]
    assert len(commands) >= 2
    for argv in commands:
        build_parser().parse_args(argv[1:])


@pytest.mark.diamond
def test_python_quick_start_runs(tmp_path, dna, make_record, monkeypatch) -> None:
    (block,) = [b for b in code_blocks("python") if "RegionFinder(" in b]
    orfs = [dna.orf(150) for _ in range(3)]
    genbank(tmp_path / "island.gbk", "island", dna, make_record, orfs)
    fasta(tmp_path / "target.fna", "contig_1", dna, orfs)
    monkeypatch.chdir(tmp_path)
    namespace: dict[str, object] = {}
    exec(block, namespace)  # noqa: S102 - the README's own example
    results = namespace["results"]
    assert results.best_candidate is not None
    assert (tmp_path / "results" / "summary.tsv").is_file()
