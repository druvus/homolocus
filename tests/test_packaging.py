"""Packaging: the type-hint marker, the console script and dependency floors."""

import re
from importlib.metadata import entry_points, requires
from importlib.resources import files

from homolocus import cli


def test_package_ships_py_typed_marker() -> None:
    assert files("homolocus").joinpath("py.typed").is_file()


def test_homolocus_command_runs_the_cli_main() -> None:
    (script,) = entry_points(group="console_scripts", name="homolocus")
    assert script.load() is cli.main


def test_every_runtime_dependency_has_a_tested_floor() -> None:
    runtime = [r for r in requires("homolocus") or [] if "extra ==" not in r]
    assert {re.split(r"[ <>=!~;]", r)[0].lower() for r in runtime} == {
        "biopython",
        "pandas",
        "matplotlib",
        "pyrodigal",
    }
    for requirement in runtime:
        assert re.fullmatch(r"[A-Za-z0-9_.-]+>=\d+(\.\d+)*", requirement), requirement
