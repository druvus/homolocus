"""Packaging: the type-hint marker and the console script."""

from importlib.metadata import entry_points
from importlib.resources import files

from homolocus import cli


def test_package_ships_py_typed_marker() -> None:
    assert files("homolocus").joinpath("py.typed").is_file()


def test_homolocus_command_runs_the_cli_main() -> None:
    (script,) = entry_points(group="console_scripts", name="homolocus")
    assert script.load() is cli.main
