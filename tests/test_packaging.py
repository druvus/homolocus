"""Packaging: the package declares its inline type hints (PEP 561)."""

from importlib.resources import files


def test_package_ships_py_typed_marker() -> None:
    assert files("homolocus").joinpath("py.typed").is_file()
