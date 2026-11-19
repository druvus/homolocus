"""Genome presentation names and order from a mapping file (spec 5.3, 12, 13).

A mapping file has one genome per line: the genome's file name, a tab and its
presentation name. Blank lines and lines starting with ``#`` are ignored.
Reading it is the I/O step; applying it to genomes is pure.
"""

import warnings
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path

from homolocus.models.genome import Genome


def read_genome_names(path: str | Path) -> dict[str, str]:
    """File name to presentation name, in file order.

    Raises ValueError naming the line for a malformed line, or a file name or
    presentation name listed twice.
    """
    names: dict[str, str] = {}
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    for number, line in enumerate(lines, start=1):
        if line.strip() and not line.startswith("#"):
            file, name = _fields(line, number)
            _check_new(names, file, name, number)
            names[file] = name
    return names


def apply_genome_names(genomes: Sequence[Genome], names: Mapping[str, str]) -> list[Genome]:
    """Genomes in mapping order with their presentation names; unlisted ones after.

    A genome matches a line by its source file name, ignoring the path. Lines
    that match no genome are reported with a warning.
    """
    by_file = {_file_name(genome): genome for genome in genomes}
    for file in names:
        if file not in by_file:
            warnings.warn(f"names file lists {file}, which is not a --genome", stacklevel=2)
    listed = [replace(by_file[f], display_name=name) for f, name in names.items() if f in by_file]
    unlisted = [genome for genome in genomes if _file_name(genome) not in names]
    return [*listed, *unlisted]


def _fields(line: str, number: int) -> tuple[str, str]:
    fields = line.split("\t")
    if len(fields) != 2 or not all(field.strip() for field in fields):
        raise ValueError(f"names file line {number}: expected a file name, a tab and a name")
    return fields[0].strip(), fields[1].strip()


def _check_new(names: Mapping[str, str], file: str, name: str, number: int) -> None:
    if file in names:
        raise ValueError(f"names file line {number}: file {file} is listed twice")
    if name in names.values():
        raise ValueError(f"names file line {number}: name {name!r} is used twice")


def _file_name(genome: Genome) -> str:
    return Path(genome.source or "").name
