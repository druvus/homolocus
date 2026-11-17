"""Shared fixtures: synthetic GenBank files written with Biopython."""

import random
import shutil
from collections.abc import Callable, Sequence
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Skip DIAMOND tests with a clear message when DIAMOND is not installed (spec 15)."""
    if shutil.which("diamond") is not None:
        return
    skip = pytest.mark.skip(reason="DIAMOND executable not found on PATH; install DIAMOND 2.0+")
    for item in items:
        if "diamond" in item.keywords:
            item.add_marker(skip)


RecordFactory = Callable[[str, str, Sequence[SeqFeature]], SeqRecord]
GenBankWriter = Callable[[Sequence[SeqRecord]], Path]


@pytest.fixture
def make_record() -> RecordFactory:
    """Build a DNA SeqRecord with the given ID, sequence and features."""

    def factory(record_id: str, sequence: str, features: Sequence[SeqFeature]) -> SeqRecord:
        record = SeqRecord(Seq(sequence), id=record_id, name=record_id, description="synthetic")
        record.annotations["molecule_type"] = "DNA"
        record.features = list(features)
        return record

    return factory


@pytest.fixture
def write_genbank(tmp_path: Path) -> GenBankWriter:
    """Write records to a GenBank file in a temporary directory and return its path."""

    def writer(records: Sequence[SeqRecord]) -> Path:
        path = tmp_path / "input.gbk"
        SeqIO.write(list(records), path, "genbank")
        return path

    return writer


_SENSE_CODONS = tuple(
    a + b + c
    for a in "ACGT"
    for b in "ACGT"
    for c in "ACGT"
    if a + b + c not in {"TAA", "TAG", "TGA"}
)


class SyntheticDna:
    """Deterministic random DNA with planted open reading frames."""

    def __init__(self, seed: int = 42) -> None:
        self._rng = random.Random(seed)

    def random(self, length: int) -> str:
        """Uniform random sequence."""
        return "".join(self._rng.choice("ACGT") for _ in range(length))

    def orf(self, codons: int) -> str:
        """ATG, ``codons`` random sense codons, then TAA."""
        return "ATG" + "".join(self._rng.choice(_SENSE_CODONS) for _ in range(codons)) + "TAA"


@pytest.fixture
def dna() -> SyntheticDna:
    """A seeded synthetic-DNA generator."""
    return SyntheticDna()


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


class SyntheticProteins:
    """Deterministic random proteins and mutated copies."""

    def __init__(self, seed: int = 4) -> None:
        self._rng = random.Random(seed)

    def random(self, length: int) -> str:
        """Uniform random protein sequence."""
        return "".join(self._rng.choice(AMINO_ACIDS) for _ in range(length))

    def mutate(self, sequence: str, identity: float) -> str:
        """Copy of ``sequence`` where each residue is kept with probability ``identity``."""
        return "".join(
            residue if self._rng.random() < identity else self._rng.choice(AMINO_ACIDS)
            for residue in sequence
        )


@pytest.fixture
def proteins() -> SyntheticProteins:
    """A seeded synthetic-protein generator."""
    return SyntheticProteins()
