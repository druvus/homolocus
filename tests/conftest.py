"""Shared fixtures: synthetic GenBank files written with Biopython."""

import os
import random
import shutil
from collections.abc import Callable, Sequence
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from genome_files import genbank
from homolocus import RegionFinder
from sample_world import results as sample_results

MISSING_DIAMOND = "DIAMOND executable not found on PATH; install DIAMOND 2.0+"


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Handle DIAMOND tests when DIAMOND is not installed (spec 15).

    They are skipped with a clear message, unless HOMOLOCUS_REQUIRE_DIAMOND=1
    (set in CI), which makes a missing DIAMOND an error.
    """
    needs_diamond = [item for item in items if "diamond" in item.keywords]
    if not needs_diamond or shutil.which("diamond") is not None:
        return
    if os.environ.get("HOMOLOCUS_REQUIRE_DIAMOND") == "1":
        raise pytest.UsageError(f"HOMOLOCUS_REQUIRE_DIAMOND=1 but {MISSING_DIAMOND}")
    for item in needs_diamond:
        item.add_marker(pytest.mark.skip(reason=MISSING_DIAMOND))


@pytest.fixture
def searches(monkeypatch: pytest.MonkeyPatch) -> list[bool]:
    """Stand-in for DIAMOND in RegionFinder that records each search and finds nothing.

    Each entry is True for a self-search and False for a main search.
    """
    calls: list[bool] = []

    def fake_search(index, config, path, *, self_search=False):
        path.write_text("")
        calls.append(self_search)
        return ()

    monkeypatch.setattr("homolocus.finder.search", fake_search)
    monkeypatch.setattr("homolocus.finder.diamond_version", lambda: "2.2.9")
    return calls


@pytest.fixture
def inputs(tmp_path, dna, make_record):
    """Command-line arguments for a readable GenBank reference and target."""
    reference = genbank(tmp_path / "ref.gbk", "ref", dna, make_record)
    target = genbank(tmp_path / "target.gbk", "t1", dna, make_record)
    output = tmp_path / "out"
    return ["find", "--reference", str(reference), "--genome", str(target), "--output", str(output)]


@pytest.fixture
def finder_calls(monkeypatch):
    """Replace the search with the sample results and record each RegionFinder used."""
    calls = []

    def fake_find(finder, genome):
        calls.append((finder, genome))
        return sample_results()

    monkeypatch.setattr(RegionFinder, "find", fake_find)
    return calls


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
