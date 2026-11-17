"""Shared fixtures: synthetic GenBank files written with Biopython."""

from collections.abc import Callable, Sequence
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

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
