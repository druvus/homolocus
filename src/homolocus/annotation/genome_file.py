"""Target genome files: format detection and annotation choice (spec 2.2).

GenBank with supplied CDS keeps its annotation under ``"auto"``. FASTA,
GenBank without CDS, and anything read with ``"pyrodigal"`` get genes
predicted by Pyrodigal.
"""

from pathlib import Path
from typing import Literal

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from homolocus.annotation.genbank import (
    check_unique_ids,
    genbank_genome,
    read_records,
    record_contig,
    record_id,
    record_sequence,
)
from homolocus.annotation.prodigal import predict_genes
from homolocus.models.genome import Genome

Annotation = Literal["auto", "pyrodigal"]
FileFormat = Literal["genbank", "fasta"]

ANNOTATION_MODES: tuple[Annotation, ...] = ("auto", "pyrodigal")


def detect_format(path: str | Path) -> FileFormat:
    """Detect GenBank (``LOCUS``) or FASTA (``>``) from the first non-empty line."""
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("LOCUS"):
                return "genbank"
            if stripped.startswith(">"):
                return "fasta"
            if stripped:
                raise ValueError(f"unrecognised format in {path}; expected GenBank or FASTA")
    raise ValueError(f"{path} is empty")


def read_genome(path: str | Path, annotation: Annotation = "auto") -> Genome:
    """Read a target genome from a GenBank or FASTA file (spec 2.2)."""
    if annotation not in ANNOTATION_MODES:
        raise ValueError(f"unknown annotation mode {annotation!r}; expected 'auto' or 'pyrodigal'")
    records = _read_fasta(path) if detect_format(path) == "fasta" else read_records(path)
    if annotation == "auto" and _has_cds(records):
        return genbank_genome(records, source=str(path))
    return _predicted_genome(records, source=str(path))


def _read_fasta(path: str | Path) -> list[SeqRecord]:
    with open(path, encoding="utf-8") as handle:
        return list(SeqIO.parse(handle, "fasta"))


def _has_cds(records: list[SeqRecord]) -> bool:
    return any(feature.type == "CDS" for record in records for feature in record.features)


def _predicted_genome(records: list[SeqRecord], source: str) -> Genome:
    check_unique_ids(records, source)
    predicted = predict_genes({record_id(r): record_sequence(r) for r in records})
    contigs = {record_id(r): record_contig(r, predicted[record_id(r)]) for r in records}
    return Genome(contigs=contigs, source=source)
