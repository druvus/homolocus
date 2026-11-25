"""Target genome files: format detection and annotation choice (spec 2.2).

GenBank with supplied CDS keeps its annotation under ``"auto"``. FASTA,
GenBank without CDS, and anything read with ``"pyrodigal"`` get genes
predicted by Pyrodigal.
"""

import logging
from pathlib import Path
from typing import Literal

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from homolocus.annotation.genbank import (
    check_unique_ids,
    genbank_genome,
    genome_name,
    open_text,
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

logger = logging.getLogger(__name__)


def detect_format(path: str | Path) -> FileFormat:
    """Detect GenBank (``LOCUS``) or FASTA (``>``) from the first non-empty line."""
    with open_text(path) as handle:
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
    """Read a target genome from a GenBank or FASTA file (spec 2.2).

    The genome is named after the file, without its extension (spec 5.3).
    """
    if annotation not in ANNOTATION_MODES:
        raise ValueError(f"unknown annotation mode {annotation!r}; expected 'auto' or 'pyrodigal'")
    records = _read_fasta(path) if detect_format(path) == "fasta" else read_records(path)
    if annotation == "auto" and _has_cds(records):
        genome = genbank_genome(records, source=str(path))
    else:
        genome = _predicted_genome(records, source=str(path))
    n_genes = sum(contig.n_genes for contig in genome.contigs.values())
    logger.info("read genome %s, contigs: %d, genes: %d", genome.name, len(genome.contigs), n_genes)
    return genome


def _read_fasta(path: str | Path) -> list[SeqRecord]:
    with open_text(path) as handle:
        return list(SeqIO.parse(handle, "fasta"))


def _has_cds(records: list[SeqRecord]) -> bool:
    return any(feature.type == "CDS" for record in records for feature in record.features)


def _predicted_genome(records: list[SeqRecord], source: str) -> Genome:
    check_unique_ids(records, source)
    predicted = predict_genes({record_id(r): record_sequence(r) for r in records})
    contigs = {record_id(r): record_contig(r, predicted[record_id(r)]) for r in records}
    return Genome(contigs=contigs, name=genome_name(source), source=source)
