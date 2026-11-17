"""GenBank input (spec 2.1, 2.2, 3). The I/O boundary for annotated genomes.

Features are checked and translated here, once. Everything built from this
module is trusted by the rest of the package.
"""

import warnings
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    SeqFeature,
    SimpleLocation,
)
from Bio.SeqRecord import SeqRecord

from homolocus.annotation.translation import (
    InternalStop,
    InvalidSupplied,
    Protein,
    Translation,
    Untranslatable,
    check_supplied_translation,
    translate_cds,
)
from homolocus.coordinates import to_zero_based
from homolocus.models.contig import Contig
from homolocus.models.gene import Gene, as_strand
from homolocus.models.genome import Genome
from homolocus.models.reference import ReferenceRegion
from homolocus.models.skipped_feature import SkippedFeature, SkipReason

Region = tuple[int, int]
Location = SimpleLocation | CompoundLocation


class AnnotationWarning(UserWarning):
    """Summary of CDS features that were skipped, untranslatable, trimmed or unparsed."""


@dataclass(frozen=True)
class CdsFeature:
    """A CDS feature checked once at the boundary: typed location and plain qualifiers."""

    contig_id: str
    location: Location
    qualifiers: dict[str, list[str]]
    nucleotides: str

    @classmethod
    def from_biopython(
        cls, feature: SeqFeature, location: Location, contig_id: str, record: SeqRecord
    ) -> "CdsFeature":
        """Copy what the parser needs out of a Biopython feature."""
        qualifiers = {key: [str(v) for v in values] for key, values in feature.qualifiers.items()}
        return cls(contig_id, location, qualifiers, str(feature.extract(record.seq)))

    def qualifier(self, key: str) -> str | None:
        """First value of a qualifier, or None."""
        values = self.qualifiers.get(key)
        return values[0] if values else None


@dataclass
class ParsedRecord:
    """Genes and skipped features from one GenBank record."""

    genes: list[Gene] = field(default_factory=list)
    skipped: list[SkippedFeature] = field(default_factory=list)
    trimmed: int = 0
    unparsed: int = 0


def read_genbank(path: str | Path) -> Genome:
    """Read every record of a GenBank file into a Genome, one contig per record."""
    contigs: dict[str, Contig] = {}
    skipped: list[SkippedFeature] = []
    for record in read_records(path):
        contig_id = record_id(record)
        if contig_id in contigs:
            raise ValueError(f"duplicate record ID {contig_id!r} in {path}")
        parsed = parse_record(record)
        contigs[contig_id] = record_contig(record, parsed.genes)
        skipped.extend(parsed.skipped)
    return Genome(contigs=contigs, source=str(path), skipped_features=tuple(skipped))


def read_reference_region(
    path: str | Path,
    *,
    contig: str | None = None,
    start: int | None = None,
    end: int | None = None,
) -> ReferenceRegion:
    """Read a reference region; ``start``/``end`` are 1-based inclusive (spec 2.1)."""
    record = _select_record(read_records(path), contig, path)
    sequence = record_sequence(record)
    region = _region(start, end, len(sequence))
    parsed = parse_record(record, region)
    return ReferenceRegion(
        contig_id=record_id(record),
        start=region[0],
        end=region[1],
        sequence=sequence[region[0] : region[1]],
        genes=tuple(parsed.genes),
        source_genome=str(path),
        skipped_features=tuple(parsed.skipped),
    )


def _select_record(records: list[SeqRecord], contig: str | None, path: str | Path) -> SeqRecord:
    if contig is None:
        if len(records) != 1:
            raise ValueError(f"{path} contains {len(records)} records; choose one with contig")
        return records[0]
    matches = [record for record in records if record_id(record) == contig]
    if not matches:
        raise ValueError(f"contig {contig!r} not found in {path}")
    if len(matches) > 1:
        raise ValueError(f"duplicate record ID {contig!r} in {path}")
    return matches[0]


def _region(start: int | None, end: int | None, record_length: int) -> Region:
    if start is None and end is None:
        return 0, record_length
    if start is None or end is None:
        raise ValueError("start and end must be given together")
    region = to_zero_based(start, end)
    if region[1] > record_length:
        raise ValueError(f"region end {end} lies beyond the record length {record_length}")
    return region


def read_records(path: str | Path) -> list[SeqRecord]:
    """Read all records of a GenBank file, closing the file before returning."""
    with open(path, encoding="utf-8") as handle:
        return list(SeqIO.parse(handle, "genbank"))


def record_id(record: SeqRecord) -> str:
    """The record's ID, which becomes the contig ID. Raises ValueError if missing."""
    if not record.id:
        raise ValueError("GenBank record without an ID")
    return str(record.id)


def record_sequence(record: SeqRecord) -> str:
    """The record's nucleotide sequence. Raises ValueError if it has none."""
    if record.seq is None:
        raise ValueError(f"GenBank record {record_id(record)!r} has no sequence")
    return str(record.seq)


def record_contig(record: SeqRecord, genes: Iterable[Gene]) -> Contig:
    """Build a Contig from a GenBank record and its parsed genes."""
    return Contig(
        id=record_id(record),
        sequence=record_sequence(record),
        genes=tuple(genes),
        description=record.description or None,
    )


def parse_record(record: SeqRecord, region: Region | None = None) -> ParsedRecord:
    """Parse the CDS features of a record, optionally only those inside ``region``.

    ``region`` is 0-based half-open. Issues one AnnotationWarning summarising
    anything skipped, untranslatable, trimmed or unparsed.
    """
    contig_id = record_id(record)
    parsed = ParsedRecord()
    for feature in record.features:
        if feature.type != "CDS":
            continue
        if feature.location is None:
            parsed.unparsed += 1
        elif _inside(feature.location, region):
            _add_cds(
                parsed, CdsFeature.from_biopython(feature, feature.location, contig_id, record)
            )
    _warn_summary(contig_id, parsed)
    return parsed


def _inside(location: Location, region: Region | None) -> bool:
    if region is None:
        return True
    return region[0] <= int(location.start) and int(location.end) <= region[1]


def _add_cds(parsed: ParsedRecord, cds: CdsFeature) -> None:
    reason = _skip_reason_before_translation(cds)
    translation = None if reason else _translation(cds)
    if isinstance(translation, InternalStop):
        reason = SkipReason.INTERNAL_STOP
    if reason:
        parsed.skipped.append(_skipped(cds, reason))
        return
    protein = translation if isinstance(translation, Protein) else None
    parsed.trimmed += bool(protein and protein.trimmed)
    parsed.genes.append(_gene(cds, protein))


def _skip_reason_before_translation(cds: CdsFeature) -> SkipReason | None:
    if isinstance(cds.location, CompoundLocation):
        return SkipReason.COMPOUND_LOCATION
    if "pseudo" in cds.qualifiers or "pseudogene" in cds.qualifiers:
        return SkipReason.PSEUDOGENE
    return None


def _translation(cds: CdsFeature) -> Translation:
    supplied = cds.qualifier("translation")
    if supplied is not None:
        checked = check_supplied_translation(supplied)
        if not isinstance(checked, InvalidSupplied):
            return checked
    table = _int_qualifier(cds, "transl_table", default=11)
    codon_start = _int_qualifier(cds, "codon_start", default=1)
    if table is None or codon_start is None:
        return Untranslatable("malformed /transl_table or /codon_start")
    return translate_cds(
        cds.nucleotides,
        table=table,
        codon_start=codon_start,
        complete_5prime=not _five_prime_partial(cds.location),
    )


def _gene(cds: CdsFeature, protein: Protein | None) -> Gene:
    return Gene(
        contig_id=cds.contig_id,
        start=int(cds.location.start),
        end=int(cds.location.end),
        strand=as_strand(cds.location.strand),
        nucleotide_sequence=cds.nucleotides,
        protein_sequence=protein.sequence if protein else None,
        gene_name=cds.qualifier("gene"),
        product=cds.qualifier("product"),
        locus_tag=cds.qualifier("locus_tag"),
        protein_id=cds.qualifier("protein_id"),
        partial=_partial(cds.location),
        qualifiers=cds.qualifiers,
    )


def _skipped(cds: CdsFeature, reason: SkipReason) -> SkippedFeature:
    return SkippedFeature(
        contig_id=cds.contig_id,
        start=int(cds.location.start),
        end=int(cds.location.end),
        locus_tag=cds.qualifier("locus_tag"),
        reason=reason,
    )


def _int_qualifier(cds: CdsFeature, key: str, *, default: int) -> int | None:
    """Integer qualifier value, ``default`` when absent, None when malformed."""
    value = cds.qualifier(key)
    if value is None:
        return default
    return int(value) if value.isdigit() else None


def _partial(location: Location) -> bool:
    return isinstance(location.start, BeforePosition) or isinstance(location.end, AfterPosition)


def _five_prime_partial(location: Location) -> bool:
    if location.strand == -1:
        return isinstance(location.end, AfterPosition)
    return isinstance(location.start, BeforePosition)


def _warn_summary(contig_id: str, parsed: ParsedRecord) -> None:
    parts = [f"{reason}: {count}" for reason, count in _skip_counts(parsed)]
    untranslatable = sum(gene.protein_sequence is None for gene in parsed.genes)
    if untranslatable:
        parts.append(f"{untranslatable} untranslatable")
    if parsed.trimmed:
        parts.append(f"{parsed.trimmed} trimmed to whole codons")
    if parsed.unparsed:
        parts.append(f"{parsed.unparsed} with an unparseable location")
    if parts:
        message = f"contig {contig_id}: " + "; ".join(parts)
        warnings.warn(message, AnnotationWarning, stacklevel=3)


def _skip_counts(parsed: ParsedRecord) -> list[tuple[str, int]]:
    return sorted(Counter(skipped.reason.value for skipped in parsed.skipped).items())
