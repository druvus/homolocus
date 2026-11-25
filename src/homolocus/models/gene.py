"""Protein-coding gene model (spec 4.1, 5.1)."""

from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Literal

from homolocus.coordinates import strand_char, to_one_based

Strand = Literal[1, -1]


def as_strand(value: object) -> Strand:
    """Parse an outside strand value; only 1 and -1 are accepted."""
    if value == 1:
        return 1
    if value == -1:
        return -1
    raise ValueError(f"invalid strand {value!r}; expected 1 or -1")


@dataclass(frozen=True)
class Gene:
    """A protein-coding CDS on a contig.

    Coordinates are 0-based half-open on the original contig. The nucleotide
    sequence is stored in coding orientation. Invalid genes cannot be built.
    """

    contig_id: str
    start: int
    end: int
    strand: Strand
    nucleotide_sequence: str
    protein_sequence: str | None

    gene_name: str | None = None
    product: str | None = None
    locus_tag: str | None = None
    protein_id: str | None = None
    partial: bool = False
    qualifiers: Mapping[str, list[str]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.contig_id:
            raise ValueError("gene contig ID must not be empty")
        to_one_based(self.start, self.end)
        strand_char(self.strand)
        self._check_sequences()

    def _check_sequences(self) -> None:
        if len(self.nucleotide_sequence) != self.length:
            raise ValueError(
                f"nucleotide sequence length {len(self.nucleotide_sequence)} "
                f"does not match gene length {self.length}"
            )
        if self.protein_sequence == "":
            raise ValueError("protein sequence must be None, not empty")

    @property
    def id(self) -> str:
        """Canonical gene ID ``{contig_id}|{start_1based}-{end_1based}|{+|-}``."""
        start_1, end_1 = to_one_based(self.start, self.end)
        return f"{self.contig_id}|{start_1}-{end_1}|{strand_char(self.strand)}"

    @property
    def length(self) -> int:
        """Length of the gene in base pairs."""
        return self.end - self.start


def genomic_order(gene: Gene) -> tuple[int, int, int]:
    """Sort key for genes on a contig: start, then end, then strand (spec 4)."""
    return gene.start, gene.end, gene.strand


def sorted_genes(
    genes: Iterable[Gene], *, contig_id: str, start: int, end: int, label: str
) -> tuple[Gene, ...]:
    """Check genes lie on one contig within ``[start, end)`` with unique IDs; sort them.

    ``label`` names the container in error messages.
    """
    checked = tuple(genes)
    for gene in checked:
        if gene.contig_id != contig_id:
            raise ValueError(f"gene {gene.id} belongs to contig {gene.contig_id}, not {contig_id}")
        if gene.start < start or gene.end > end:
            raise ValueError(f"gene {gene.id} extends past the bounds of {label}")
    duplicates = _duplicate_ids(checked)
    if duplicates:
        raise ValueError(f"{label} has duplicate gene IDs: {', '.join(duplicates)}")
    return tuple(sorted(checked, key=genomic_order))


def _duplicate_ids(genes: Iterable[Gene]) -> list[str]:
    """Gene IDs that occur more than once, sorted."""
    counts = Counter(gene.id for gene in genes)
    return sorted(gene_id for gene_id, count in counts.items() if count > 1)
