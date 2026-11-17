"""GeneMatch: aggregated evidence for one reference-target protein pair (spec 5.7)."""

from collections.abc import Iterable
from dataclasses import dataclass

from homolocus.intervals import Interval, union_length
from homolocus.models.alignment import Alignment


@dataclass(frozen=True)
class GeneMatch:
    """All HSPs between one reference protein and one target protein.

    Every HSP is kept. Coverage, identity, bit score and e-value are derived
    from the HSPs, so they cannot disagree with them.
    """

    reference_gene_id: str
    target_gene_id: str
    alignments: tuple[Alignment, ...]

    def __post_init__(self) -> None:
        if not self.alignments:
            raise ValueError("a gene match needs at least one alignment")
        for hsp in self.alignments:
            if hsp.query_id != self.reference_gene_id:
                raise ValueError(f"alignment query {hsp.query_id} is not {self.reference_gene_id}")
            if hsp.subject_id != self.target_gene_id:
                raise ValueError(f"alignment subject {hsp.subject_id} is not {self.target_gene_id}")
        if len({hsp.query_length for hsp in self.alignments}) > 1:
            raise ValueError("alignments disagree on the query length")
        if len({hsp.subject_length for hsp in self.alignments}) > 1:
            raise ValueError("alignments disagree on the subject length")

    @property
    def query_length(self) -> int:
        """Length of the reference protein."""
        return self.alignments[0].query_length

    @property
    def subject_length(self) -> int:
        """Length of the target protein."""
        return self.alignments[0].subject_length

    @property
    def query_coverage(self) -> float:
        """Fraction of the reference protein covered by the union of HSPs."""
        return _coverage((hsp.query_interval for hsp in self.alignments), self.query_length)

    @property
    def subject_coverage(self) -> float:
        """Fraction of the target protein covered by the union of HSPs."""
        return _coverage((hsp.subject_interval for hsp in self.alignments), self.subject_length)

    @property
    def identity(self) -> float:
        """Alignment-length-weighted identity, in percent."""
        weighted = sum(hsp.identity * hsp.alignment_length for hsp in self.alignments)
        return weighted / sum(hsp.alignment_length for hsp in self.alignments)

    @property
    def bit_score(self) -> float:
        """Sum of the HSP bit scores."""
        return sum(hsp.bit_score for hsp in self.alignments)

    @property
    def evalue(self) -> float:
        """Best (smallest) HSP e-value."""
        return min(hsp.evalue for hsp in self.alignments)

    def passes(self, *, max_evalue: float, min_identity: float) -> bool:
        """Pair-level thresholds on the aggregated evidence; coverage is not used."""
        return self.evalue <= max_evalue and self.identity >= min_identity


def _coverage(intervals: Iterable[Interval], length: int) -> float:
    return union_length(intervals) / length


def aggregate(alignments: Iterable[Alignment]) -> tuple[GeneMatch, ...]:
    """Group HSPs into one GeneMatch per protein pair.

    Pairs appear in the order they are first seen; HSPs keep input order.
    """
    groups: dict[tuple[str, str], list[Alignment]] = {}
    for hsp in alignments:
        groups.setdefault((hsp.query_id, hsp.subject_id), []).append(hsp)
    return tuple(
        GeneMatch(query, subject, tuple(hsps)) for (query, subject), hsps in groups.items()
    )
