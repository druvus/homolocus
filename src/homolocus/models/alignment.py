"""Alignment model: one DIAMOND HSP (spec 5.6)."""

from dataclasses import dataclass

from homolocus.coordinates import to_one_based
from homolocus.intervals import Interval


@dataclass(frozen=True)
class Alignment:
    """One high-scoring segment pair between a reference and a target protein.

    IDs are canonical gene IDs. Intervals are 0-based half-open protein
    positions. Identity is a percentage. Invalid HSPs cannot be built.
    """

    query_id: str
    subject_id: str

    identity: float
    alignment_length: int

    query_start: int
    query_end: int
    subject_start: int
    subject_end: int

    query_length: int
    subject_length: int

    evalue: float
    bit_score: float

    def __post_init__(self) -> None:
        if not self.query_id:
            raise ValueError("alignment query ID must not be empty")
        if not self.subject_id:
            raise ValueError("alignment subject ID must not be empty")
        _check_span("query", self.query_interval, self.query_length)
        _check_span("subject", self.subject_interval, self.subject_length)
        self._check_scores()

    def _check_scores(self) -> None:
        if not 0 <= self.identity <= 100:
            raise ValueError(f"identity {self.identity} is not a percentage")
        longest_span = max(self.query_end - self.query_start, self.subject_end - self.subject_start)
        if self.alignment_length < longest_span:
            raise ValueError(
                f"alignment length {self.alignment_length} is shorter than an aligned span"
            )
        if self.evalue < 0:
            raise ValueError(f"e-value {self.evalue} is negative")
        if self.bit_score < 0:
            raise ValueError(f"bit score {self.bit_score} is negative")

    @property
    def query_interval(self) -> Interval:
        """Aligned reference-protein positions, 0-based half-open."""
        return self.query_start, self.query_end

    @property
    def subject_interval(self) -> Interval:
        """Aligned target-protein positions, 0-based half-open."""
        return self.subject_start, self.subject_end


def _check_span(side: str, interval: Interval, length: int) -> None:
    if length <= 0:
        raise ValueError(f"{side} length {length} must be positive")
    to_one_based(*interval)
    if interval[1] > length:
        raise ValueError(f"{side} interval {interval} extends past the {side} length {length}")
