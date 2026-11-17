"""CDS features excluded during annotation parsing (spec 3.2, 5.10)."""

from dataclasses import dataclass
from enum import Enum

from homolocus.coordinates import to_one_based


class SkipReason(str, Enum):
    """Why a CDS feature was excluded."""

    PSEUDOGENE = "pseudogene"
    INTERNAL_STOP = "internal_stop"
    COMPOUND_LOCATION = "compound_location"


@dataclass(frozen=True)
class SkippedFeature:
    """A CDS excluded during parsing, kept for the run metadata.

    ``start`` and ``end`` are 0-based half-open and span the whole location.
    """

    contig_id: str
    start: int
    end: int
    locus_tag: str | None
    reason: SkipReason

    def __post_init__(self) -> None:
        if not self.contig_id:
            raise ValueError("skipped feature contig ID must not be empty")
        to_one_based(self.start, self.end)
