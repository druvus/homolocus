"""Candidate generation on one contig: chains, windows and padding (spec 7.3, 7.4, 10).

Events are matching-event genes. Gaps and spans use the furthest end reached
so far, so nested or overlapping genes never shorten a candidate.
"""

from collections.abc import Iterable
from dataclasses import dataclass

from homolocus.models.contig import Contig
from homolocus.models.gene import Gene, genomic_order


@dataclass(frozen=True)
class EventWindow:
    """Matching events on one contig evaluated together as one candidate."""

    contig_id: str
    events: tuple[Gene, ...]

    def __post_init__(self) -> None:
        if not self.events:
            raise ValueError("a window needs at least one event")
        if any(event.contig_id != self.contig_id for event in self.events):
            raise ValueError(f"window events must lie on contig {self.contig_id}")
        keys = [genomic_order(event) for event in self.events]
        if any(left >= right for left, right in zip(keys, keys[1:], strict=False)):
            raise ValueError("window events must be distinct and in genomic order")

    @property
    def start(self) -> int:
        """Start of the first event (0-based)."""
        return self.events[0].start

    @property
    def end(self) -> int:
        """Furthest end of any event (exclusive)."""
        return max(event.end for event in self.events)

    @property
    def span(self) -> int:
        """Unpadded length in base pairs, the value checked against max_region_bp."""
        return self.end - self.start

    @property
    def event_ids(self) -> tuple[str, ...]:
        """Canonical IDs of the events."""
        return tuple(event.id for event in self.events)


@dataclass
class _Chain:
    events: list[Gene]
    reach: int


def chains(
    contig: Contig,
    events: Iterable[Gene],
    *,
    max_intervening_genes: int,
    max_inter_match_bp_gap: int,
) -> tuple[tuple[Gene, ...], ...]:
    """Single-linkage chains of events in genomic order (spec 7.3).

    Consecutive events are linked when at most ``max_intervening_genes``
    contig genes lie between them and the next start is at most
    ``max_inter_match_bp_gap`` past the furthest end so far (minimum 0).
    """
    position = {gene.id: index for index, gene in enumerate(contig.genes)}
    groups: list[_Chain] = []
    for event in sorted(events, key=genomic_order):
        if event.id not in position:
            raise ValueError(f"event {event.id} is not a gene of contig {contig.id}")
        last = groups[-1] if groups else None
        if last and _joins(last, event, position, max_intervening_genes, max_inter_match_bp_gap):
            last.events.append(event)
            last.reach = max(last.reach, event.end)
        else:
            groups.append(_Chain([event], event.end))
    return tuple(tuple(chain.events) for chain in groups)


def _joins(
    chain: _Chain, event: Gene, position: dict[str, int], max_genes: int, max_gap: int
) -> bool:
    intervening = position[event.id] - position[chain.events[-1].id] - 1
    return intervening <= max_genes and max(event.start - chain.reach, 0) <= max_gap
