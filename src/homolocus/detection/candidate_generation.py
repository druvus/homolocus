"""Candidate generation on one contig: chains, windows and padding (spec 7.3, 7.4, 10).

Events are matching-event genes. Gaps and spans use the furthest end reached
so far, so nested or overlapping genes never shorten a candidate.
"""

from collections.abc import Iterable
from dataclasses import dataclass

from homolocus.config import RegionDetection
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


def windows(
    contig: Contig, events: Iterable[Gene], config: RegionDetection
) -> tuple[EventWindow, ...]:
    """Candidate windows on one contig, ordered by start then end (spec 7.3).

    A chain whose span is at most ``max_region_bp`` is one window. A longer
    chain is split: each event anchors a window that extends right while the
    span stays within the limit, and windows contained in another are dropped.
    """
    found: list[tuple[Gene, ...]] = []
    for chain in chains(
        contig,
        events,
        max_intervening_genes=config.max_intervening_genes,
        max_inter_match_bp_gap=config.max_inter_match_bp_gap,
    ):
        found.extend(_split(chain, config.max_region_bp))
    kept = [EventWindow(contig.id, group) for group in _maximal(found)]
    return tuple(sorted(kept, key=lambda window: (window.start, window.end, window.event_ids)))


def _split(chain: tuple[Gene, ...], max_region_bp: int) -> list[tuple[Gene, ...]]:
    groups: list[tuple[Gene, ...]] = []
    for anchor, first in enumerate(chain):
        group: tuple[Gene, ...] = (first,)
        for event in chain[anchor + 1 :]:
            if _span((*group, event)) > max_region_bp:
                break
            group = (*group, event)
        groups.append(group)
    return groups


def _span(events: tuple[Gene, ...]) -> int:
    return max(event.end for event in events) - events[0].start


def _maximal(groups: list[tuple[Gene, ...]]) -> list[tuple[Gene, ...]]:
    unique = list({frozenset(g.id for g in group): group for group in groups}.items())
    return [group for ids, group in unique if not any(ids < other for other, _ in unique)]


def padded_interval(
    window: EventWindow, contig: Contig, *, padding_genes: int, padding_bp: int
) -> tuple[int, int]:
    """Displayed interval of a window: the event span widened by padding (spec 10).

    ``padding_genes`` adds whole genes from the contig's sorted gene list on
    each side of the outermost events; ``padding_bp`` adds base pairs. When
    both are given, the wider extent wins on each side. The result is clipped
    to the contig. Padding never changes which events form the window.
    """
    if window.contig_id != contig.id:
        raise ValueError(f"window is on contig {window.contig_id}, not {contig.id}")
    position = {gene.id: index for index, gene in enumerate(contig.genes)}
    for event in window.events:
        if event.id not in position:
            raise ValueError(f"event {event.id} is not a gene of contig {contig.id}")
    first, last = position[window.events[0].id], position[window.events[-1].id]
    around = contig.genes[max(first - padding_genes, 0) : last + 1 + padding_genes]
    start = min(window.start - padding_bp, *(gene.start for gene in around))
    end = max(window.end + padding_bp, *(gene.end for gene in around))
    return max(start, 0), min(end, contig.length)
