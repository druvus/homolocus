"""Candidate flags: contig edge, potential splits and fusions (spec 9, 10). Pure functions."""

from collections.abc import Iterable, Mapping, Sequence

from homolocus.config import RegionDetection
from homolocus.detection.candidate_generation import EventWindow
from homolocus.intervals import Interval, merge_intervals, union_length
from homolocus.models.cluster import SplitEvent
from homolocus.models.contig import Contig
from homolocus.models.gene import Gene
from homolocus.models.gene_match import GeneMatch

Part = tuple[str, tuple[Interval, ...]]
"""A labelled set of merged intervals: a split fragment or a fused part."""


def touches_contig_edge(window: EventWindow, contig: Contig, *, edge_margin_bp: int) -> bool:
    """Whether a candidate may be truncated by the end of its contig (spec 10).

    True when a boundary event (the first, or the one reaching furthest right)
    is partial, or lies within ``edge_margin_bp`` (inclusive) of a contig end.
    """
    first = window.events[0]
    last = max(window.events, key=lambda event: event.end)
    near_start = first.start <= edge_margin_bp
    near_end = contig.length - last.end <= edge_margin_bp
    return first.partial or last.partial or near_start or near_end


def split_events(
    contig: Contig, matches: Iterable[GeneMatch], detection: RegionDetection
) -> tuple[SplitEvent, ...]:
    """Reference proteins covered by complementary target fragments (spec 9.2, 9.4).

    ``matches`` are a candidate's passing evidence matches. For each reference
    gene, in order of first appearance, fragments whose unique contribution is
    at least ``min_split_contribution`` are chosen greedily (largest first,
    overlap at most ``max_fragment_overlap`` of the shorter). Two or more chosen
    fragments whose union reaches ``min_reference_gene_coverage`` form a split.
    """
    pieces: dict[str, dict[str, list[Interval]]] = {}
    lengths: dict[str, int] = {}
    for match in matches:
        target_pieces = pieces.setdefault(match.reference_gene_id, {})
        spans = target_pieces.setdefault(match.target_gene_id, [])
        spans.extend(hsp.query_interval for hsp in match.alignments)
        lengths[match.reference_gene_id] = match.query_length
    genes = {gene.id: gene for gene in contig.genes}
    position = {gene.id: index for index, gene in enumerate(contig.genes)}
    found = []
    for reference_id, by_target in pieces.items():
        merged = {target: merge_intervals(spans) for target, spans in by_target.items()}
        fragments = _split_fragments(merged, lengths[reference_id], position, detection)
        if fragments:
            found.append(_split_event(reference_id, fragments, genes, position))
    return tuple(found)


def _split_fragments(
    merged: Mapping[str, tuple[Interval, ...]],
    length: int,
    position: Mapping[str, int],
    detection: RegionDetection,
) -> list[str]:
    unique = {target: _unique_length(target, merged) for target in merged}
    candidates = sorted(
        (
            part
            for part in merged.items()
            if unique[part[0]] / length >= detection.min_split_contribution
        ),
        key=lambda part: (-unique[part[0]], position[part[0]]),
    )
    chosen = _choose(candidates, detection.max_fragment_overlap)
    covered = union_length(interval for _, spans in chosen for interval in spans)
    if len(chosen) < 2 or covered / length < detection.min_reference_gene_coverage:
        return []
    return sorted((target for target, _ in chosen), key=position.__getitem__)


def _split_event(
    reference_id: str, fragments: list[str], genes: Mapping[str, Gene], position: Mapping[str, int]
) -> SplitEvent:
    first, last = position[fragments[0]], position[fragments[-1]]
    return SplitEvent(
        reference_gene_id=reference_id,
        fragment_ids=tuple(fragments),
        adjacent=last - first == len(fragments) - 1,
        same_strand=len({genes[fragment].strand for fragment in fragments}) == 1,
    )


def _unique_length(target: str, merged: Mapping[str, tuple[Interval, ...]]) -> int:
    everything = [interval for spans in merged.values() for interval in spans]
    others = [interval for key, spans in merged.items() if key != target for interval in spans]
    return union_length(everything) - union_length(others)


def _choose(parts: Sequence[Part], max_overlap: float) -> list[Part]:
    """Greedy choice in the given order, skipping a part that overlaps a chosen one too much."""
    chosen: list[Part] = []
    for part in parts:
        if all(_overlap_fraction(part[1], other[1]) <= max_overlap for other in chosen):
            chosen.append(part)
    return chosen


def _overlap_fraction(first: tuple[Interval, ...], second: tuple[Interval, ...]) -> float:
    """Shared length of two interval sets as a fraction of the shorter set."""
    size_first, size_second = union_length(first), union_length(second)
    shared = size_first + size_second - union_length((*first, *second))
    return shared / min(size_first, size_second)
