"""Candidate flags (spec 9, 10). Split and fusion flags follow in milestone 9."""

from homolocus.detection.candidate_generation import EventWindow
from homolocus.models.contig import Contig


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
