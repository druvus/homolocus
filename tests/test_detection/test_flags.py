"""Candidate flags: contig edge (spec 10)."""

import pytest

from homolocus import Contig, Gene
from homolocus.detection.candidate_generation import EventWindow
from homolocus.detection.flags import touches_contig_edge

LENGTH = 10_000


def gene(start: int, end: int, partial: bool = False) -> Gene:
    return Gene("t1", start, end, 1, "A" * (end - start), "M", partial=partial)


def flag(*events: Gene, margin: int = 500) -> bool:
    contig = Contig("t1", "A" * LENGTH, events)
    return touches_contig_edge(EventWindow("t1", events), contig, edge_margin_bp=margin)


def test_candidate_in_the_middle_does_not_touch_an_edge() -> None:
    assert not flag(gene(2000, 3000), gene(4000, 5000))


@pytest.mark.parametrize(("start", "touches"), [(500, True), (501, False)])
def test_margin_at_contig_start_is_inclusive(start: int, touches: bool) -> None:
    assert flag(gene(start, start + 1000), gene(4000, 5000)) is touches


@pytest.mark.parametrize(("end", "touches"), [(9500, True), (9499, False)])
def test_margin_at_contig_end_is_inclusive(end: int, touches: bool) -> None:
    assert flag(gene(2000, 3000), gene(end - 1000, end)) is touches


def test_right_edge_uses_the_furthest_reaching_event() -> None:
    long_event, nested = gene(2000, 9600), gene(3000, 3500)
    assert flag(long_event, nested)


@pytest.mark.parametrize("position", [0, 1])
def test_partial_boundary_event_touches_an_edge(position: int) -> None:
    events = [gene(2000, 3000), gene(4000, 5000)]
    events[position] = gene(events[position].start, events[position].end, partial=True)
    assert flag(*events)


def test_partial_inner_event_does_not_count() -> None:
    assert not flag(gene(2000, 3000), gene(3500, 3800, partial=True), gene(4000, 5000))


def test_zero_margin_flags_only_genes_at_the_ends() -> None:
    assert flag(gene(0, 1000), gene(4000, 5000), margin=0)
    assert not flag(gene(1, 1000), gene(4000, 5000), margin=0)
