"""Chains, windows and padding on one contig (spec 7.3, 7.4, 10)."""

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from homolocus import Contig, Gene, RegionDetection
from homolocus.detection.candidate_generation import EventWindow, chains, padded_interval, windows


def gene(start: int, end: int, contig_id: str = "t1", strand: int = 1) -> Gene:
    return Gene(contig_id, start, end, strand, "A" * (end - start), "M")


def contig_of(*genes: Gene, length: int = 200_000) -> Contig:
    return Contig("t1", "A" * length, genes)


def chain_ids(contig: Contig, events, *, genes: int = 10, gap: int = 20_000) -> list[list[str]]:
    found = chains(contig, events, max_intervening_genes=genes, max_inter_match_bp_gap=gap)
    return [[g.id for g in chain] for chain in found]


# EventWindow


def test_window_span_reaches_the_furthest_end() -> None:
    long_gene, nested = gene(0, 10_000), gene(100, 500)
    window = EventWindow("t1", (long_gene, nested))
    assert (window.start, window.end, window.span) == (0, 10_000, 10_000)
    assert window.event_ids == (long_gene.id, nested.id)


def test_window_span_is_measured_from_its_start() -> None:
    window = EventWindow("t1", (gene(3000, 4000), gene(6000, 7500)))
    assert (window.start, window.end, window.span) == (3000, 7500, 4500)


@pytest.mark.parametrize(
    ("contig_id", "events", "message"),
    [
        ("t1", (), "at least one event"),
        ("t1", (gene(500, 600), gene(0, 100)), "genomic order"),
        ("t1", (gene(0, 100), gene(0, 100)), "genomic order"),
        ("t2", (gene(0, 100),), "contig"),
    ],
)
def test_illegal_windows_are_rejected(contig_id, events, message) -> None:
    with pytest.raises(ValueError, match=message):
        EventWindow(contig_id, events)


# Chains


def test_events_within_both_limits_form_one_chain() -> None:
    a, b, c = gene(0, 1000), gene(5000, 6000), gene(9000, 10_000)
    assert chain_ids(contig_of(a, b, c), (a, b, c)) == [[a.id, b.id, c.id]]


def test_bp_gap_limit_is_inclusive() -> None:
    a, b = gene(0, 1000), gene(21_000, 22_000)
    assert len(chain_ids(contig_of(a, b), (a, b), gap=20_000)) == 1
    assert len(chain_ids(contig_of(a, b), (a, b), gap=19_999)) == 2


def test_intervening_gene_limit_is_inclusive() -> None:
    a, x1, x2, b = gene(0, 100), gene(200, 300), gene(400, 500), gene(600, 700)
    contig = contig_of(a, x1, x2, b)
    assert len(chain_ids(contig, (a, b), genes=2)) == 1
    assert len(chain_ids(contig, (a, b), genes=1)) == 2


def test_overlapping_events_have_zero_gap() -> None:
    a, b = gene(0, 1000), gene(900, 2000)
    assert len(chain_ids(contig_of(a, b), (a, b), gap=0)) == 1


def test_gap_is_measured_from_the_furthest_end_so_far() -> None:
    long_gene, nested, nxt = gene(0, 10_000), gene(100, 500), gene(12_000, 13_000)
    contig = contig_of(long_gene, nested, nxt)
    assert len(chain_ids(contig, (long_gene, nested, nxt), gap=2_000)) == 1
    assert len(chain_ids(contig, (long_gene, nested, nxt), gap=1_999)) == 2


def test_insertions_count_as_intervening_genes_and_copies_are_events() -> None:
    a, b1, x, b2, c, d = (gene(n * 1000, n * 1000 + 800) for n in range(6))
    contig = contig_of(a, b1, x, b2, c, d)
    assert chain_ids(contig, (a, b1, b2, c, d), genes=1) == [[a.id, b1.id, b2.id, c.id, d.id]]
    assert chain_ids(contig, (a, b1, b2, c, d), genes=0) == [[a.id, b1.id], [b2.id, c.id, d.id]]


def test_events_are_taken_in_genomic_order() -> None:
    a, b = gene(0, 1000), gene(2000, 3000)
    assert chain_ids(contig_of(a, b), (b, a)) == [[a.id, b.id]]


def test_no_events_give_no_chains() -> None:
    assert chain_ids(contig_of(gene(0, 100)), ()) == []


def test_event_must_be_a_gene_of_the_contig() -> None:
    with pytest.raises(ValueError, match="not a gene of contig"):
        chains(
            contig_of(gene(0, 100)),
            (gene(500, 600),),
            max_intervening_genes=1,
            max_inter_match_bp_gap=1,
        )


# Windows


def window_ids(contig: Contig, events, **config) -> list[tuple[str, ...]]:
    return [w.event_ids for w in windows(contig, events, RegionDetection(**config))]


def spaced(count: int, step: int = 3000, length: int = 1000) -> list[Gene]:
    return [gene(n * step, n * step + length) for n in range(count)]


def test_chain_within_max_region_is_one_window() -> None:
    events = spaced(4)
    assert window_ids(contig_of(*events), events, max_region_bp=10_000) == [
        tuple(g.id for g in events)
    ]


def test_max_region_is_inclusive() -> None:
    events = spaced(4)  # span 0..10_000
    assert len(window_ids(contig_of(*events), events, max_region_bp=10_000)) == 1
    assert len(window_ids(contig_of(*events), events, max_region_bp=9_999)) > 1


def test_long_chain_is_split_into_maximal_windows() -> None:
    a, b, c, d, e = spaced(5)  # events every 3 kb, 1 kb long
    found = window_ids(contig_of(a, b, c, d, e), (a, b, c, d, e), max_region_bp=7_000)
    assert found == [(a.id, b.id, c.id), (b.id, c.id, d.id), (c.id, d.id, e.id)]


def test_windows_contained_in_another_window_are_dropped() -> None:
    a, b, c = gene(0, 1000), gene(2000, 9000), gene(9500, 10_000)
    found = window_ids(contig_of(a, b, c), (a, b, c), max_region_bp=9_000)
    assert found == [(a.id, b.id), (b.id, c.id)]


def test_single_event_longer_than_max_region_is_its_own_window() -> None:
    huge, small = gene(0, 60_000), gene(61_000, 62_000)
    found = window_ids(contig_of(huge, small), (huge, small), max_region_bp=50_000)
    assert found == [(huge.id,), (small.id,)]


def test_separate_chains_give_separate_windows_in_order() -> None:
    a, b, c, d = gene(0, 1000), gene(2000, 3000), gene(60_000, 61_000), gene(62_000, 63_000)
    found = window_ids(contig_of(a, b, c, d), (d, c, b, a), max_inter_match_bp_gap=20_000)
    assert found == [(a.id, b.id), (c.id, d.id)]


def test_chain_limits_come_from_the_configuration() -> None:
    a, x, b = gene(0, 1000), gene(1500, 1800), gene(2000, 3000)
    contig = contig_of(a, x, b)
    assert len(window_ids(contig, (a, b), max_intervening_genes=0)) == 2
    assert len(window_ids(contig, (a, b), max_intervening_genes=1)) == 1


def test_no_events_give_no_windows() -> None:
    assert window_ids(contig_of(gene(0, 100)), ()) == []


@st.composite
def contig_with_events(draw):
    lengths = draw(st.lists(st.integers(100, 3000), min_size=1, max_size=25))
    gaps = draw(st.lists(st.integers(-50, 8000), min_size=len(lengths), max_size=len(lengths)))
    genes, position = [], 0
    for length, gap in zip(lengths, gaps, strict=True):
        start = max(position + gap, 0)
        genes.append(gene(start, start + length))
        position = start + length
    flags = draw(st.lists(st.booleans(), min_size=len(genes), max_size=len(genes)))
    events = [g for g, chosen in zip(genes, flags, strict=True) if chosen]
    unique = list({g.id: g for g in genes}.values())
    return contig_of(*unique, length=position + 1), [g for g in events if g in unique]


@settings(max_examples=150, deadline=None)
@given(
    contig_with_events(),
    st.integers(500, 20_000),
    st.integers(0, 5),
    st.integers(0, 10_000),
)
def test_window_properties(case, max_region, max_genes, max_gap) -> None:
    contig, events = case
    config = RegionDetection(
        max_region_bp=max_region,
        max_intervening_genes=max_genes,
        max_inter_match_bp_gap=max_gap,
    )
    found = windows(contig, events, config)
    event_sets = [set(w.event_ids) for w in found]
    assert {g.id for g in events} == set().union(*event_sets)
    assert all(w.span <= max_region or len(w.events) == 1 for w in found)
    assert not any(
        i != j and left <= right
        for i, left in enumerate(event_sets)
        for j, right in enumerate(event_sets)
    )
    assert [(w.start, w.end) for w in found] == sorted((w.start, w.end) for w in found)


# Padding

# Genes at 1-2 kb, 4-5 kb, [events 10-11 kb and 13-14 kb], 20-21 kb, 30-31 kb on a 40 kb contig.
LEFT_2, LEFT_1, EVENT_1, EVENT_2, RIGHT_1, RIGHT_2 = (
    gene(1000, 2000),
    gene(4000, 5000),
    gene(10_000, 11_000),
    gene(13_000, 14_000),
    gene(20_000, 21_000),
    gene(30_000, 31_000),
)
PADDED_CONTIG = contig_of(LEFT_2, LEFT_1, EVENT_1, EVENT_2, RIGHT_1, RIGHT_2, length=40_000)
WINDOW = EventWindow("t1", (EVENT_1, EVENT_2))


def pad(genes: int = 0, bp: int = 0, window: EventWindow = WINDOW) -> tuple[int, int]:
    return padded_interval(window, PADDED_CONTIG, padding_genes=genes, padding_bp=bp)


def test_no_padding_is_the_event_span() -> None:
    assert pad() == (10_000, 14_000)


@pytest.mark.parametrize(
    ("genes", "interval"),
    [(1, (4000, 21_000)), (2, (1000, 31_000)), (5, (1000, 31_000))],
)
def test_gene_padding_adds_whole_neighbouring_genes(genes: int, interval) -> None:
    assert pad(genes=genes) == interval


@pytest.mark.parametrize(
    ("bp", "interval"),
    [(500, (9500, 14_500)), (30_000, (0, 40_000))],
)
def test_bp_padding_is_clipped_to_the_contig(bp: int, interval) -> None:
    assert pad(bp=bp) == interval


def test_larger_extent_wins_on_each_side_when_both_are_given() -> None:
    assert pad(genes=1, bp=8000) == (2000, 22_000)


def test_gene_padding_never_shrinks_a_window_with_a_long_nested_event() -> None:
    long_event, nested = gene(0, 9000), gene(500, 800)
    contig = contig_of(long_event, nested, gene(9500, 9800), length=10_000)
    window = EventWindow("t1", (long_event, nested))
    assert padded_interval(window, contig, padding_genes=1, padding_bp=0) == (0, 9800)
    assert padded_interval(window, contig, padding_genes=0, padding_bp=0) == (0, 9000)


def test_gene_padding_counts_from_the_last_event() -> None:
    window = EventWindow("t1", (LEFT_1, EVENT_1, EVENT_2))
    assert pad(genes=1, window=window) == (1000, 21_000)


def test_window_events_must_be_genes_of_the_contig() -> None:
    other = EventWindow("t1", (gene(50_000, 51_000),))
    with pytest.raises(ValueError, match=r"event t1\|50001-51000\|\+ is not a gene of contig t1"):
        pad(window=other)


def test_window_must_be_on_the_same_contig() -> None:
    other = EventWindow("t2", (gene(10_000, 11_000, contig_id="t2"),))
    with pytest.raises(ValueError, match="window is on contig t2, not t1"):
        pad(window=other)
