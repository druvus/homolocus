"""Chains, windows and padding on one contig (spec 7.3, 7.4, 10)."""

import pytest

from homolocus import Contig, Gene
from homolocus.detection.candidate_generation import EventWindow, chains


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
