"""Candidate flags: contig edge, potential splits and fusions (spec 9, 10)."""

import pytest

from homolocus import Alignment, Contig, Gene, GeneMatch, RegionDetection
from homolocus.detection.candidate_generation import EventWindow
from homolocus.detection.flags import split_events, touches_contig_edge
from homolocus.models.cluster import SplitEvent

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


# Potential splits (spec 9.2, 9.4)

R1, R2 = "ref|1-300|+", "ref|401-700|+"


def target_gene(slot: int, strand: int = 1) -> Gene:
    start = 1000 + slot * 1000
    return Gene("t1", start, start + 300, strand, "A" * 300, "M" * 100)


T1, T2, T3, BETWEEN = target_gene(0), target_gene(1), target_gene(3), target_gene(2)


def piece(reference: str, target: Gene, *spans: tuple[int, int], length: int = 100) -> GeneMatch:
    hsps = tuple(
        Alignment(
            reference,
            target.id,
            80.0,
            end - start,
            start,
            end,
            0,
            end - start,
            length,
            100,
            1e-30,
            100.0,
        )
        for start, end in spans
    )
    return GeneMatch(reference, target.id, hsps)


def splits(*matches: GeneMatch, genes=(T1, T2, BETWEEN, T3), **settings) -> tuple[SplitEvent, ...]:
    contig = Contig("t1", "A" * LENGTH, genes)
    detection = RegionDetection(**settings)
    return split_events(contig, matches, detection)


def test_two_complementary_halves_are_a_split() -> None:
    (event,) = splits(piece(R1, T1, (0, 50)), piece(R1, T2, (50, 100)), genes=(T1, T2))
    assert event == SplitEvent(R1, (T1.id, T2.id), adjacent=True, same_strand=True)


def test_full_length_duplicate_copies_are_not_a_split() -> None:
    assert splits(piece(R1, T1, (0, 100)), piece(R1, T2, (0, 100))) == ()


def test_one_matching_gene_is_not_a_split() -> None:
    assert splits(piece(R1, T1, (0, 50))) == ()


@pytest.mark.parametrize(("start", "flagged"), [(80, True), (81, False)])
def test_unique_contribution_threshold_is_inclusive(start: int, flagged: bool) -> None:
    found = splits(piece(R1, T1, (0, start)), piece(R1, T2, (start - 5, 100)))
    assert bool(found) is flagged


@pytest.mark.parametrize(("overlap", "flagged"), [(6, True), (7, False)])
def test_fragment_overlap_limit_is_inclusive(overlap: int, flagged: bool) -> None:
    found = splits(piece(R1, T1, (0, 30)), piece(R1, T2, (30 - overlap, 100)))
    assert bool(found) is flagged


@pytest.mark.parametrize(("end", "flagged"), [(50, True), (49, False)])
def test_fragments_together_must_reach_the_represented_threshold(end: int, flagged: bool) -> None:
    found = splits(piece(R1, T1, (0, 25)), piece(R1, T2, (25, end)))
    assert bool(found) is flagged


def test_three_fragments_in_genomic_order() -> None:
    matches = (piece(R1, T3, (60, 100)), piece(R1, T1, (0, 30)), piece(R1, T2, (30, 60)))
    (event,) = splits(*matches)
    assert event.fragment_ids == (T1.id, T2.id, T3.id)


def test_hsps_of_one_fragment_are_merged() -> None:
    (event,) = splits(piece(R1, T1, (0, 20), (20, 50)), piece(R1, T2, (50, 100)))
    assert event.fragment_ids == (T1.id, T2.id)


def test_adjacency_and_strand_are_recorded() -> None:
    (apart,) = splits(piece(R1, T1, (0, 50)), piece(R1, T3, (50, 100)))
    assert not apart.adjacent
    reversed_t2 = target_gene(1, strand=-1)
    (mixed,) = splits(
        piece(R1, T1, (0, 50)), piece(R1, reversed_t2, (50, 100)), genes=(T1, reversed_t2)
    )
    assert mixed.adjacent
    assert not mixed.same_strand


def test_greedy_choice_keeps_the_best_non_overlapping_fragments() -> None:
    matches = (
        piece(R1, T1, (0, 40)),
        piece(R1, T2, (25, 70)),
        piece(R1, T3, (70, 100)),
    )
    (event,) = splits(*matches)
    assert event.fragment_ids == (T2.id, T3.id)


def test_each_reference_gene_is_judged_on_its_own() -> None:
    matches = (
        piece(R2, T1, (0, 50)),
        piece(R2, T2, (50, 100)),
        piece(R1, T3, (0, 100)),
        piece(R1, T1, (0, 100)),
    )
    (event,) = splits(*matches)
    assert event.reference_gene_id == R2


def test_events_follow_the_order_of_the_matches() -> None:
    matches = (
        piece(R2, T1, (0, 50)),
        piece(R2, T2, (50, 100)),
        piece(R1, T2, (0, 50)),
        piece(R1, T3, (50, 100)),
    )
    assert [e.reference_gene_id for e in splits(*matches)] == [R2, R1]


def test_custom_thresholds_are_used() -> None:
    matches = (piece(R1, T1, (0, 90)), piece(R1, T2, (85, 100)))
    assert splits(*matches) == ()
    assert splits(*matches, min_split_contribution=0.1, max_fragment_overlap=0.4) != ()


def test_adjacency_spans_all_fragments() -> None:
    matches = (piece(R1, T1, (0, 30)), piece(R1, T2, (30, 60)), piece(R1, T3, (60, 100)))
    (event,) = splits(*matches)
    assert not event.adjacent
    (inner,) = splits(piece(R1, T2, (0, 50)), piece(R1, BETWEEN, (50, 100)))
    assert inner.adjacent


def test_adjacency_is_false_when_the_first_fragments_are_set_apart() -> None:
    # Contig order: T1, T2, BETWEEN, T3. Fragments sit at positions 0, 2 and 3.
    matches = (piece(R1, T1, (0, 30)), piece(R1, BETWEEN, (30, 60)), piece(R1, T3, (60, 100)))
    (event,) = splits(*matches)
    assert event.fragment_ids == (T1.id, BETWEEN.id, T3.id)
    assert not event.adjacent


def test_small_separate_fragment_below_the_threshold_is_not_a_fragment() -> None:
    assert splits(piece(R1, T1, (0, 90)), piece(R1, T2, (90, 100))) == ()


# Potential fusions (spec 9.3, 9.4)

from homolocus import GeneFamily, ReferenceFamilies  # noqa: E402
from homolocus.detection.flags import fusion_events  # noqa: E402
from homolocus.models.cluster import FusionEvent  # noqa: E402

R3 = "ref|801-1100|+"
FAMILY_SET = ReferenceFamilies(
    (GeneFamily("F001", (R1,)), GeneFamily("F002", (R2,)), GeneFamily("F003", (R3,)))
)
PARALOGS = ReferenceFamilies((GeneFamily("F001", (R1, R2)), GeneFamily("F002", (R3,))))


def part(reference: str, target: Gene, *spans: tuple[int, int]) -> GeneMatch:
    hsps = tuple(
        Alignment(
            reference,
            target.id,
            80.0,
            end - start,
            0,
            end - start,
            start,
            end,
            100,
            100,
            1e-30,
            100.0,
        )
        for start, end in spans
    )
    return GeneMatch(reference, target.id, hsps)


def fusions(*matches: GeneMatch, families=FAMILY_SET, **settings) -> tuple[FusionEvent, ...]:
    contig = Contig("t1", "A" * LENGTH, (T1, T2, BETWEEN, T3))
    return fusion_events(contig, matches, families, RegionDetection(**settings))


def test_two_families_on_separate_parts_are_a_fusion() -> None:
    (event,) = fusions(part(R1, T1, (0, 50)), part(R2, T1, (50, 100)))
    assert event == FusionEvent(T1.id, ("F001", "F002"))


def test_families_are_listed_in_protein_order() -> None:
    (event,) = fusions(part(R1, T1, (50, 100)), part(R2, T1, (0, 50)))
    assert event.family_ids == ("F002", "F001")


def test_members_of_one_family_are_not_a_fusion() -> None:
    matches = (part(R1, T1, (0, 50)), part(R2, T1, (50, 100)))
    assert fusions(*matches, families=PARALOGS) == ()


def test_shared_domain_is_not_a_fusion() -> None:
    assert fusions(part(R1, T1, (0, 60)), part(R2, T1, (10, 70))) == ()


@pytest.mark.parametrize(("start", "flagged"), [(40, True), (39, False)])
def test_part_overlap_limit_is_inclusive(start: int, flagged: bool) -> None:
    found = fusions(part(R1, T1, (0, 50)), part(R2, T1, (start, 100)))
    assert bool(found) is flagged


@pytest.mark.parametrize(("end", "flagged"), [(100, True), (99, False)])
def test_part_coverage_threshold_is_inclusive(end: int, flagged: bool) -> None:
    found = fusions(part(R1, T1, (0, 60)), part(R2, T1, (80, end)))
    assert bool(found) is flagged


def test_a_third_overlapping_family_does_not_hide_a_fusion() -> None:
    matches = (part(R3, T1, (0, 45)), part(R1, T1, (0, 50)), part(R2, T1, (50, 100)))
    (event,) = fusions(*matches)
    assert event.family_ids == ("F001", "F002")


def test_hsps_of_one_family_are_merged_into_one_part() -> None:
    matches = (part(R1, T1, (0, 30)), part(R2, T1, (30, 60)), part(R3, T1, (60, 100)))
    (event,) = fusions(*matches, families=PARALOGS)
    assert event.family_ids == ("F001", "F002")


def test_one_family_alone_is_not_a_fusion() -> None:
    assert fusions(part(R1, T1, (0, 50), (50, 100))) == ()


def test_fusions_are_reported_in_genomic_order() -> None:
    matches = (
        part(R1, T3, (0, 50)),
        part(R2, T3, (50, 100)),
        part(R1, T1, (0, 50)),
        part(R3, T1, (50, 100)),
    )
    assert [e.target_gene_id for e in fusions(*matches)] == [T1.id, T3.id]


def test_custom_part_coverage_is_used() -> None:
    matches = (part(R1, T1, (0, 70)), part(R2, T1, (85, 100)))
    assert fusions(*matches) == ()
    assert fusions(*matches, min_fusion_part_coverage=0.15) != ()


def test_equal_parts_in_conflict_resolve_by_protein_position() -> None:
    matches = (part(R1, T1, (30, 80)), part(R2, T1, (0, 50)), part(R3, T1, (80, 100)))
    (event,) = fusions(*matches)
    assert event.family_ids == ("F002", "F003")


def test_order_follows_the_protein_not_the_part_size() -> None:
    (event,) = fusions(part(R1, T1, (0, 30)), part(R2, T1, (30, 100)))
    assert event.family_ids == ("F001", "F002")
