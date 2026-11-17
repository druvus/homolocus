"""Score components (spec 8)."""

import pytest
from hypothesis import given
from hypothesis import strategies as st

from homolocus import Gene, ScoreComponents
from homolocus.detection.scoring import FamilyPairing, components


def gene(contig: str, start: int, strand: int = 1, length: int = 900) -> Gene:
    return Gene(contig, start, start + length, strand, "A" * length, "M")


def pairing(
    ref_start: int,
    target_start: int,
    *,
    ref_strand: int = 1,
    target_strand: int = 1,
    identity: float = 80.0,
    coverage: float = 1.0,
) -> FamilyPairing:
    return FamilyPairing(
        member=gene("ref", ref_start, ref_strand),
        copy=gene("tgt", target_start, target_strand),
        identity=identity,
        coverage=coverage,
    )


def score(pairings, *, total: int | None = None, span: int | None = None) -> ScoreComponents:
    total = len(pairings) if total is None else total
    if span is None:
        span = (
            max(p.copy.end for p in pairings) - min(p.copy.start for p in pairings)
            if pairings
            else 1
        )
    return components(pairings, total_families=total, candidate_span=span)


COLINEAR = [pairing(n * 1000, 50_000 + n * 1000) for n in range(4)]


def test_colinear_region_scores_one_on_structure() -> None:
    result = score(COLINEAR)
    assert (result.content, result.order, result.orientation, result.spacing) == (
        1.0,
        1.0,
        1.0,
        1.0,
    )
    assert result.similarity == pytest.approx(0.8)


def test_content_counts_represented_families() -> None:
    assert score(COLINEAR[:3], total=4).content == 0.75


def test_similarity_is_mean_identity_times_coverage() -> None:
    pairings = [
        pairing(0, 0, identity=90.0, coverage=0.5),
        pairing(1000, 1000, identity=60.0, coverage=1.0),
    ]
    assert score(pairings).similarity == pytest.approx((0.45 + 0.6) / 2)


def test_reversed_region_keeps_order_and_orientation() -> None:
    reversed_copy = [pairing(n * 1000, 60_000 - n * 1000, target_strand=-1) for n in range(4)]
    result = score(reversed_copy)
    assert (result.order, result.orientation) == (1.0, 1.0)


def test_one_rearranged_family_lowers_order() -> None:
    swapped = [pairing(0, 0), pairing(1000, 3000), pairing(2000, 1000), pairing(3000, 2000)]
    assert score(swapped).order == 0.75


def test_order_takes_the_longer_of_increasing_and_decreasing() -> None:
    mostly_reversed = [pairing(0, 3000), pairing(1000, 2000), pairing(2000, 0), pairing(3000, 1000)]
    assert score(mostly_reversed).order == 0.75


def test_fused_copy_shared_by_two_families_keeps_order() -> None:
    fused = gene("tgt", 1000, length=1800)
    pairings = [
        pairing(0, 0),
        FamilyPairing(member=gene("ref", 1000), copy=fused, identity=80.0, coverage=1.0),
        FamilyPairing(member=gene("ref", 2000), copy=fused, identity=80.0, coverage=1.0),
        pairing(3000, 4000),
    ]
    assert score(pairings).order == 1.0


def test_orientation_counts_agreeing_neighbour_pairs() -> None:
    flipped = [
        pairing(0, 0),
        pairing(1000, 1000, target_strand=-1),
        pairing(2000, 2000),
        pairing(3000, 3000),
    ]
    assert score(flipped).orientation == pytest.approx(1 / 3)


def test_orientation_compares_relative_strands() -> None:
    mixed = [
        pairing(0, 0, ref_strand=-1, target_strand=1),
        pairing(1000, 1000, ref_strand=1, target_strand=-1),
    ]
    assert score(mixed).orientation == 1.0


def test_insertion_lowers_only_spacing() -> None:
    plain, inserted = score(COLINEAR), score(COLINEAR, span=7_800)
    assert inserted.spacing == pytest.approx(3_900 / 7_800)
    assert (inserted.content, inserted.similarity, inserted.order, inserted.orientation) == (
        plain.content,
        plain.similarity,
        plain.order,
        plain.orientation,
    )


def test_compressed_region_lowers_spacing() -> None:
    assert score(COLINEAR, span=1_950).spacing == pytest.approx(1_950 / 3_900)


def test_reference_order_follows_the_members_not_the_input() -> None:
    assert score(list(reversed(COLINEAR))).order == 1.0
    shuffled = [COLINEAR[2], COLINEAR[0], COLINEAR[3], COLINEAR[1]]
    assert score(shuffled).orientation == 1.0


@pytest.mark.parametrize("count", [0, 1])
def test_fewer_than_two_families_give_zero_structure(count: int) -> None:
    result = score(COLINEAR[:count], total=4, span=900)
    assert (result.order, result.orientation, result.spacing) == (0.0, 0.0, 0.0)
    assert result.similarity == (0.0 if count == 0 else pytest.approx(0.8))


@given(
    st.lists(
        st.tuples(st.integers(0, 50), st.sampled_from([1, -1]), st.sampled_from([1, -1])),
        min_size=2,
        max_size=8,
        unique_by=lambda t: t[0],
    )
)
def test_mirroring_the_target_keeps_order_and_orientation(spec) -> None:
    forward = [
        pairing(n * 1000, slot * 1000, ref_strand=ref_strand, target_strand=strand)
        for n, (slot, ref_strand, strand) in enumerate(spec)
    ]
    mirrored = [
        pairing(n * 1000, (100 - slot) * 1000, ref_strand=ref_strand, target_strand=-strand)
        for n, (slot, ref_strand, strand) in enumerate(spec)
    ]
    assert score(forward).order == score(mirrored).order
    assert score(forward).orientation == score(mirrored).orientation


def test_reference_span_is_measured_from_the_first_member() -> None:
    offset = [pairing(10_000 + n * 1000, 50_000 + n * 1000) for n in range(3)]
    assert score(offset, span=2_900).spacing == 1.0
    assert score(offset, span=5_800).spacing == pytest.approx(0.5)
