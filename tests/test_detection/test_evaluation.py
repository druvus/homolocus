"""Candidate evaluation: coverage, families, ranking (spec 7.7, 7.9, 8)."""

import pytest

from homolocus import Alignment, GeneFamily, GeneMatch
from homolocus.detection.evaluation import family_evidence

R1, R2, R3 = "ref|1-300|+", "ref|401-700|+", "ref|801-1100|-"
T1, T2, T3 = "t|1-300|+", "t|401-700|+", "t|801-1100|+"
TARGET_ORDER = {T1: 0, T2: 1, T3: 2}


def hit(
    reference: str,
    target: str,
    *spans: tuple[int, int],
    bit_score: float = 100.0,
    query_length: int = 100,
) -> GeneMatch:
    hsps = tuple(
        Alignment(
            reference,
            target,
            80.0,
            end - start,
            start,
            end,
            0,
            end - start,
            query_length,
            200,
            1e-30,
            bit_score,
        )
        for start, end in spans
    )
    return GeneMatch(reference, target, hsps)


def evidence_for(family: GeneFamily, *matches: GeneMatch, threshold: float = 0.5):
    return family_evidence(
        family, matches, target_order=TARGET_ORDER, min_reference_gene_coverage=threshold
    )


def test_split_fragments_combine_into_one_coverage() -> None:
    family = GeneFamily("F001", (R1,))
    result = evidence_for(family, hit(R1, T1, (0, 45)), hit(R1, T2, (45, 100)))
    assert result.coverage == 1.0
    assert result.represented
    assert result.copies == (T1, T2)


def test_overlapping_hsps_are_counted_once() -> None:
    family = GeneFamily("F001", (R1,))
    result = evidence_for(family, hit(R1, T1, (0, 40), (20, 60)), hit(R1, T2, (50, 70)))
    assert result.coverage == pytest.approx(0.7)


@pytest.mark.parametrize(("end", "represented"), [(50, True), (49, False)])
def test_represented_threshold_is_inclusive(end: int, represented: bool) -> None:
    result = evidence_for(GeneFamily("F001", (R1,)), hit(R1, T1, (0, end)))
    assert result.represented is represented
    assert result.best_copy_id == T1


def test_best_member_has_the_highest_coverage() -> None:
    family = GeneFamily("F001", (R1, R2))
    result = evidence_for(family, hit(R1, T1, (0, 30)), hit(R2, T2, (0, 90)))
    assert (result.best_member_id, result.coverage) == (R2, 0.9)
    assert result.best_copy_id == T2
    assert result.copies == (T1, T2)


def test_tied_members_resolve_by_genomic_order() -> None:
    family = GeneFamily("F001", (R1, R2))
    result = evidence_for(family, hit(R2, T1, (0, 80)), hit(R1, T2, (0, 80)))
    assert result.best_member_id == R1


def test_best_copy_has_the_highest_bit_score_to_the_best_member() -> None:
    family = GeneFamily("F001", (R1,))
    result = evidence_for(
        family,
        hit(R1, T1, (0, 60), bit_score=150.0),
        hit(R1, T3, (0, 90), bit_score=400.0),
        hit(R1, T2, (0, 90), bit_score=200.0),
    )
    assert result.best_copy_id == T3
    assert result.copies == (T1, T2, T3)
    assert result.copy_count == 3


def test_tied_copies_resolve_by_target_order() -> None:
    family = GeneFamily("F001", (R1,))
    result = evidence_for(family, hit(R1, T3, (0, 90)), hit(R1, T2, (0, 90)))
    assert result.best_copy_id == T2


def test_copies_of_other_members_count_but_do_not_pair() -> None:
    family = GeneFamily("F001", (R1, R2))
    result = evidence_for(
        family, hit(R1, T1, (0, 90), bit_score=50.0), hit(R2, T3, (0, 40), bit_score=900.0)
    )
    assert (result.best_member_id, result.best_copy_id) == (R1, T1)
    assert result.copies == (T1, T3)


def test_matches_to_other_families_are_ignored() -> None:
    result = evidence_for(GeneFamily("F001", (R1,)), hit(R3, T1, (0, 90)))
    assert (result.coverage, result.represented) == (0.0, False)
    assert (result.best_member_id, result.best_copy_id, result.copies) == (None, None, ())


def test_family_without_matches_is_missing() -> None:
    result = evidence_for(GeneFamily("F002", (R2,)))
    assert result.family_id == "F002"
    assert not result.represented


def test_matches_may_be_given_as_a_one_shot_iterator() -> None:
    family = GeneFamily("F001", (R1, R2))
    matches = iter([hit(R1, T1, (0, 30)), hit(R2, T2, (0, 90))])
    result = family_evidence(
        family, matches, target_order=TARGET_ORDER, min_reference_gene_coverage=0.5
    )
    assert (result.best_member_id, result.copies) == (R2, (T1, T2))
