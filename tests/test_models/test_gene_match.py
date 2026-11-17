"""GeneMatch: aggregated evidence for one protein pair (spec 5.7)."""

import pytest

from homolocus import Alignment, GeneMatch
from homolocus.models.gene_match import aggregate

REF = "ref|1-300|+"
TARGET = "tgt|1-330|-"


def hsp(
    query: tuple[int, int],
    subject: tuple[int, int],
    *,
    identity: float = 50.0,
    evalue: float = 1e-10,
    bit_score: float = 100.0,
    query_id: str = REF,
    subject_id: str = TARGET,
    query_length: int = 100,
    subject_length: int = 200,
) -> Alignment:
    length = max(query[1] - query[0], subject[1] - subject[0])
    return Alignment(
        query_id,
        subject_id,
        identity,
        length,
        *query,
        *subject,
        query_length,
        subject_length,
        evalue,
        bit_score,
    )


def match(*alignments: Alignment) -> GeneMatch:
    return GeneMatch(REF, TARGET, alignments)


def test_single_hsp() -> None:
    gene_match = match(hsp((0, 50), (0, 100), identity=40.0, evalue=1e-8, bit_score=80.0))
    assert gene_match.query_coverage == 0.5
    assert gene_match.subject_coverage == 0.5
    assert gene_match.identity == 40.0
    assert gene_match.bit_score == 80.0
    assert gene_match.evalue == 1e-8
    assert (gene_match.query_length, gene_match.subject_length) == (100, 200)


def test_disjoint_hsps_add_up() -> None:
    gene_match = match(hsp((0, 30), (0, 30)), hsp((50, 80), (100, 130)))
    assert gene_match.query_coverage == pytest.approx(0.6)
    assert gene_match.subject_coverage == pytest.approx(0.3)


def test_overlapping_hsps_are_not_double_counted() -> None:
    gene_match = match(hsp((0, 60), (0, 60)), hsp((40, 100), (40, 100)), hsp((10, 20), (10, 20)))
    assert gene_match.query_coverage == 1.0
    assert gene_match.subject_coverage == pytest.approx(0.5)


def test_identity_is_weighted_by_alignment_length() -> None:
    gene_match = match(hsp((0, 30), (0, 30), identity=90.0), hsp((30, 90), (30, 90), identity=30.0))
    assert gene_match.identity == pytest.approx((90 * 30 + 30 * 60) / 90)


def test_bit_scores_are_summed_and_best_evalue_kept() -> None:
    gene_match = match(
        hsp((0, 30), (0, 30), evalue=1e-3, bit_score=30.0),
        hsp((50, 80), (50, 80), evalue=1e-12, bit_score=70.0),
    )
    assert gene_match.bit_score == 100.0
    assert gene_match.evalue == 1e-12


def test_weak_secondary_hsp_is_kept_and_counts_toward_coverage() -> None:
    strong = hsp((0, 50), (0, 50), evalue=1e-30)
    weak = hsp((60, 90), (60, 90), evalue=0.5)
    gene_match = match(strong, weak)
    assert gene_match.alignments == (strong, weak)
    assert gene_match.query_coverage == pytest.approx(0.8)
    assert gene_match.passes(max_evalue=1e-5, min_identity=30)


@pytest.mark.parametrize(
    ("evalue", "identity", "passes"),
    [
        (1e-5, 30.0, True),
        (1e-6, 31.0, True),
        (1.1e-5, 30.0, False),
        (1e-5, 29.9, False),
    ],
)
def test_pair_level_thresholds_are_inclusive(evalue: float, identity: float, passes: bool) -> None:
    gene_match = match(hsp((0, 50), (0, 50), evalue=evalue, identity=identity))
    assert gene_match.passes(max_evalue=1e-5, min_identity=30) is passes


def test_coverage_is_not_a_pair_level_threshold() -> None:
    tiny = match(hsp((0, 1), (0, 1), evalue=1e-50, identity=100.0))
    assert tiny.query_coverage == 0.01
    assert tiny.passes(max_evalue=1e-5, min_identity=30)


@pytest.mark.parametrize(
    ("alignments", "message"),
    [
        ((), "at least one"),
        ((hsp((0, 10), (0, 10), query_id="other|1-3|+"),), "query"),
        ((hsp((0, 10), (0, 10), subject_id="other|1-3|+"),), "subject"),
        ((hsp((0, 10), (0, 10)), hsp((0, 10), (0, 10), query_length=99)), "query length"),
        ((hsp((0, 10), (0, 10)), hsp((0, 10), (0, 10), subject_length=99)), "subject length"),
    ],
)
def test_illegal_matches_are_rejected(alignments: tuple[Alignment, ...], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        GeneMatch(REF, TARGET, alignments)


def test_aggregate_groups_hsps_by_protein_pair_in_first_seen_order() -> None:
    a1 = hsp((0, 10), (0, 10), subject_id="t1|1-3|+")
    b1 = hsp((0, 10), (0, 10), subject_id="t2|1-3|+")
    a2 = hsp((20, 30), (20, 30), subject_id="t1|1-3|+")
    c1 = hsp((0, 10), (0, 10), query_id="r2|1-3|+", subject_id="t1|1-3|+")
    matches = aggregate([a1, b1, a2, c1])
    assert [(m.reference_gene_id, m.target_gene_id) for m in matches] == [
        (REF, "t1|1-3|+"),
        (REF, "t2|1-3|+"),
        ("r2|1-3|+", "t1|1-3|+"),
    ]
    assert matches[0].alignments == (a1, a2)
    assert sum(len(m.alignments) for m in matches) == 4


def test_aggregate_of_nothing_is_empty() -> None:
    assert aggregate([]) == ()
