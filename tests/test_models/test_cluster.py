"""Cluster result model (spec 5.8, 8, 8.1)."""

from typing import Any

import pytest

from homolocus import (
    Alignment,
    Cluster,
    FamilyEvidence,
    Gene,
    GeneMatch,
    ScoreComponents,
    ScoringWeights,
)


def gene(start: int, end: int, partial: bool = False, contig_id: str = "t1") -> Gene:
    return Gene(contig_id, start, end, 1, "A" * (end - start), "M", partial=partial)


LEFT_PAD, EVENT_A, INSERT, EVENT_B, RIGHT_PAD = (
    gene(0, 300),
    gene(1000, 1300),
    gene(1500, 1800),
    gene(2000, 2300),
    gene(3000, 3300),
)
REF_1, REF_2, REF_3 = "ref|1-300|+", "ref|401-700|+", "ref|801-1100|+"


def gene_match(reference_id: str, target: Gene) -> GeneMatch:
    hsp = Alignment(reference_id, target.id, 80.0, 100, 0, 100, 0, 100, 100, 100, 1e-40, 200.0)
    return GeneMatch(reference_id, target.id, (hsp,))


def evidence(
    family_id: str,
    member: str | None,
    copies: tuple[str, ...],
    coverage: float,
    *,
    represented: bool,
) -> FamilyEvidence:
    return FamilyEvidence(
        family_id=family_id,
        coverage=coverage,
        represented=represented,
        best_member_id=member,
        best_copy_id=copies[0] if copies else None,
        copies=copies,
    )


FAMILIES = (
    evidence("F001", REF_1, (EVENT_A.id,), 0.9, represented=True),
    evidence("F002", REF_2, (EVENT_B.id, INSERT.id), 0.6, represented=True),
    evidence("F003", None, (), 0.0, represented=False),
)
COMPONENTS = ScoreComponents(content=0.5, similarity=1.0, order=1.0, orientation=0.0, spacing=0.5)


def make_cluster(**overrides: Any) -> Cluster:
    fields: dict[str, Any] = {
        "contig_id": "t1",
        "start": 1000,
        "end": 2300,
        "display_start": 0,
        "display_end": 3300,
        "genes": (LEFT_PAD, EVENT_A, INSERT, EVENT_B, RIGHT_PAD),
        "event_ids": (EVENT_A.id, EVENT_B.id),
        "matches": (
            gene_match(REF_1, EVENT_A),
            gene_match(REF_2, EVENT_B),
            gene_match(REF_2, INSERT),
        ),
        "families": FAMILIES,
        "components": COMPONENTS,
        "weights": ScoringWeights(),
        "touches_contig_edge": False,
    }
    fields.update(overrides)
    return Cluster(**fields)


# ScoreComponents


@pytest.mark.parametrize("field", ["content", "similarity", "order", "orientation", "spacing"])
@pytest.mark.parametrize("value", [-0.01, 1.01])
def test_components_must_be_fractions(field: str, value: float) -> None:
    values = {"content": 0.5, "similarity": 0.5, "order": 0.5, "orientation": 0.5, "spacing": 0.5}
    values[field] = value
    with pytest.raises(ValueError, match=field):
        ScoreComponents(**values)


# FamilyEvidence


def test_family_evidence_counts_copies() -> None:
    assert FAMILIES[1].copy_count == 2
    assert FAMILIES[2].copy_count == 0


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"family_id": ""}, "family ID"),
        ({"coverage": 1.1}, "coverage"),
        ({"coverage": -0.1}, "coverage"),
        ({"best_copy_id": "x|1-3|+"}, "best copy"),
        ({"best_member_id": None}, "best member"),
        ({"best_copy_id": None}, "represented"),
        ({"copies": (EVENT_A.id, EVENT_A.id)}, "duplicate"),
    ],
)
def test_illegal_family_evidence_is_rejected(overrides: dict[str, Any], message: str) -> None:
    fields = {
        "family_id": "F001",
        "coverage": 0.9,
        "represented": True,
        "best_member_id": REF_1,
        "best_copy_id": EVENT_A.id,
        "copies": (EVENT_A.id,),
    }
    fields.update(overrides)
    with pytest.raises(ValueError, match=message):
        FamilyEvidence(**fields)


def test_unrepresented_family_without_evidence_has_zero_coverage() -> None:
    with pytest.raises(ValueError, match="coverage"):
        evidence("F003", None, (), 0.2, represented=False)


# Cluster


def test_roles_follow_the_events() -> None:
    cluster = make_cluster()
    assert cluster.role(EVENT_A.id) == "matched"
    assert cluster.role(INSERT.id) == "additional"
    assert cluster.role(LEFT_PAD.id) == "additional"
    with pytest.raises(KeyError):
        cluster.role("t1|9-10|+")


def test_evidence_genes_lie_fully_inside_the_unpadded_span() -> None:
    assert [g.id for g in make_cluster().evidence_genes] == [EVENT_A.id, INSERT.id, EVENT_B.id]


def test_coverage_metrics() -> None:
    cluster = make_cluster()
    assert [f.family_id for f in cluster.represented_families] == ["F001", "F002"]
    assert [f.family_id for f in cluster.missing_families] == ["F003"]
    assert cluster.reference_coverage == pytest.approx(2 / 3)
    assert cluster.mean_protein_coverage == pytest.approx(0.75)
    assert cluster.span == 1300


def test_no_represented_family_gives_zero_protein_coverage() -> None:
    family = (evidence("F001", None, (), 0.0, represented=False),)
    cluster = make_cluster(families=family, matches=())
    assert (cluster.reference_coverage, cluster.mean_protein_coverage) == (0.0, 0.0)


def test_score_is_the_normalised_weighted_sum() -> None:
    weights = ScoringWeights(content=2, similarity=2, order=0, orientation=0, spacing=0)
    assert make_cluster(weights=weights).score == pytest.approx(0.75)
    assert make_cluster().score == pytest.approx(0.4 * 0.5 + 0.3 + 0.15 + 0.05 * 0.5)


def test_partial_evidence_gene_is_reported() -> None:
    partial_insert = gene(1500, 1800, partial=True)
    genes = (LEFT_PAD, EVENT_A, partial_insert, EVENT_B, RIGHT_PAD)
    assert make_cluster(genes=genes).contains_partial_genes
    assert not make_cluster().contains_partial_genes
    padded_partial = (gene(0, 300, partial=True), EVENT_A, INSERT, EVENT_B, RIGHT_PAD)
    assert not make_cluster(genes=padded_partial).contains_partial_genes


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"start": 2300}, "interval"),
        ({"display_start": 1001}, "display"),
        ({"display_end": 2299}, "display"),
        ({"genes": (EVENT_A, EVENT_B, gene(0, 300, contig_id="t2"))}, "contig"),
        ({"event_ids": ()}, "at least one event"),
        ({"event_ids": (LEFT_PAD.id,)}, "event"),
        ({"matches": (gene_match(REF_1, RIGHT_PAD),)}, "match"),
        ({"families": ()}, "at least one family"),
        ({"families": (FAMILIES[0], FAMILIES[0])}, "duplicate family"),
        ({"families": (evidence("F001", REF_1, (RIGHT_PAD.id,), 0.9, represented=True),)}, "copy"),
    ],
)
def test_illegal_clusters_are_rejected(overrides: dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        make_cluster(**overrides)


@pytest.mark.parametrize("coverage", [0.0, 1.0])
def test_coverage_bounds_are_allowed(coverage: float) -> None:
    assert evidence("F001", REF_1, (EVENT_A.id,), coverage, represented=True).coverage == coverage


def test_cluster_without_padding_starting_at_zero() -> None:
    first, second = gene(0, 300), gene(500, 800)
    families = (evidence("F001", REF_1, (first.id, second.id), 1.0, represented=True),)
    cluster = make_cluster(
        start=0,
        end=800,
        display_start=0,
        display_end=800,
        genes=(first, second),
        event_ids=(first.id, second.id),
        matches=(gene_match(REF_1, first), gene_match(REF_1, second)),
        families=families,
    )
    assert (cluster.display_start, cluster.display_end) == (cluster.start, cluster.end)
    assert [g.id for g in cluster.evidence_genes] == [first.id, second.id]


def test_score_uses_every_weighted_component() -> None:
    parts = ScoreComponents(content=0.2, similarity=0.4, order=0.6, orientation=0.8, spacing=0.5)
    expected = 0.40 * 0.2 + 0.30 * 0.4 + 0.15 * 0.6 + 0.10 * 0.8 + 0.05 * 0.5
    assert make_cluster(components=parts).score == pytest.approx(expected)
