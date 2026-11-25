"""Search configuration (spec 6.2)."""

from typing import Any

import pytest

from homolocus import DiamondSearch, FamilyGrouping, PlotOptions, RegionDetection, ScoringWeights


def test_defaults_match_the_specification() -> None:
    config = DiamondSearch()
    assert config.evalue == 1e-5
    assert config.hsp_evalue == 1e-3
    assert config.min_identity == 30.0
    assert config.sensitivity == "sensitive"
    assert config.max_target_seqs == 100
    assert config.max_hsps == 0
    assert config.threads is None


@pytest.mark.parametrize(
    "overrides",
    [
        {"evalue": 0.0, "hsp_evalue": 0.0},
        {"evalue": 1e-3, "hsp_evalue": 1e-3},
        {"min_identity": 0.0},
        {"min_identity": 100.0},
        {"max_target_seqs": 0},
        {"max_hsps": 1},
        {"threads": 1},
        {"sensitivity": "default"},
        {"sensitivity": "ultra-sensitive"},
    ],
)
def test_boundary_values_are_allowed(overrides: dict[str, Any]) -> None:
    config = DiamondSearch(**overrides)
    assert all(getattr(config, key) == value for key, value in overrides.items())


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"evalue": -1e-9}, "evalue"),
        ({"evalue": 1e-2}, "hsp_evalue"),
        ({"min_identity": -0.1}, "min_identity"),
        ({"min_identity": 100.1}, "min_identity"),
        ({"sensitivity": "turbo"}, "sensitivity"),
        ({"max_target_seqs": -1}, "max_target_seqs"),
        ({"max_hsps": -1}, "max_hsps"),
        ({"threads": 0}, "threads"),
    ],
)
def test_invalid_settings_are_rejected(overrides: dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        DiamondSearch(**overrides)


def test_family_grouping_defaults() -> None:
    config = FamilyGrouping()
    assert (config.family_min_identity, config.family_min_coverage) == (90.0, 0.8)


@pytest.mark.parametrize(
    ("identity", "coverage"),
    [(None, 0.8), (0.0, 0.0), (100.0, 1.0)],
)
def test_family_grouping_boundaries_are_allowed(identity: float | None, coverage: float) -> None:
    config = FamilyGrouping(family_min_identity=identity, family_min_coverage=coverage)
    assert (config.family_min_identity, config.family_min_coverage) == (identity, coverage)


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"family_min_identity": -0.1}, "family_min_identity"),
        ({"family_min_identity": 100.1}, "family_min_identity"),
        ({"family_min_coverage": -0.01}, "family_min_coverage"),
        ({"family_min_coverage": 1.01}, "family_min_coverage"),
    ],
)
def test_invalid_family_grouping_is_rejected(overrides: dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        FamilyGrouping(**overrides)


def test_grouping_is_disabled_when_identity_is_none() -> None:
    assert not FamilyGrouping(family_min_identity=None).enabled
    assert FamilyGrouping().enabled


def test_region_detection_defaults() -> None:
    config = RegionDetection()
    assert (config.min_matches, config.min_reference_coverage) == (3, 0.3)
    assert (config.min_reference_gene_coverage, config.min_event_coverage) == (0.5, 0.5)
    assert (config.max_intervening_genes, config.max_inter_match_bp_gap) == (10, 20_000)
    assert (config.max_region_bp, config.max_candidate_overlap) == (50_000, 0.5)
    assert (config.padding_genes, config.padding_bp, config.edge_margin_bp) == (0, 0, 500)
    assert (
        config.min_split_contribution,
        config.max_fragment_overlap,
        config.min_fusion_part_coverage,
    ) == (0.2, 0.2, 0.2)


@pytest.mark.parametrize(
    "overrides",
    [
        {"min_matches": 2},
        {"min_reference_coverage": 0.0, "min_reference_gene_coverage": 1.0},
        {"min_event_coverage": 0.0, "max_candidate_overlap": 1.0},
        {"max_intervening_genes": 0, "max_inter_match_bp_gap": 0},
        {"max_region_bp": 1, "padding_genes": 0, "padding_bp": 0, "edge_margin_bp": 0},
        {"min_split_contribution": 0.0, "max_fragment_overlap": 1.0},
        {"min_fusion_part_coverage": 1.0, "max_fragment_overlap": 0.0},
    ],
)
def test_region_detection_boundaries_are_allowed(overrides: dict[str, Any]) -> None:
    config = RegionDetection(**overrides)
    assert all(getattr(config, key) == value for key, value in overrides.items())


@pytest.mark.parametrize(
    "field",
    [
        "min_reference_coverage",
        "min_reference_gene_coverage",
        "min_event_coverage",
        "max_candidate_overlap",
        "min_split_contribution",
        "max_fragment_overlap",
        "min_fusion_part_coverage",
    ],
)
@pytest.mark.parametrize("value", [-0.01, 1.01])
def test_region_detection_fractions_are_checked(field: str, value: float) -> None:
    with pytest.raises(ValueError, match=field):
        RegionDetection(**{field: value})


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("min_matches", 1),
        ("max_intervening_genes", -1),
        ("max_inter_match_bp_gap", -1),
        ("max_region_bp", 0),
        ("padding_genes", -1),
        ("padding_bp", -1),
        ("edge_margin_bp", -1),
    ],
)
def test_region_detection_counts_are_checked(field: str, value: int) -> None:
    with pytest.raises(ValueError, match=field):
        RegionDetection(**{field: value})


def test_scoring_weight_defaults() -> None:
    weights = ScoringWeights()
    assert (weights.content, weights.similarity, weights.order) == (0.40, 0.30, 0.15)
    assert (weights.orientation, weights.spacing) == (0.10, 0.05)


def test_weights_are_normalised_to_sum_to_one() -> None:
    normalized = ScoringWeights(
        content=2, similarity=1, order=1, orientation=0, spacing=0
    ).normalized()
    assert normalized == ScoringWeights(
        content=0.5, similarity=0.25, order=0.25, orientation=0.0, spacing=0.0
    )
    assert sum(vars(ScoringWeights().normalized()).values()) == pytest.approx(1.0)


@pytest.mark.parametrize("field", ["content", "similarity", "order", "orientation", "spacing"])
def test_negative_weights_are_rejected(field: str) -> None:
    with pytest.raises(ValueError, match=field):
        ScoringWeights(**{field: -0.01})


def test_all_zero_weights_are_rejected() -> None:
    with pytest.raises(ValueError, match="at least one weight"):
        ScoringWeights(content=0, similarity=0, order=0, orientation=0, spacing=0)


def test_plot_option_defaults() -> None:
    assert PlotOptions() == PlotOptions(candidates=3, lines_per_page=30)


@pytest.mark.parametrize(("candidates", "lines_per_page"), [(0, 2), (1, 2), (50, 500)])
def test_plot_option_boundaries_are_allowed(candidates: int, lines_per_page: int) -> None:
    PlotOptions(candidates=candidates, lines_per_page=lines_per_page)


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"candidates": -1}, "candidates -1 must be 0 .all. or more"),
        ({"lines_per_page": 1}, "lines_per_page 1 must be at least 2"),
    ],
)
def test_invalid_plot_options_are_rejected(overrides: dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        PlotOptions(**overrides)
