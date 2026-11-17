"""Search configuration (spec 6.2)."""

from typing import Any

import pytest

from homolocus import DiamondSearch, FamilyGrouping


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
