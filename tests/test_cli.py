"""Command-line interface: ``homolocus find`` (spec 13)."""

import dataclasses

import pytest

from homolocus import DiamondSearch, FamilyGrouping, RegionDetection, ScoringWeights
from homolocus.cli import FinderSettings, build_parser, finder_settings

REQUIRED = ["find", "--reference", "ref.gbk", "--genome", "target.fna", "--output", "out"]
CONFIGS = (
    (DiamondSearch, ""),
    (FamilyGrouping, ""),
    (RegionDetection, ""),
    (ScoringWeights, "weight-"),
)


def settings(*options: str) -> FinderSettings:
    return finder_settings(build_parser().parse_args([*REQUIRED, *options]))


def find_help(capsys) -> str:
    with pytest.raises(SystemExit) as stopped:
        build_parser().parse_args(["find", "--help"])
    assert stopped.value.code == 0
    return str(capsys.readouterr().out)


def test_only_required_options_give_the_default_settings() -> None:
    assert settings() == FinderSettings(
        DiamondSearch(), FamilyGrouping(), RegionDetection(), ScoringWeights()
    )


def test_every_config_field_has_a_documented_flag(capsys) -> None:
    text = find_help(capsys)
    for config, prefix in CONFIGS:
        for field in dataclasses.fields(config):
            assert f"--{prefix}{field.name.replace('_', '-')} " in text
    assert "--no-family-grouping" in text


def test_help_shows_the_dataclass_defaults(capsys) -> None:
    text = " ".join(find_help(capsys).split())
    assert "(default: 30.0)" in text
    assert "(default: sensitive)" in text
    assert "(default: 0.15)" in text
    assert "(default: None)" not in text


def test_typed_flags_reach_their_settings() -> None:
    chosen = settings(
        "--min-identity",
        "40",
        "--sensitivity",
        "fast",
        "--threads",
        "2",
        "--max-target-seqs",
        "0",
        "--padding-bp",
        "500",
        "--min-reference-coverage",
        "0.6",
        "--family-min-coverage",
        "0.7",
        "--weight-order",
        "0.5",
    )
    assert chosen.search == DiamondSearch(
        min_identity=40.0, sensitivity="fast", threads=2, max_target_seqs=0
    )
    assert chosen.detection == RegionDetection(padding_bp=500, min_reference_coverage=0.6)
    assert chosen.family_grouping == FamilyGrouping(family_min_coverage=0.7)
    assert chosen.scoring == ScoringWeights(order=0.5)


def test_family_grouping_can_be_disabled_or_retuned() -> None:
    assert settings("--no-family-grouping").family_grouping.family_min_identity is None
    assert settings("--family-min-identity", "95").family_grouping.family_min_identity == 95.0


@pytest.mark.parametrize(
    "options",
    [
        ["--sensitivity", "bogus"],
        ["--min-identity", "high"],
        ["--padding-bp", "1.5"],
        ["--no-family-grouping", "--family-min-identity", "95"],
    ],
)
def test_malformed_flags_are_usage_errors(options, capsys) -> None:
    with pytest.raises(SystemExit) as stopped:
        build_parser().parse_args([*REQUIRED, *options])
    assert stopped.value.code == 2
    assert "error:" in capsys.readouterr().err


@pytest.mark.parametrize(
    "argv",
    [
        [],
        ["find", "--genome", "t.fna", "--output", "out"],
        ["find", "--reference", "r.gbk", "--output", "out"],
        ["find", "--reference", "r.gbk", "--genome", "t.fna"],
    ],
)
def test_a_command_and_the_input_and_output_paths_are_required(argv, capsys) -> None:
    with pytest.raises(SystemExit) as stopped:
        build_parser().parse_args(argv)
    assert stopped.value.code == 2
    assert "required" in capsys.readouterr().err


def test_top_level_help_describes_the_tool(capsys) -> None:
    with pytest.raises(SystemExit):
        build_parser().parse_args(["--help"])
    assert "homolocus find" in capsys.readouterr().out
