"""Command-line interface: ``homolocus find`` (spec 13)."""

import dataclasses
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from homolocus import (
    DiamondSearch,
    FamilyGrouping,
    RegionDetection,
    RegionFinder,
    ScoringWeights,
)
from homolocus.cli import FinderSettings, build_parser, finder_settings, main
from sample_world import results as sample_results

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


def genbank(path, record_id, dna, make_record, n_genes=3):
    """A record of ``n_genes`` forward ORFs separated by 60 random bases."""
    parts, features, position = [], [], 60
    parts.append(dna.random(60))
    for _ in range(n_genes):
        orf = dna.orf(150)
        features.append(SeqFeature(FeatureLocation(position, position + len(orf), 1), type="CDS"))
        parts += [orf, dna.random(60)]
        position += len(orf) + 60
    SeqIO.write([make_record(record_id, "".join(parts), features)], path, "genbank")
    return path


@pytest.fixture
def inputs(tmp_path, dna, make_record):
    """Command-line arguments for a readable GenBank reference and target."""
    reference = genbank(tmp_path / "ref.gbk", "ref", dna, make_record)
    target = genbank(tmp_path / "target.gbk", "t1", dna, make_record)
    output = tmp_path / "out"
    return ["find", "--reference", str(reference), "--genome", str(target), "--output", str(output)]


@pytest.fixture
def finder_calls(monkeypatch):
    """Replace the search with the sample results and record each RegionFinder used."""
    calls = []

    def fake_find(finder, genome):
        calls.append((finder, genome))
        return sample_results()

    monkeypatch.setattr(RegionFinder, "find", fake_find)
    return calls


def test_find_reads_inputs_runs_the_finder_and_writes_the_output(inputs, finder_calls) -> None:
    assert main([*inputs, "--min-identity", "40", "--no-family-grouping"]) == 0
    ((finder, genome),) = finder_calls
    assert finder.reference.contig_id == "ref"
    assert len(finder.reference.genes) == 3
    assert finder.search.min_identity == 40.0
    assert finder.family_grouping.family_min_identity is None
    assert len(list(genome.genes())) == 3
    output = inputs[-1]
    for name in ("summary.tsv", "analysis.json", "plots/regions.png"):
        assert (Path(output) / name).is_file()


def test_reference_coordinates_are_one_based(inputs, finder_calls) -> None:
    assert main([*inputs, "--contig", "ref", "--start", "61", "--end", "516"]) == 0
    ((finder, _),) = finder_calls
    assert (finder.reference.start, finder.reference.end) == (60, 516)
    assert len(finder.reference.genes) == 1


def test_start_without_end_is_a_usage_error(inputs, finder_calls, capsys) -> None:
    with pytest.raises(SystemExit) as stopped:
        main([*inputs, "--start", "61"])
    assert stopped.value.code == 2
    assert "--start and --end" in capsys.readouterr().err
    assert finder_calls == []


def test_an_invalid_setting_is_a_usage_error(inputs, finder_calls, capsys) -> None:
    with pytest.raises(SystemExit) as stopped:
        main([*inputs, "--min-matches", "1"])
    assert stopped.value.code == 2
    assert "min_matches 1 must be at least 2" in capsys.readouterr().err
    assert finder_calls == []


def test_an_unreadable_input_fails_with_a_one_line_error(inputs, tmp_path, capsys) -> None:
    missing = str(tmp_path / "missing.gbk")
    assert main([*inputs, "--reference", missing]) == 1
    err = capsys.readouterr().err
    assert err.startswith("homolocus: error: ")
    assert "missing.gbk" in err
    assert "Traceback" not in err
    assert len(err.splitlines()) == 1


def test_an_invalid_region_fails_with_a_one_line_error(inputs, capsys) -> None:
    assert main([*inputs, "--start", "500", "--end", "100"]) == 1
    assert capsys.readouterr().err.startswith("homolocus: error: ")


def test_missing_diamond_fails_with_a_one_line_error(inputs, monkeypatch, capsys) -> None:
    monkeypatch.setenv("PATH", "")
    assert main(inputs) == 1
    err = capsys.readouterr().err
    assert err.startswith("homolocus: error: ")
    assert "DIAMOND" in err
