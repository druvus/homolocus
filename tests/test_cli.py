"""Command-line interface: ``homolocus find`` (spec 13)."""

import dataclasses
from dataclasses import replace
from pathlib import Path

import pytest

from genome_files import fasta, genbank
from homolocus import (
    DiamondSearch,
    FamilyGrouping,
    GenomeSetResults,
    RegionDetection,
    RegionFinder,
    ScoringWeights,
)
from homolocus.cli import FinderSettings, build_parser, finder_settings, main
from sample_world import GENOME as SAMPLE_GENOME
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


def test_find_prints_one_line_per_candidate_then_the_output(inputs, finder_calls, capsys) -> None:
    assert main(inputs) == 0
    assert capsys.readouterr().out.splitlines() == [
        "#1 candidate_001  t1:10001-14500  score 0.827  coverage 0.75",
        "#2 candidate_002  t1:100001-102100  score 0.740  coverage 0.50",
        f"results written to {inputs[-1]}",
    ]


def test_candidate_flags_end_their_line(inputs, monkeypatch, capsys) -> None:
    found = sample_results()
    flagged = dataclasses.replace(found.candidates[0], touches_contig_edge=True)
    flagged_results = dataclasses.replace(found, candidates=(flagged,))
    monkeypatch.setattr(RegionFinder, "find", lambda finder, genome: flagged_results)
    assert main(inputs) == 0
    assert capsys.readouterr().out.splitlines()[0].endswith("coverage 0.75  edge")


def test_no_candidates_is_said_plainly(inputs, monkeypatch, capsys) -> None:
    empty = sample_results(candidates=())
    monkeypatch.setattr(RegionFinder, "find", lambda finder, genome: empty)
    assert main(inputs) == 0
    assert capsys.readouterr().out.splitlines() == [
        "no candidates",
        f"results written to {inputs[-1]}",
    ]


def test_errors_print_nothing_on_stdout(inputs, tmp_path, capsys) -> None:
    assert main([*inputs, "--reference", str(tmp_path / "missing.gbk")]) == 1
    assert capsys.readouterr().out == ""


@pytest.mark.diamond
def test_find_end_to_end_with_a_predicted_fasta_target(tmp_path, dna, make_record, capsys) -> None:
    orfs = [dna.orf(150) for _ in range(3)]
    reference = genbank(tmp_path / "ref.gbk", "ref", dna, make_record, orfs)
    target = fasta(tmp_path / "target.fna", "t1", dna, orfs)
    output = tmp_path / "out"
    argv = ["find", "--reference", str(reference), "--genome", str(target), "--output", str(output)]
    assert main([*argv, "--threads", "1"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith("#1 candidate_001  t1:")
    assert "coverage 1.00" in lines[0]
    assert lines[-1] == f"results written to {output}"
    for name in ("summary.tsv", "candidates/candidate_001.tsv", "plots/candidate_001.png"):
        assert (output / name).is_file()


def test_genome_takes_several_paths_or_repeats() -> None:
    argv = ["find", "--reference", "r.gbk", "--genome", "target.fna", "b.fna", "--genome", "c.gbk"]
    args = build_parser().parse_args([*argv, "--output", "out"])
    assert args.genome == ["target.fna", "b.fna", "c.gbk"]


@pytest.fixture
def two_targets(inputs, tmp_path, dna, make_record):
    """Command-line arguments with a second GenBank target, strain_b.gbk."""
    second = genbank(tmp_path / "strain_b.gbk", "t9", dna, make_record)
    genome_at = inputs.index("--genome") + 1
    return [*inputs[: genome_at + 1], str(second), *inputs[genome_at + 1 :]]


@pytest.fixture
def find_all_calls(monkeypatch):
    """Replace the multi-genome search with sample results, one per genome read."""
    calls = []

    def sample_for(genome):
        """The sample results under this genome's name; strain_b finds nothing."""
        named = replace(
            SAMPLE_GENOME, name=genome.name, display_name=genome.display_name, source=genome.source
        )
        if genome.name == "strain_b":
            return sample_results(genome=named, candidates=())
        return sample_results(genome=named)

    def fake_find_all(finder, genomes):
        genomes = list(genomes)
        calls.append(genomes)
        return GenomeSetResults(tuple(sample_for(g) for g in genomes))

    monkeypatch.setattr(RegionFinder, "find_all", fake_find_all)
    return calls


def test_several_genomes_are_named_from_files_and_written_per_genome(
    two_targets, find_all_calls, capsys
) -> None:
    assert main(two_targets) == 0
    ((first, second),) = find_all_calls
    assert (first.name, second.name) == ("target", "strain_b")
    output = Path(two_targets[-1])
    assert (output / "summary.tsv").read_text().startswith("genome\trank\t")
    assert (output / "plots" / "regions.png").is_file()
    assert (output / "target" / "candidates" / "candidate_001.tsv").is_file()
    assert (output / "strain_b" / "summary.tsv").is_file()
    assert capsys.readouterr().out.splitlines() == [
        "target  #1 candidate_001  t1:10001-14500  score 0.827  coverage 0.75",
        "target  #2 candidate_002  t1:100001-102100  score 0.740  coverage 0.50",
        "strain_b: no candidates",
        f"results written to {output}",
    ]


def test_an_unreadable_genome_fails_before_any_search(
    two_targets, find_all_calls, tmp_path, capsys
) -> None:
    missing = str(tmp_path / "missing.fna")
    broken = [missing if a.endswith("strain_b.gbk") else a for a in two_targets]
    assert main(broken) == 1
    assert "missing.fna" in capsys.readouterr().err
    assert find_all_calls == []
    assert not Path(broken[-1]).exists()


def test_genomes_with_the_same_file_name_are_rejected(inputs, tmp_path, dna, make_record, capsys):
    (tmp_path / "other").mkdir()
    clash = genbank(tmp_path / "other" / "target.gbk", "t2", dna, make_record)
    genome_at = inputs.index("--genome") + 1
    argv = [*inputs[: genome_at + 1], str(clash), *inputs[genome_at + 1 :]]
    assert main(argv) == 1
    err = capsys.readouterr().err
    assert err.startswith("homolocus: error: duplicate genome name 'target'")
    assert len(err.splitlines()) == 1


@pytest.mark.diamond
def test_find_end_to_end_with_two_fasta_targets(tmp_path, dna, make_record, capsys) -> None:
    orfs = [dna.orf(150) for _ in range(3)]
    reference = genbank(tmp_path / "ref.gbk", "ref", dna, make_record, orfs)
    first = fasta(tmp_path / "a.fna", "c1", dna, orfs)
    second = fasta(tmp_path / "b.fna", "c1", dna, [dna.orf(150) for _ in range(3)])
    output = tmp_path / "out"
    argv = ["find", "--reference", str(reference), "--genome", str(first), str(second)]
    assert main([*argv, "--output", str(output), "--threads", "1"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith("a  #1 candidate_001  c1:")
    assert lines[1] == "b: no candidates"
    summary = (output / "summary.tsv").read_text().splitlines()
    assert [row.split("\t")[0] for row in summary[1:]] == ["a"]
    assert (output / "a" / "plots" / "candidate_001.png").is_file()
    assert (output / "plots" / "regions.png").is_file()


def with_names(argv, tmp_path, text: str):
    names = tmp_path / "names.tsv"
    names.write_text(text, encoding="utf-8")
    return [*argv, "--names", str(names)]


def test_names_file_sets_order_and_presentation_names(
    two_targets, find_all_calls, tmp_path, capsys
) -> None:
    text = "# order for the figure\nstrain_b.gbk\tE. coli O157:H7\ntarget.gbk\tE. coli K-12\n"
    assert main(with_names(two_targets, tmp_path, text)) == 0
    ((first, second),) = find_all_calls
    assert [(g.name, g.label) for g in (first, second)] == [
        ("strain_b", "E. coli O157:H7"),
        ("target", "E. coli K-12"),
    ]
    output = Path(two_targets[-1])
    assert capsys.readouterr().out.splitlines() == [
        "E. coli O157:H7: no candidates",
        "E. coli K-12  #1 candidate_001  t1:10001-14500  score 0.827  coverage 0.75",
        "E. coli K-12  #2 candidate_002  t1:100001-102100  score 0.740  coverage 0.50",
        f"results written to {output}",
    ]
    assert (output / "target" / "summary.tsv").is_file()
    assert (output / "summary.tsv").read_text().splitlines()[1].startswith("E. coli K-12\t1\t")


def test_unlisted_genomes_follow_and_unmatched_lines_warn(
    two_targets, find_all_calls, tmp_path
) -> None:
    text = "missing.gbk\tE. coli CFT073\nstrain_b.gbk\tE. coli O157:H7\n"
    with pytest.warns(UserWarning, match="names file lists missing.gbk, which is not a --genome"):
        assert main(with_names(two_targets, tmp_path, text)) == 0
    ((first, second),) = find_all_calls
    assert [g.label for g in (first, second)] == ["E. coli O157:H7", "target"]


def test_a_malformed_names_file_fails_before_any_search(
    two_targets, find_all_calls, tmp_path, capsys
) -> None:
    argv = with_names(two_targets, tmp_path, "target.gbk E. coli K-12\n")
    assert main(argv) == 1
    err = capsys.readouterr().err
    assert err.startswith("homolocus: error: names file line 1: expected a file name")
    assert len(err.splitlines()) == 1
    assert find_all_calls == []


@pytest.mark.diamond
def test_names_order_a_real_two_genome_run(tmp_path, dna, make_record, capsys) -> None:
    orfs = [dna.orf(150) for _ in range(3)]
    reference = genbank(tmp_path / "ref.gbk", "ref", dna, make_record, orfs)
    first = fasta(tmp_path / "a.fna", "c1", dna, orfs)
    second = fasta(tmp_path / "b.fna", "c1", dna, orfs)
    argv = ["find", "--reference", str(reference), "--genome", str(first), str(second)]
    argv = with_names(
        [*argv, "--output", str(tmp_path / "out")], tmp_path, "b.fna\tSecond\na.fna\tFirst\n"
    )
    assert main([*argv, "--threads", "1"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith("Second  #1 candidate_001")
    assert lines[1].startswith("First  #1 candidate_001")
