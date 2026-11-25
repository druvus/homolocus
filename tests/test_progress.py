"""Progress through the ``homolocus`` logger, -q/-v and --version (spec 13)."""

import logging
from dataclasses import replace

import pytest

from fake_tools import fake_diamond
from homolocus import RegionFinder
from homolocus.annotation.genome_file import read_genome
from homolocus.cli import main
from homolocus.search.diamond import DiamondNotFoundError, diamond_version
from sample_world import GENOME, REFERENCE
from sample_world import results as sample_results


def messages(caplog: pytest.LogCaptureFixture) -> list[str]:
    return [record.getMessage() for record in caplog.records if record.name.startswith("homolocus")]


def test_find_all_reports_the_reference_families_and_each_genome(searches, caplog) -> None:
    genomes = [replace(GENOME, name="a"), replace(GENOME, name="b", display_name="E. coli")]
    with caplog.at_level(logging.INFO, logger="homolocus"):
        RegionFinder(REFERENCE).find_all(genomes)
    assert messages(caplog) == [
        "reference ref:1-4500, genes: 4",
        "reference families: 4",
        "searching genome 1/2: a",
        "a: 0 candidates",
        "searching genome 2/2: E. coli",
        "E. coli: 0 candidates",
    ]


def test_find_reports_its_single_genome(searches, caplog) -> None:
    with caplog.at_level(logging.INFO, logger="homolocus"):
        RegionFinder(REFERENCE).find(GENOME)
    assert messages(caplog)[2:] == ["searching genome 1/1: target", "target: 0 candidates"]


def test_reading_a_genome_reports_its_contigs_and_genes(tmp_path, dna, caplog) -> None:
    path = tmp_path / "strain_A.fna"
    path.write_text(f">c1\n{dna.random(500)}{dna.orf(300)}{dna.random(500)}\n", encoding="utf-8")
    with caplog.at_level(logging.INFO, logger="homolocus"):
        genome = read_genome(path)
    n_genes = len(list(genome.genes()))
    assert messages(caplog) == [f"read genome strain_A, contigs: 1, genes: {n_genes}"]


def test_diamond_command_lines_are_reported_at_debug(tmp_path, caplog) -> None:
    executable = fake_diamond(tmp_path, 'echo "diamond version 2.1.0"\n')
    with caplog.at_level(logging.DEBUG, logger="homolocus"):
        diamond_version(executable)
    assert messages(caplog) == [f"running: {executable} version"]


def test_the_package_adds_no_log_handlers(inputs, finder_calls) -> None:
    assert main(inputs) == 0
    assert logging.getLogger("homolocus").handlers == []


def debug_finder(monkeypatch):
    def fake_find(finder, genome):
        logging.getLogger("homolocus.search.diamond").debug("running: diamond blastp")
        return sample_results()

    monkeypatch.setattr(RegionFinder, "find", fake_find)


def test_find_shows_progress_on_stderr_and_keeps_stdout(inputs, finder_calls, capsys) -> None:
    assert main(inputs) == 0
    captured = capsys.readouterr()
    assert "homolocus: read genome target, contigs: 1, genes: 3\n" in captured.err
    assert f"homolocus: writing results to {inputs[-1]}\n" in captured.err
    assert "homolocus:" not in captured.out


@pytest.mark.parametrize("flag", ["-q", "--quiet"])
def test_quiet_shows_no_progress(inputs, finder_calls, capsys, flag) -> None:
    assert main([*inputs, flag]) == 0
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize(("flags", "shown"), [((), False), (("-v",), True), (("--verbose",), True)])
def test_verbose_adds_diamond_command_lines(inputs, monkeypatch, capsys, flags, shown) -> None:
    debug_finder(monkeypatch)
    assert main([*inputs, *flags]) == 0
    assert ("homolocus: running: diamond blastp\n" in capsys.readouterr().err) is shown


def test_quiet_and_verbose_cannot_be_combined(inputs, capsys) -> None:
    with pytest.raises(SystemExit) as stop:
        main([*inputs, "-q", "-v"])
    assert stop.value.code == 2


def test_version_names_homolocus_and_diamond(monkeypatch, capsys) -> None:
    monkeypatch.setattr("homolocus.cli.diamond_version", lambda: "2.2.9")
    with pytest.raises(SystemExit) as stop:
        main(["--version"])
    assert stop.value.code == 0
    out = capsys.readouterr().out
    assert out.startswith("homolocus ")
    assert out.endswith(" (DIAMOND 2.2.9)\n")


def test_version_without_diamond_still_succeeds(monkeypatch, capsys) -> None:
    def missing():
        raise DiamondNotFoundError("not found")

    monkeypatch.setattr("homolocus.cli.diamond_version", missing)
    with pytest.raises(SystemExit) as stop:
        main(["--version"])
    assert stop.value.code == 0
    assert capsys.readouterr().out.endswith(" (DIAMOND not found)\n")


def test_one_candidate_is_reported_in_the_singular(searches, monkeypatch, caplog) -> None:
    (first, *_) = sample_results().candidates
    monkeypatch.setattr("homolocus.finder.detect", lambda context, genome, matches: (first,))
    with caplog.at_level(logging.INFO, logger="homolocus"):
        RegionFinder(REFERENCE).find(GENOME)
    assert messages(caplog)[-1] == "target: 1 candidate"
