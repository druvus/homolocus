"""Genome presentation names and order from a mapping file (spec 5.3, 12, 13)."""

import warnings
from dataclasses import replace

import pytest

from homolocus.genome_names import apply_genome_names, read_genome_names
from sample_world import GENOME


def mapping(tmp_path, text: str):
    path = tmp_path / "names.tsv"
    path.write_text(text, encoding="utf-8")
    return path


def test_lines_map_file_names_to_presentation_names_in_file_order(tmp_path) -> None:
    text = (
        "# file\tname\n"
        "sample_031.gbk\tE. coli K-12 MG1655\n"
        "\n"
        "sample_007.gbk\tE. coli O157:H7 Sakai\n"
        "sample_112.gbk\tE. coli CFT073 / UPEC\n"
    )
    names = read_genome_names(mapping(tmp_path, text))
    assert list(names.items()) == [
        ("sample_031.gbk", "E. coli K-12 MG1655"),
        ("sample_007.gbk", "E. coli O157:H7 Sakai"),
        ("sample_112.gbk", "E. coli CFT073 / UPEC"),
    ]


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("sample_031.gbk E. coli K-12\n", "line 1: expected a file name, a tab and a name"),
        ("a.gbk\tA\tB\n", "line 1: expected a file name, a tab and a name"),
        ("a.gbk\t\n", "line 1: expected a file name, a tab and a name"),
        ("\tA\n", "line 1: expected a file name, a tab and a name"),
        ("a.gbk\tA\n# x\na.gbk\tB\n", "line 3: file a.gbk is listed twice"),
        ("a.gbk\tA\nb.gbk\tA\n", "line 2: name 'A' is used twice"),
    ],
)
def test_malformed_or_duplicated_lines_are_rejected_with_their_line(
    tmp_path, text: str, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        read_genome_names(mapping(tmp_path, text))


def genome(file: str):
    return replace(GENOME, name=file.split(".", maxsplit=1)[0], source=f"data/run_1/{file}")


def test_listed_genomes_come_first_in_mapping_order_with_their_names() -> None:
    a, b, c = genome("a.gbk"), genome("b.gbk"), genome("c.fna")
    ordered = apply_genome_names([a, b, c], {"c.fna": "Gamma", "a.gbk": "Alpha"})
    assert [(g.name, g.label) for g in ordered] == [("c", "Gamma"), ("a", "Alpha"), ("b", "b")]
    assert ordered[0].source == "data/run_1/c.fna"


def test_lines_that_match_no_genome_warn() -> None:
    with pytest.warns(UserWarning, match="names file lists missing.gbk, which is not a --genome"):
        ordered = apply_genome_names([genome("a.gbk")], {"missing.gbk": "X", "a.gbk": "A"})
    assert [g.label for g in ordered] == ["A"]


def test_an_empty_mapping_changes_nothing() -> None:
    genomes = [genome("a.gbk"), genome("b.gbk")]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert apply_genome_names(genomes, {}) == genomes
