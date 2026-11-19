"""Several target genomes in one run: find_all and GenomeSetResults (spec 2.2, 12)."""

from dataclasses import replace

import pytest
from matplotlib.figure import Figure

from homolocus import DiamondSearch, GenomeSetResults, ReferenceFamilies, RegionFinder
from homolocus.results import SUMMARY_COLUMNS
from homolocus.visualization.layout import genome_set_tracks, overview_tracks
from homolocus.visualization.region_plot import overview_figure
from sample_world import FAMILIES, GENOME, REFERENCE
from sample_world import results as sample_results
from scenario_world import Slot, World, genome, run


def named(name: str, **overrides):
    """The sample results for a genome called ``name``."""
    return sample_results(genome=replace(GENOME, name=name), **overrides)


def test_results_are_kept_in_input_order_and_found_by_name() -> None:
    first, second = named("a"), named("b", candidates=())
    found = GenomeSetResults((first, second))
    assert list(found) == [first, second]
    assert len(found) == 2
    assert found.names == ("a", "b")
    assert found["b"] is second
    with pytest.raises(KeyError):
        found["c"]


def test_combined_frame_has_the_genome_first_then_rank() -> None:
    frame = GenomeSetResults((named("a"), named("b", candidates=()), named("c"))).to_dataframe()
    assert list(frame.columns) == ["genome", *SUMMARY_COLUMNS]
    assert list(zip(frame["genome"], frame["rank"], strict=True)) == [
        ("a", 1),
        ("a", 2),
        ("c", 1),
        ("c", 2),
    ]


def test_combined_frame_without_candidates_still_has_columns() -> None:
    frame = GenomeSetResults((named("a", candidates=()),)).to_dataframe()
    assert frame.empty
    assert list(frame.columns) == ["genome", *SUMMARY_COLUMNS]


@pytest.mark.parametrize(
    ("names", "message"),
    [
        (("a", "a"), "duplicate genome name 'a'"),
        ((None,), "has no name"),
        (("",), "has no name"),
        (("x/y",), "not usable as a directory"),
        ((".hidden",), "not usable as a directory"),
    ],
)
def test_genome_names_must_be_unique_and_directory_safe(names, message) -> None:
    with pytest.raises(ValueError, match=message):
        GenomeSetResults(tuple(named(n) for n in names))


@pytest.fixture
def searches(monkeypatch):
    """Stand-in for DIAMOND that records each search and finds nothing."""
    calls: list[bool] = []

    def fake_search(index, config, path, *, self_search=False):
        path.write_text("")
        calls.append(self_search)
        return ()

    monkeypatch.setattr("homolocus.finder.search", fake_search)
    monkeypatch.setattr("homolocus.finder.diamond_version", lambda: "2.2.9")
    return calls


def test_find_all_runs_the_self_search_once_and_one_search_per_genome(searches) -> None:
    genomes = [replace(GENOME, name=n) for n in ("a", "b", "c")]
    found = RegionFinder(REFERENCE).find_all(genomes)
    assert found.names == ("a", "b", "c")
    assert searches == [True, False, False, False]
    assert all(results.raw_self == found["a"].raw_self for results in found)


def test_find_all_checks_names_before_any_search(searches) -> None:
    with pytest.raises(ValueError, match="duplicate genome name"):
        RegionFinder(REFERENCE).find_all([GENOME, GENOME])
    assert searches == []


@pytest.mark.diamond
def test_each_genome_matches_a_run_on_its_own() -> None:
    world = World(4)
    reference = world.reference()
    strong = world.contig("t1", [Slot(world.copy(i, 0.95)) for i in range(4)])
    weak = world.contig("t1", [Slot(world.copy(i, 0.7)) for i in (0, 2, 3)])
    unrelated = world.contig("t1", [Slot(world.unrelated()) for _ in range(4)])
    genomes = [
        replace(genome(c), name=n) for c, n in ((strong, "a"), (weak, "b"), (unrelated, "c"))
    ]
    together = RegionFinder(reference, search=DiamondSearch(threads=1)).find_all(genomes)
    for alone in genomes:
        assert together[alone.name] == run(reference, alone)


def three_genomes() -> GenomeSetResults:
    """Genome a and c have two candidates each; b has none."""
    return GenomeSetResults((named("a"), named("b", candidates=()), named("c")))


def test_combined_tracks_are_the_reference_then_each_genome_in_rank_order() -> None:
    tracks = genome_set_tracks(three_genomes())
    single = overview_tracks(named("a"))
    assert tracks[0].title == single[0].title
    assert [t.title.split(" · ")[:2] for t in tracks[1:]] == [
        ["a", "#1 candidate_001"],
        ["a", "#2 candidate_002"],
        ["c", "#1 candidate_001"],
        ["c", "#2 candidate_002"],
    ]
    assert tracks[1].title == "a · " + single[1].title
    assert tracks[1].glyphs == single[1].glyphs


def test_genomes_without_candidates_are_named_on_the_reference_line() -> None:
    assert genome_set_tracks(three_genomes())[0].note == "no candidates: b"
    assert genome_set_tracks(GenomeSetResults((named("a"),)))[0].note == ""
    nothing = GenomeSetResults((named("a", candidates=()), named("b", candidates=())))
    assert genome_set_tracks(nothing)[0].note == "no candidates: a, b"


def test_combined_overview_has_one_axis_per_line_and_one_legend() -> None:
    found = three_genomes()
    figure = found.plot()
    assert isinstance(figure, Figure)
    assert figure.canvas.manager is None
    tracks = genome_set_tracks(found)
    assert [a.get_title(loc="left") for a in figure.axes] == [t.title for t in tracks]
    (legend,) = figure.legends
    assert [t.get_text() for t in legend.get_texts()] == [
        t.get_text() for t in overview_figure(named("a")).legends[0].get_texts()
    ]
    assert "no candidates: b" in [t.get_text() for t in figure.axes[0].texts]


def test_a_genome_set_needs_at_least_one_genome() -> None:
    with pytest.raises(ValueError, match="at least one genome"):
        GenomeSetResults(())


def test_genomes_in_a_set_share_one_reference() -> None:
    other = named("b", families=ReferenceFamilies(FAMILIES.families[:3]))
    with pytest.raises(ValueError, match="same reference families"):
        GenomeSetResults((named("a"), other))
    assert GenomeSetResults((named("a"), named("b"))).families == FAMILIES


def test_items_pair_each_name_with_its_results() -> None:
    found = three_genomes()
    assert [(name, r.genome.name) for name, r in found.items()] == [
        ("a", "a"),
        ("b", "b"),
        ("c", "c"),
    ]


def files(root) -> dict[str, bytes]:
    return {
        str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()
    }


def test_write_gives_each_genome_its_folder_plus_a_combined_summary_and_plot(tmp_path) -> None:
    three_genomes().write(tmp_path / "out")
    written = files(tmp_path / "out")
    top = sorted(name for name in written if "/" not in name or name.startswith("plots/"))
    assert top == ["plots/regions.png", "summary.tsv"]
    assert sorted({name.split("/")[0] for name in written} - {"plots", "summary.tsv"}) == [
        "a",
        "b",
        "c",
    ]
    named("a").write(tmp_path / "alone")
    alone = files(tmp_path / "alone")
    assert {n[2:]: b for n, b in written.items() if n.startswith("a/")} == alone
    assert "b/summary.tsv" in written
    assert not [n for n in written if n.startswith("b/candidates/candidate_")]


def test_combined_summary_matches_the_combined_frame(tmp_path) -> None:
    found = three_genomes()
    found.write(tmp_path)
    lines = (tmp_path / "summary.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ["genome", *SUMMARY_COLUMNS]
    assert [line.split("\t")[:3] for line in lines[1:]] == [
        ["a", "1", "candidate_001"],
        ["a", "2", "candidate_002"],
        ["c", "1", "candidate_001"],
        ["c", "2", "candidate_002"],
    ]


def test_writing_a_genome_set_twice_gives_identical_bytes(tmp_path) -> None:
    three_genomes().write(tmp_path / "x")
    first = files(tmp_path / "x")
    three_genomes().write(tmp_path / "x")
    three_genomes().write(tmp_path / "y")
    assert files(tmp_path / "x") == first == files(tmp_path / "y")


def labelled(name: str, label: str, **overrides):
    found = named(name, **overrides)
    return replace(found, genome=replace(found.genome, display_name=label))


def labelled_set() -> GenomeSetResults:
    return GenomeSetResults(
        (
            labelled("sample_112", "E. coli CFT073"),
            labelled("sample_007", "E. coli O157:H7", candidates=()),
            named("sample_031"),
        )
    )


def test_labels_are_shown_while_names_stay_the_keys(tmp_path) -> None:
    found = labelled_set()
    assert found.names == ("sample_112", "sample_007", "sample_031")
    assert found.labels == ("E. coli CFT073", "E. coli O157:H7", "sample_031")
    assert found["E. coli CFT073"] is found["sample_112"]
    assert list(dict.fromkeys(found.to_dataframe()["genome"])) == ["E. coli CFT073", "sample_031"]
    tracks = genome_set_tracks(found)
    assert tracks[0].note == "no candidates: E. coli O157:H7"
    assert tracks[1].title.startswith("E. coli CFT073 · #1 candidate_001")
    found.write(tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir() if p.is_dir()) == [
        "plots",
        "sample_007",
        "sample_031",
        "sample_112",
    ]
    rows = (tmp_path / "summary.tsv").read_text().splitlines()[1:]
    assert rows[0].startswith("E. coli CFT073\t1\t")


def test_labels_must_be_unique_in_a_set() -> None:
    with pytest.raises(ValueError, match="duplicate genome label 'E. coli'"):
        GenomeSetResults((labelled("a", "E. coli"), labelled("b", "E. coli")))
    with pytest.raises(ValueError, match="duplicate genome label 'b'"):
        GenomeSetResults((labelled("a", "b"), named("b")))


def test_find_all_checks_labels_before_any_search(searches) -> None:
    twins = [replace(GENOME, name=n, display_name="E. coli") for n in ("a", "b")]
    with pytest.raises(ValueError, match="duplicate genome label 'E. coli'"):
        RegionFinder(REFERENCE).find_all(twins)
    assert searches == []
