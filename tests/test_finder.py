"""RegionFinder: the whole pipeline from reference and genome to results (spec 12)."""

import warnings

import pytest

from homolocus import (
    Contig,
    DiamondSearch,
    FamilyGrouping,
    Gene,
    Genome,
    ReferenceRegion,
    RegionDetection,
    RegionFinder,
    RegionResults,
)
from homolocus.search.diamond import DiamondNotFoundError


def build_world(proteins):
    """Reference A B C D (paralogs C and D) and a target with a reversed copy of A B C."""
    family = proteins.random(300)
    sequences = [proteins.random(300), proteins.random(300), family, proteins.mutate(family, 0.97)]
    ref_genes = tuple(
        Gene("ref", n * 1200, n * 1200 + 900, 1, "A" * 900, seq) for n, seq in enumerate(sequences)
    )
    reference = ReferenceRegion("ref", 0, 4500, "A" * 4500, ref_genes, source_genome="ref.gbk")
    target_genes = (
        Gene("t1", 20_000, 20_900, -1, "A" * 900, proteins.mutate(sequences[2], 0.8)),
        Gene("t1", 21_200, 22_100, -1, "A" * 900, proteins.mutate(sequences[1], 0.8)),
        Gene("t1", 22_400, 23_300, -1, "A" * 900, proteins.mutate(sequences[0], 0.8)),
        Gene("t1", 30_000, 30_900, 1, "A" * 900, proteins.random(300)),
    )
    genome = Genome(
        contigs={"t1": Contig("t1", "A" * 40_000, target_genes)},
        name="target",
        source="target.fna",
    )
    return reference, genome


def test_missing_diamond_gives_a_clear_error(proteins, monkeypatch) -> None:
    reference, genome = build_world(proteins)
    monkeypatch.setenv("PATH", "")
    with pytest.raises(DiamondNotFoundError, match="not found"):
        RegionFinder(reference).find(genome)


def test_defaults_are_the_configuration_defaults(proteins) -> None:
    reference, _ = build_world(proteins)
    finder = RegionFinder(reference)
    assert finder.search == DiamondSearch()
    assert finder.family_grouping == FamilyGrouping()
    assert finder.detection == RegionDetection()


@pytest.mark.diamond
def test_find_runs_the_whole_pipeline(proteins) -> None:
    reference, genome = build_world(proteins)
    results = RegionFinder(reference).find(genome)
    assert isinstance(results, RegionResults)
    assert [f.gene_ids for f in results.families] == [
        (reference.genes[0].id,),
        (reference.genes[1].id,),
        (reference.genes[2].id, reference.genes[3].id),
    ]
    best = results.best_candidate
    assert best is not None
    assert best.contig_id == "t1"
    assert best.reference_coverage == 1.0
    assert best.components.order == 1.0
    assert results.diamond_version.startswith("2.")
    assert results.raw_main.startswith("r0000")
    assert results.raw_self
    assert results.matches
    assert results.warnings == ()


@pytest.mark.diamond
def test_written_output_is_identical_across_runs(tmp_path, proteins) -> None:
    reference, genome = build_world(proteins)
    finder = RegionFinder(reference, search=DiamondSearch(threads=1))
    finder.find(genome).write(tmp_path / "first")
    finder.find(genome).write(tmp_path / "second")

    def tree(root):
        return {p.relative_to(root): p.read_bytes() for p in root.rglob("*") if p.is_file()}

    first = tree(tmp_path / "first")
    assert first == tree(tmp_path / "second")
    assert any(str(path).startswith("candidates/") for path in first)


@pytest.mark.diamond
def test_disabled_grouping_skips_the_self_search(proteins) -> None:
    reference, genome = build_world(proteins)
    finder = RegionFinder(reference, family_grouping=FamilyGrouping(family_min_identity=None))
    results = finder.find(genome)
    assert results.raw_self == ""
    assert len(results.families) == 4


@pytest.mark.diamond
def test_min_matches_above_the_family_count_is_an_error(proteins) -> None:
    reference, genome = build_world(proteins)
    finder = RegionFinder(reference, detection=RegionDetection(min_matches=4))
    with pytest.raises(ValueError, match="min_matches 4 exceeds the 3 reference families"):
        finder.find(genome)


@pytest.mark.diamond
def test_warnings_during_a_run_are_kept_and_still_shown(proteins, monkeypatch) -> None:
    import homolocus.finder as finder_module  # noqa: PLC0415

    original = finder_module.detect

    def noisy_detect(*args, **kwargs):
        warnings.warn("something to note", UserWarning, stacklevel=2)
        return original(*args, **kwargs)

    monkeypatch.setattr(finder_module, "detect", noisy_detect)
    reference, genome = build_world(proteins)
    with pytest.warns(UserWarning, match="something to note"):
        results = RegionFinder(reference).find(genome)
    assert results.warnings == ("something to note",)
