"""Biological scenarios end to end with real DIAMOND (spec 1, 17).

Each test builds a small reference region and a target from synthetic
proteins, runs RegionFinder with real DIAMOND and default thresholds, and
checks what a user would see in the results.
"""

import pytest

from scenario_world import Slot, World, genome, run

pytestmark = pytest.mark.diamond


def test_complete_colinear_copy_is_one_full_candidate() -> None:
    world = World(4)
    target = world.contig("t1", [Slot(world.copy(i)) for i in range(4)])
    found = run(world.reference(), genome(target))
    (candidate,) = found.candidates
    assert candidate.reference_coverage == 1.0
    assert (candidate.components.order, candidate.components.orientation) == (1.0, 1.0)
    assert candidate.event_ids == tuple(gene.id for gene in target.genes)


def test_reversed_copy_scores_like_the_forward_one_on_real_coordinates() -> None:
    world = World(4)
    copies = [world.copy(i) for i in range(4)]
    forward = world.contig("t1", [Slot(p) for p in copies])
    reverse = world.contig("t1", [Slot(p, strand=-1) for p in reversed(copies)])
    forward_hit = run(world.reference(), genome(forward)).candidates[0]
    reverse_hit = run(world.reference(), genome(reverse)).candidates[0]
    assert reverse_hit.components == forward_hit.components
    assert (reverse_hit.start, reverse_hit.end) == (reverse.genes[0].start, reverse.genes[-1].end)
    assert all(gene.strand == -1 for gene in reverse_hit.evidence_genes)


def test_unrelated_genome_gives_no_candidates_but_still_writes_output(tmp_path) -> None:
    world = World(4)
    target = world.contig("t1", [Slot(world.unrelated()) for _ in range(6)])
    found = run(world.reference(), genome(target))
    assert found.candidates == ()
    assert found.best_candidate is None
    found.write(tmp_path)
    assert (tmp_path / "summary.tsv").read_text().startswith("rank\tcandidate_id\t")
    assert (tmp_path / "plots" / "regions.png").is_file()


def test_missing_reference_genes_are_listed_and_lower_coverage() -> None:
    world = World(6)
    present = (0, 1, 2, 4)
    target = world.contig("t1", [Slot(world.copy(i)) for i in present])
    (candidate,) = run(world.reference(), genome(target)).candidates
    assert [f.family_id for f in candidate.missing_families] == ["F004", "F006"]
    assert candidate.reference_coverage == pytest.approx(4 / 6)
    assert candidate.components.order == 1.0


def test_inserted_unrelated_genes_are_additional_and_lower_spacing() -> None:
    world = World(4)
    slots = [
        Slot(world.copy(0)),
        Slot(world.unrelated()),
        Slot(world.unrelated()),
        Slot(world.copy(1)),
        Slot(world.copy(2)),
        Slot(world.unrelated(), strand=-1),
        Slot(world.copy(3)),
    ]
    target = world.contig("t1", slots)
    (candidate,) = run(world.reference(), genome(target)).candidates
    inserts = [target.genes[i].id for i in (1, 2, 5)]
    assert [candidate.role(gene_id) for gene_id in inserts] == ["additional"] * 3
    assert candidate.reference_coverage == 1.0
    assert candidate.components.spacing < 1.0


def test_rearranged_order_lowers_only_the_order_component() -> None:
    world = World(4)
    copies = [world.copy(i) for i in range(4)]
    colinear = world.contig("t1", [Slot(p) for p in copies])
    shuffled = world.contig("t1", [Slot(copies[i]) for i in (2, 0, 3, 1)])
    straight = run(world.reference(), genome(colinear)).candidates[0].components
    (candidate,) = run(world.reference(), genome(shuffled)).candidates
    assert candidate.reference_coverage == 1.0
    assert candidate.components.order < straight.order == 1.0
    assert candidate.components.content == straight.content
    assert candidate.components.similarity == pytest.approx(straight.similarity)
