"""Biological scenarios end to end with real DIAMOND (spec 1, 17).

Each test builds a small reference region and a target from synthetic
proteins, runs RegionFinder with real DIAMOND and default thresholds, and
checks what a user would see in the results.
"""

import csv

import pytest

from homolocus import RegionDetection
from scenario_world import Slot, World, blocks, genome, run

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


def test_duplicated_homologue_keeps_both_copies_and_counts_the_family_once() -> None:
    world = World(4)
    slots = [Slot(world.copy(0)), Slot(world.copy(1, 0.95)), Slot(world.copy(1, 0.7))]
    target = world.contig("t1", [*slots, Slot(world.copy(2)), Slot(world.copy(3))])
    (candidate,) = run(world.reference(), genome(target)).candidates
    strong, weak = target.genes[1].id, target.genes[2].id
    assert {strong, weak} <= set(candidate.event_ids)
    (family,) = [f for f in candidate.families if f.family_id == "F002"]
    assert set(family.copies) >= {strong, weak}
    assert family.best_copy_id == strong
    assert candidate.reference_coverage == 1.0


def test_close_paralogs_form_one_family_counted_once() -> None:
    world = World(4)
    world.make_paralog(3, of=1)
    found = run(
        world.reference(), genome(world.contig("t1", [Slot(world.copy(i)) for i in range(3)]))
    )
    reference = world.reference()
    assert [f.gene_ids for f in found.families] == [
        (reference.genes[0].id,),
        (reference.genes[1].id, reference.genes[3].id),
        (reference.genes[2].id,),
    ]
    (candidate,) = found.candidates
    assert candidate.reference_coverage == 1.0
    assert candidate.missing_families == ()


def test_two_copies_on_one_contig_are_ranked_by_strength() -> None:
    world = World(4)
    strong = [Slot(world.copy(i, 0.95)) for i in range(4)]
    weak = [Slot(world.copy(i, 0.7)) for i in range(4)]
    target = blocks("t1", (200_000, weak), (10_000, strong), length=300_000)
    first, second = run(world.reference(), genome(target)).candidates
    assert first.start == 10_000
    assert second.start == 200_000
    assert first.score > second.score


def gene_flags(output, gene_id: str) -> set[str]:
    """Flags of one gene in the written candidate_001.tsv."""
    with (output / "candidates" / "candidate_001.tsv").open() as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["target_gene_id"] == gene_id
        ]
    return {flag for row in rows for flag in row["flags"].split(";") if flag}


def test_split_protein_is_flagged_and_still_represents_its_gene(tmp_path) -> None:
    world = World(4)
    halves = [Slot(world.part(1, 0, 150)), Slot(world.part(1, 150, 300))]
    target = world.contig(
        "t1", [Slot(world.copy(0)), *halves, Slot(world.copy(2)), Slot(world.copy(3))]
    )
    found = run(world.reference(), genome(target))
    (candidate,) = found.candidates
    assert candidate.potential_split
    assert candidate.reference_coverage == 1.0
    found.write(tmp_path)
    assert "split" in found.to_dataframe()["flags"][0].split(";")
    for fragment in target.genes[1:3]:
        assert "split_fragment" in gene_flags(tmp_path, fragment.id)


def test_fused_protein_is_flagged_on_the_candidate_and_the_gene(tmp_path) -> None:
    world = World(4)
    fused = Slot(world.copy(1) + world.copy(2))
    target = world.contig("t1", [Slot(world.copy(0)), fused, Slot(world.copy(3))])
    found = run(world.reference(), genome(target))
    (candidate,) = found.candidates
    assert candidate.potential_fusion
    assert candidate.reference_coverage == 1.0
    found.write(tmp_path)
    assert "fusion" in found.to_dataframe()["flags"][0].split(";")
    assert "fused" in gene_flags(tmp_path, target.genes[1].id)


def test_region_at_a_contig_end_is_flagged_edge() -> None:
    world = World(4)
    slots = [Slot(world.copy(i)) for i in range(4)]
    at_start = world.contig("t1", slots, offset=100)
    inside = world.contig("t1", slots)
    assert run(world.reference(), genome(at_start)).candidates[0].touches_contig_edge
    assert not run(world.reference(), genome(inside)).candidates[0].touches_contig_edge


def test_halves_on_two_contigs_are_never_combined() -> None:
    world = World(4)
    left = world.contig("t1", [Slot(world.copy(0)), Slot(world.copy(1))])
    right = world.contig("t2", [Slot(world.copy(2)), Slot(world.copy(3))])
    reference = world.reference()
    assert run(reference, genome(left, right)).candidates == ()
    # min_matches=2 lets each half pass on its own contig.
    halves = run(reference, genome(left, right), detection=RegionDetection(min_matches=2))
    assert sorted(c.contig_id for c in halves.candidates) == ["t1", "t2"]
    for candidate in halves.candidates:
        assert {gene.contig_id for gene in candidate.genes} == {candidate.contig_id}
    assert [c.reference_coverage for c in halves.candidates] == [0.5, 0.5]
