"""Candidate evaluation: coverage, families, ranking (spec 7.7, 7.9, 8)."""

from dataclasses import replace

import pytest

from homolocus import Alignment, GeneFamily, GeneMatch
from homolocus.detection.evaluation import family_evidence

R1, R2, R3 = "ref|1-300|+", "ref|401-700|+", "ref|801-1100|-"
T1, T2, T3 = "t|1-300|+", "t|401-700|+", "t|801-1100|+"
TARGET_ORDER = {T1: 0, T2: 1, T3: 2}


def hit(
    reference: str,
    target: str,
    *spans: tuple[int, int],
    bit_score: float = 100.0,
    query_length: int = 100,
) -> GeneMatch:
    hsps = tuple(
        Alignment(
            reference,
            target,
            80.0,
            end - start,
            start,
            end,
            0,
            end - start,
            query_length,
            200,
            1e-30,
            bit_score,
        )
        for start, end in spans
    )
    return GeneMatch(reference, target, hsps)


def evidence_for(family: GeneFamily, *matches: GeneMatch, threshold: float = 0.5):
    return family_evidence(
        family, matches, target_order=TARGET_ORDER, min_reference_gene_coverage=threshold
    )


def test_split_fragments_combine_into_one_coverage() -> None:
    family = GeneFamily("F001", (R1,))
    result = evidence_for(family, hit(R1, T1, (0, 45)), hit(R1, T2, (45, 100)))
    assert result.coverage == 1.0
    assert result.represented
    assert result.copies == (T1, T2)


def test_overlapping_hsps_are_counted_once() -> None:
    family = GeneFamily("F001", (R1,))
    result = evidence_for(family, hit(R1, T1, (0, 40), (20, 60)), hit(R1, T2, (50, 70)))
    assert result.coverage == pytest.approx(0.7)


@pytest.mark.parametrize(("end", "represented"), [(50, True), (49, False)])
def test_represented_threshold_is_inclusive(end: int, represented: bool) -> None:
    result = evidence_for(GeneFamily("F001", (R1,)), hit(R1, T1, (0, end)))
    assert result.represented is represented
    assert result.best_copy_id == T1


def test_best_member_has_the_highest_coverage() -> None:
    family = GeneFamily("F001", (R1, R2))
    result = evidence_for(family, hit(R1, T1, (0, 30)), hit(R2, T2, (0, 90)))
    assert (result.best_member_id, result.coverage) == (R2, 0.9)
    assert result.best_copy_id == T2
    assert result.copies == (T1, T2)


def test_tied_members_resolve_by_genomic_order() -> None:
    family = GeneFamily("F001", (R1, R2))
    result = evidence_for(family, hit(R2, T1, (0, 80)), hit(R1, T2, (0, 80)))
    assert result.best_member_id == R1


def test_best_copy_has_the_highest_bit_score_to_the_best_member() -> None:
    family = GeneFamily("F001", (R1,))
    result = evidence_for(
        family,
        hit(R1, T1, (0, 60), bit_score=150.0),
        hit(R1, T3, (0, 90), bit_score=400.0),
        hit(R1, T2, (0, 90), bit_score=200.0),
    )
    assert result.best_copy_id == T3
    assert result.copies == (T1, T2, T3)
    assert result.copy_count == 3


def test_tied_copies_resolve_by_target_order() -> None:
    family = GeneFamily("F001", (R1,))
    result = evidence_for(family, hit(R1, T3, (0, 90)), hit(R1, T2, (0, 90)))
    assert result.best_copy_id == T2


def test_copies_of_other_members_count_but_do_not_pair() -> None:
    family = GeneFamily("F001", (R1, R2))
    result = evidence_for(
        family, hit(R1, T1, (0, 90), bit_score=50.0), hit(R2, T3, (0, 40), bit_score=900.0)
    )
    assert (result.best_member_id, result.best_copy_id) == (R1, T1)
    assert result.copies == (T1, T3)


def test_matches_to_other_families_are_ignored() -> None:
    result = evidence_for(GeneFamily("F001", (R1,)), hit(R3, T1, (0, 90)))
    assert (result.coverage, result.represented) == (0.0, False)
    assert (result.best_member_id, result.best_copy_id, result.copies) == (None, None, ())


def test_family_without_matches_is_missing() -> None:
    result = evidence_for(GeneFamily("F002", (R2,)))
    assert result.family_id == "F002"
    assert not result.represented


def test_evidence_carries_the_family_name_whether_present_or_missing() -> None:
    assert evidence_for(GeneFamily("F002", (R2,), name="sitB")).family_name == "sitB"
    present = evidence_for(GeneFamily("F001", (R1,), name="sitA"), hit(R1, T1, (0, 100)))
    assert (present.represented, present.family_name) == (True, "sitA")
    assert evidence_for(GeneFamily("F001", (R1,))).family_name is None


def test_matches_may_be_given_as_a_one_shot_iterator() -> None:
    family = GeneFamily("F001", (R1, R2))
    matches = iter([hit(R1, T1, (0, 30)), hit(R2, T2, (0, 90))])
    result = family_evidence(
        family, matches, target_order=TARGET_ORDER, min_reference_gene_coverage=0.5
    )
    assert (result.best_member_id, result.copies) == (R2, (T1, T2))


# Whole-candidate evaluation and detection

from homolocus import (  # noqa: E402
    Contig,
    DiamondSearch,
    Gene,
    Genome,
    ReferenceFamilies,
    ReferenceRegion,
    RegionDetection,
    ScoringWeights,
)
from homolocus.detection.evaluation import DetectionContext, detect  # noqa: E402

PROTEIN = 300
REF_GENES = tuple(
    Gene("ref", n * 1200, n * 1200 + 900, 1, "A" * 900, "M" * PROTEIN) for n in range(4)
)
REFERENCE = ReferenceRegion("ref", 0, 4500, "A" * 4500, REF_GENES)
FAMILIES = ReferenceFamilies(
    tuple(GeneFamily(f"F{n + 1:03d}", (g.id,)) for n, g in enumerate(REF_GENES))
)


def target(contig_id: str, start: int, strand: int = 1, length: int = 900) -> Gene:
    return Gene(contig_id, start, start + length, strand, "A" * length, "M" * 300)


def match(ref: Gene, tgt: Gene, *, query=(0, PROTEIN), identity=80.0, evalue=1e-40, bits=300.0):
    span = query[1] - query[0]
    hsp = Alignment(ref.id, tgt.id, identity, span, *query, 0, span, PROTEIN, PROTEIN, evalue, bits)
    return GeneMatch(ref.id, tgt.id, (hsp,))


def genome(*contigs: Contig) -> Genome:
    return Genome(contigs={c.id: c for c in contigs})


def contig(contig_id: str, *genes: Gene, length: int = 200_000) -> Contig:
    return Contig(contig_id, "A" * length, genes)


def context(**detection) -> DetectionContext:
    settings = {"min_matches": 2, "min_reference_coverage": 0.3, **detection}
    return DetectionContext(
        reference=REFERENCE,
        families=FAMILIES,
        search=DiamondSearch(),
        detection=RegionDetection(**settings),
        weights=ScoringWeights(),
    )


def colinear_copy(contig_id: str, offset: int, strand: int = 1):
    genes = [target(contig_id, offset + n * 1200, strand) for n in range(4)]
    if strand == -1:
        genes = [target(contig_id, offset + (3 - n) * 1200, -1) for n in range(4)]
    return genes, [match(r, t) for r, t in zip(REF_GENES, genes, strict=True)]


def test_colinear_copy_is_one_full_candidate() -> None:
    genes, matches = colinear_copy("t1", 10_000)
    (cluster,) = detect(context(), genome(contig("t1", *genes)), matches)
    assert (cluster.contig_id, cluster.start, cluster.end) == ("t1", 10_000, 14_500)
    assert (cluster.display_start, cluster.display_end) == (cluster.start, cluster.end)
    assert cluster.reference_coverage == 1.0
    assert cluster.components.order == cluster.components.orientation == 1.0
    assert cluster.components.spacing == 1.0
    assert all(cluster.role(g.id) == "matched" for g in genes)
    assert not cluster.touches_contig_edge


def test_reversed_copy_scores_like_the_forward_one() -> None:
    forward_genes, forward = colinear_copy("t1", 10_000)
    reversed_genes, backward = colinear_copy("t1", 10_000, strand=-1)
    (a,) = detect(context(), genome(contig("t1", *forward_genes)), forward)
    (b,) = detect(context(), genome(contig("t1", *reversed_genes)), backward)
    assert a.score == pytest.approx(b.score)


def test_min_matches_and_min_reference_coverage_must_both_hold() -> None:
    genes, matches = colinear_copy("t1", 10_000)
    world = genome(contig("t1", *genes[:2]))
    assert len(detect(context(min_matches=2, min_reference_coverage=0.5), world, matches[:2])) == 1
    assert detect(context(min_matches=3, min_reference_coverage=0.5), world, matches[:2]) == ()
    assert detect(context(min_matches=2, min_reference_coverage=0.51), world, matches[:2]) == ()


def test_min_matches_above_the_family_count_is_an_error() -> None:
    with pytest.raises(ValueError, match="min_matches 5 exceeds the 4 reference families"):
        context(min_matches=5)


def test_family_genes_must_belong_to_the_reference() -> None:
    stray = ReferenceFamilies((*FAMILIES.families, GeneFamily("F005", ("other|1-3|+",))))
    with pytest.raises(ValueError, match="not a gene of the reference"):
        DetectionContext(
            REFERENCE, stray, DiamondSearch(), RegionDetection(min_matches=2), ScoringWeights()
        )


def test_insertion_and_weak_fragment_are_additional_evidence() -> None:
    a, b = target("t1", 10_000), target("t1", 13_000)
    insert, fragment = target("t1", 11_000, length=600), target("t1", 12_000, length=600)
    matches = [
        match(REF_GENES[0], a, query=(0, 160)),
        match(REF_GENES[0], fragment, query=(160, 300)),
        match(REF_GENES[1], b),
    ]
    (cluster,) = detect(context(), genome(contig("t1", a, insert, fragment, b)), matches)
    assert (cluster.role(insert.id), cluster.role(fragment.id)) == ("additional", "additional")
    assert cluster.families[0].coverage == 1.0
    assert cluster.families[0].copies == (a.id, fragment.id)


def test_failing_matches_are_kept_but_do_not_count() -> None:
    a, b, weak = target("t1", 10_000), target("t1", 13_000), target("t1", 11_500)
    failing = match(REF_GENES[2], weak, evalue=1.0)
    matches = [match(REF_GENES[0], a), match(REF_GENES[1], b), failing]
    (cluster,) = detect(context(), genome(contig("t1", a, weak, b)), matches)
    assert failing in cluster.matches
    assert not cluster.families[2].represented


def test_gene_straddling_the_span_is_not_evidence() -> None:
    a, b = target("t1", 10_000), target("t1", 13_000)
    straddler = target("t1", 13_500, length=1000)
    matches = [
        match(REF_GENES[0], a),
        match(REF_GENES[1], b),
        match(REF_GENES[2], straddler, query=(0, 100)),
    ]
    (cluster,) = detect(context(), genome(contig("t1", a, b, straddler)), matches)
    assert straddler.id not in {g.id for g in cluster.evidence_genes}
    assert all(m.target_gene_id != straddler.id for m in cluster.matches)


def test_padding_changes_only_the_display_interval() -> None:
    genes, matches = colinear_copy("t1", 10_000)
    left, right = target("t1", 5_000), target("t1", 20_000)
    world = genome(contig("t1", left, *genes, right))
    (plain,) = detect(context(), world, matches)
    (padded,) = detect(context(padding_genes=1), world, matches)
    assert (padded.display_start, padded.display_end) == (5_000, 20_900)
    assert (padded.start, padded.end, padded.score) == (plain.start, plain.end, plain.score)
    assert padded.role(left.id) == "additional"


def test_candidates_on_separate_contigs_are_never_combined() -> None:
    first = [target("t1", 10_000), target("t1", 11_200)]
    second = [target("t2", 10_000), target("t2", 11_200)]
    matches = [match(r, t) for r, t in zip(REF_GENES, first + second, strict=True)]
    clusters = detect(context(), genome(contig("t1", *first), contig("t2", *second)), matches)
    assert [c.contig_id for c in clusters] == ["t1", "t2"]
    assert [len(c.event_ids) for c in clusters] == [2, 2]


def test_better_candidates_rank_first() -> None:
    weak_genes = [target("t1", 10_000), target("t1", 11_200)]
    strong_genes, strong = colinear_copy("t1", 100_000)
    weak = [match(r, t) for r, t in zip(REF_GENES, weak_genes, strict=False)]
    clusters = detect(context(), genome(contig("t1", *weak_genes, *strong_genes)), weak + strong)
    assert [c.start for c in clusters] == [100_000, 10_000]


def test_ties_break_by_contig_order_then_start() -> None:
    b_late, a_late, b_early = (
        colinear_copy("b", 50_000),
        colinear_copy("a", 50_000),
        colinear_copy("b", 10_000),
    )
    world = genome(contig("b", *b_late[0], *b_early[0]), contig("a", *a_late[0]))
    clusters = detect(context(), world, b_late[1] + a_late[1] + b_early[1])
    assert [(c.contig_id, c.start) for c in clusters] == [
        ("b", 10_000),
        ("b", 50_000),
        ("a", 50_000),
    ]


def test_ties_break_by_reference_coverage_before_position() -> None:
    three = [target("t1", 10_000 + n * 1200) for n in range(3)]
    two = [target("t1", 100_000 + n * 1200) for n in range(2)]
    matches = [match(r, t) for r, t in zip(REF_GENES, three, strict=False)]
    matches += [match(r, t) for r, t in zip(REF_GENES, two, strict=False)]
    flat = ScoringWeights(content=0, similarity=1, order=0, orientation=0, spacing=0)
    ctx = DetectionContext(
        REFERENCE, FAMILIES, DiamondSearch(), RegionDetection(min_matches=2), flat
    )
    clusters = detect(ctx, genome(contig("t1", *two, *three)), matches)
    assert [c.reference_coverage for c in clusters] == [0.75, 0.5]


def test_overlapping_windows_keep_only_the_best_unless_allowed() -> None:
    genes = [target("t1", 10_000 + n * 3000) for n in range(4)]
    matches = [match(r, t) for r, t in zip(REF_GENES, genes, strict=True)]
    world = genome(contig("t1", *genes))
    limited = {"max_region_bp": 6_000, "max_inter_match_bp_gap": 5_000}
    assert len(detect(context(**limited, max_candidate_overlap=1.0), world, matches)) == 3
    assert len(detect(context(**limited, max_candidate_overlap=0.0), world, matches)) == 2


def test_edge_flag_is_carried_into_the_cluster() -> None:
    genes, matches = colinear_copy("t1", 100)
    (cluster,) = detect(context(), genome(contig("t1", *genes)), matches)
    assert cluster.touches_contig_edge


def test_split_is_carried_into_the_cluster_without_changing_the_score() -> None:
    first, half_a = target("t1", 10_000), target("t1", 11_200, length=450)
    half_b = target("t1", 11_700, length=450)
    rest = [target("t1", 12_400), target("t1", 13_600)]
    matches = [
        match(REF_GENES[0], first),
        match(REF_GENES[1], half_a, query=(0, 150)),
        match(REF_GENES[1], half_b, query=(150, 300)),
        match(REF_GENES[2], rest[0]),
        match(REF_GENES[3], rest[1]),
    ]
    (cluster,) = detect(context(), genome(contig("t1", first, half_a, half_b, *rest)), matches)
    assert cluster.potential_split
    assert not cluster.potential_fusion
    (split,) = cluster.splits
    assert (split.reference_gene_id, split.fragment_ids) == (
        REF_GENES[1].id,
        (half_a.id, half_b.id),
    )
    assert split.adjacent
    assert cluster.families[1].coverage == 1.0
    assert replace(cluster, splits=()).score == cluster.score


def test_fusion_is_carried_into_the_cluster() -> None:
    first, second, fused = target("t1", 10_000), target("t1", 11_200), target("t1", 12_400)
    # Each reference protein aligns in full to one half of the fused target protein.
    fused_c = Alignment(
        REF_GENES[2].id, fused.id, 80.0, 300, 0, 300, 0, 150, 300, 300, 1e-40, 300.0
    )
    fused_d = Alignment(
        REF_GENES[3].id, fused.id, 80.0, 300, 0, 300, 150, 300, 300, 300, 1e-40, 300.0
    )
    matches = [
        match(REF_GENES[0], first),
        match(REF_GENES[1], second),
        GeneMatch(REF_GENES[2].id, fused.id, (fused_c,)),
        GeneMatch(REF_GENES[3].id, fused.id, (fused_d,)),
    ]
    (cluster,) = detect(context(), genome(contig("t1", first, second, fused)), matches)
    assert cluster.potential_fusion
    assert cluster.fusions[0].target_gene_id == fused.id
    assert cluster.fusions[0].family_ids == ("F003", "F004")
    assert not cluster.potential_split


def test_min_matches_may_equal_the_family_count() -> None:
    assert context(min_matches=4).detection.min_matches == 4


def test_genes_outside_the_display_interval_are_not_in_the_cluster() -> None:
    genes, matches = colinear_copy("t1", 10_000)
    left, right = target("t1", 5_000), target("t1", 20_000)
    (cluster,) = detect(context(), genome(contig("t1", left, *genes, right)), matches)
    for outside in (left, right):
        with pytest.raises(KeyError):
            cluster.role(outside.id)


def test_padding_never_adds_evidence() -> None:
    genes, matches = colinear_copy("t1", 10_000)
    genes, matches = genes[:3], matches[:3]
    neighbour = target("t1", 5_000)
    weak = match(REF_GENES[3], neighbour, query=(0, 100))
    world = genome(contig("t1", neighbour, *genes))
    (cluster,) = detect(context(padding_genes=1), world, [*matches, weak])
    assert cluster.role(neighbour.id) == "additional"
    assert weak not in cluster.matches
    assert cluster.families[3].copies == ()


def test_unrepresented_family_with_a_copy_does_not_score() -> None:
    genes = [target("t1", start) for start in (10_000, 12_000, 14_000)]
    middle = target("t1", 11_000, length=600)
    matches = [match(r, t) for r, t in zip(REF_GENES, genes, strict=False)]
    partial = match(REF_GENES[3], middle, query=(0, 140))
    world = genome(contig("t1", genes[0], middle, *genes[1:]))
    (with_partial,) = detect(context(), world, [*matches, partial])
    (without,) = detect(context(), world, matches)
    assert with_partial.families[3].best_copy_id == middle.id
    assert not with_partial.families[3].represented
    assert with_partial.components == without.components


@pytest.mark.diamond
def test_end_to_end_with_diamond_finds_a_rearranged_reversed_copy(tmp_path, proteins) -> None:
    from homolocus import FamilyGrouping  # noqa: PLC0415
    from homolocus.detection.families import group_families  # noqa: PLC0415
    from homolocus.models.gene_match import aggregate  # noqa: PLC0415
    from homolocus.search.diamond import search  # noqa: PLC0415
    from homolocus.search.protein_index import ProteinIndex  # noqa: PLC0415

    sequences = [proteins.random(300) for _ in range(4)]
    ref_genes = [
        Gene("ref", n * 1200, n * 1200 + 900, 1, "A" * 900, seq) for n, seq in enumerate(sequences)
    ]
    reference = ReferenceRegion("ref", 0, 4500, "A" * 4500, tuple(ref_genes))

    def copy(contig_id: str, start: int, strand: int, protein: str) -> Gene:
        return Gene(contig_id, start, start + 900, strand, "A" * 900, protein)

    # Reference A B C D. Target t1 (reversed): C', A', insertion, B'; D is missing.
    t1 = [
        copy("t1", 20_000, -1, proteins.mutate(sequences[2], 0.8)),
        copy("t1", 21_200, -1, proteins.mutate(sequences[0], 0.8)),
        copy("t1", 22_400, 1, proteins.random(300)),
        copy("t1", 23_600, -1, proteins.mutate(sequences[1], 0.8)),
    ]
    decoy = [copy("t2", 5_000, 1, proteins.mutate(sequences[3], 0.8))]
    world = Genome(
        contigs={"t1": Contig("t1", "A" * 40_000, t1), "t2": Contig("t2", "A" * 20_000, decoy)}
    )

    index = ProteinIndex.build(ref_genes, list(world.genes()))
    self_hits = search(index, DiamondSearch(), tmp_path / "self.tsv", self_search=True)
    families = group_families(ref_genes, self_hits, FamilyGrouping())
    matches = aggregate(search(index, DiamondSearch(), tmp_path / "main.tsv"))
    ctx = DetectionContext(
        reference, families, DiamondSearch(), RegionDetection(), ScoringWeights()
    )

    clusters = detect(ctx, world, matches)
    assert len(clusters) == 1
    (best,) = clusters
    assert best.contig_id == "t1"
    assert [f.family_id for f in best.represented_families] == ["F001", "F002", "F003"]
    assert [f.family_id for f in best.missing_families] == ["F004"]
    assert best.role(t1[2].id) == "additional"
    assert best.components.content == 0.75
    assert 0 < best.components.order < 1


@pytest.mark.diamond
def test_end_to_end_with_diamond_flags_a_split_protein(tmp_path, proteins) -> None:
    from homolocus import FamilyGrouping  # noqa: PLC0415
    from homolocus.detection.families import group_families  # noqa: PLC0415
    from homolocus.models.gene_match import aggregate  # noqa: PLC0415
    from homolocus.search.diamond import search  # noqa: PLC0415
    from homolocus.search.protein_index import ProteinIndex  # noqa: PLC0415

    sequences = [proteins.random(300) for _ in range(4)]
    ref_genes = [
        Gene("ref", n * 1200, n * 1200 + 900, 1, "A" * 900, seq) for n, seq in enumerate(sequences)
    ]
    reference = ReferenceRegion("ref", 0, 4500, "A" * 4500, tuple(ref_genes))
    split_b = proteins.mutate(sequences[1], 0.85)
    proteins_on_target = [
        proteins.mutate(sequences[0], 0.85),
        split_b[:150],
        split_b[150:],
        proteins.mutate(sequences[2], 0.85),
        proteins.mutate(sequences[3], 0.85),
    ]
    genes, start = [], 10_000
    for protein in proteins_on_target:
        length = 3 * len(protein)
        genes.append(Gene("t1", start, start + length, 1, "A" * length, protein))
        start += length + 300
    world = Genome(contigs={"t1": Contig("t1", "A" * 40_000, genes)})

    index = ProteinIndex.build(ref_genes, genes)
    self_hits = search(index, DiamondSearch(), tmp_path / "self.tsv", self_search=True)
    families = group_families(ref_genes, self_hits, FamilyGrouping())
    matches = aggregate(search(index, DiamondSearch(), tmp_path / "main.tsv"))
    ctx = DetectionContext(
        reference, families, DiamondSearch(), RegionDetection(), ScoringWeights()
    )

    (cluster,) = detect(ctx, world, matches)
    assert cluster.reference_coverage == 1.0
    assert cluster.potential_split
    (split,) = cluster.splits
    assert split.reference_gene_id == ref_genes[1].id
    assert split.fragment_ids == (genes[1].id, genes[2].id)
    assert split.adjacent
    assert split.same_strand
    assert not cluster.potential_fusion
