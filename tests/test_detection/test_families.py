"""Paralogous reference families from the self-search (spec 7.6)."""

import pytest

from homolocus import Alignment, DiamondSearch, FamilyGrouping, Gene
from homolocus.detection.families import group_families
from homolocus.search.diamond import search
from homolocus.search.protein_index import ProteinIndex


def gene(number: int, protein: str | None = "M" * 100, strand: int = 1) -> Gene:
    start = number * 1000
    return Gene("ref", start, start + 300, strand, "A" * 300, protein)


A, B, C, D = gene(0), gene(1), gene(2, strand=-1), gene(3)


def hit(
    query: Gene,
    subject: Gene,
    *,
    identity: float = 95.0,
    query_span: tuple[int, int] = (0, 100),
    subject_span: tuple[int, int] = (0, 100),
) -> Alignment:
    length = max(query_span[1] - query_span[0], subject_span[1] - subject_span[0])
    return Alignment(
        query.id, subject.id, identity, length, *query_span, *subject_span, 100, 100, 1e-50, 300.0
    )


def members(families) -> list[tuple[str, tuple[str, ...]]]:
    return [(family.id, family.gene_ids) for family in families]


def test_paralogs_above_both_thresholds_form_one_family() -> None:
    families = group_families([A, B, C], [hit(A, B), hit(B, A)], FamilyGrouping())
    assert members(families) == [("F001", (A.id, B.id)), ("F002", (C.id,))]
    assert families.family_of(B.id).representative_id == A.id


def test_thresholds_are_inclusive() -> None:
    at_limit = hit(A, B, identity=90.0, query_span=(0, 80), subject_span=(20, 100))
    assert len(group_families([A, B], [at_limit], FamilyGrouping())) == 1


def test_identity_just_below_the_threshold_does_not_link() -> None:
    assert len(group_families([A, B], [hit(A, B, identity=89.9)], FamilyGrouping())) == 2


@pytest.mark.parametrize(
    ("query_span", "subject_span"),
    [((0, 79), (0, 100)), ((0, 100), (21, 100))],
)
def test_shared_domain_without_full_coverage_does_not_link(query_span, subject_span) -> None:
    domain = hit(A, B, identity=99.0, query_span=query_span, subject_span=subject_span)
    assert len(group_families([A, B], [domain], FamilyGrouping())) == 2


def test_one_passing_direction_is_enough() -> None:
    weak_reverse = hit(B, A, identity=50.0)
    assert len(group_families([A, B], [hit(A, B), weak_reverse], FamilyGrouping())) == 1


def test_coverage_of_split_hsps_is_combined() -> None:
    first = hit(A, B, query_span=(0, 45), subject_span=(0, 45))
    second = hit(A, B, query_span=(40, 90), subject_span=(40, 90))
    assert len(group_families([A, B], [first, second], FamilyGrouping())) == 1


def test_families_are_connected_components() -> None:
    families = group_families([A, B, C, D], [hit(A, B), hit(C, B)], FamilyGrouping())
    assert members(families) == [("F001", (A.id, B.id, C.id)), ("F002", (D.id,))]


def test_custom_thresholds_are_used() -> None:
    loose = FamilyGrouping(family_min_identity=40.0, family_min_coverage=0.5)
    weak = hit(A, B, identity=45.0, query_span=(0, 50), subject_span=(0, 50))
    assert len(group_families([A, B], [weak], loose)) == 1
    assert len(group_families([A, B], [weak], FamilyGrouping())) == 2


def test_disabled_grouping_gives_one_family_per_gene() -> None:
    disabled = FamilyGrouping(family_min_identity=None)
    families = group_families([A, B], [hit(A, B, identity=100.0)], disabled)
    assert members(families) == [("F001", (A.id,)), ("F002", (B.id,))]


def test_genes_without_protein_have_no_family() -> None:
    no_protein = gene(5, protein=None)
    families = group_families([A, no_protein], [], FamilyGrouping())
    assert members(families) == [("F001", (A.id,))]
    with pytest.raises(KeyError):
        families.family_of(no_protein.id)


def test_order_follows_the_genome_not_the_input() -> None:
    families = group_families([D, C, B, A], [hit(D, B)], FamilyGrouping())
    assert members(families) == [("F001", (A.id,)), ("F002", (B.id, D.id)), ("F003", (C.id,))]


def test_self_hits_and_unknown_genes_are_ignored() -> None:
    outside = gene(9)
    families = group_families(
        [A, B], [hit(A, A), hit(A, outside), hit(outside, B)], FamilyGrouping()
    )
    assert len(families) == 2


def test_no_alignments_gives_singletons_and_counts_each_family_once() -> None:
    assert len(group_families([A, B, C], [], FamilyGrouping())) == 3
    assert len(group_families([A, B, C], [hit(A, B), hit(B, C)], FamilyGrouping())) == 1


def test_family_ids_grow_beyond_three_digits() -> None:
    genes = [gene(n) for n in range(1000)]
    families = group_families(genes, [], FamilyGrouping())
    assert [f.id for f in families][-2:] == ["F999", "F1000"]


@pytest.mark.diamond
def test_self_search_groups_paralogs_but_not_shared_domains(tmp_path, proteins) -> None:
    family = proteins.random(300)
    domain = family[:100]
    sequences = [
        family,
        proteins.mutate(family, 0.97),
        domain + proteins.random(200),
        proteins.random(300),
    ]
    genes = [
        Gene("ref", n * 2000, n * 2000 + 900, 1, "A" * 900, sequence)
        for n, sequence in enumerate(sequences)
    ]
    index = ProteinIndex.build(genes, [])
    alignments = search(index, DiamondSearch(), tmp_path / "self.tsv", self_search=True)
    families = group_families(genes, alignments, FamilyGrouping())
    assert members(families) == [
        ("F001", (genes[0].id, genes[1].id)),
        ("F002", (genes[2].id,)),
        ("F003", (genes[3].id,)),
    ]


def named(number: int, gene_name: str | None = None, product: str | None = None) -> Gene:
    start = number * 1000
    return Gene(
        "ref", start, start + 300, 1, "A" * 300, "M" * 100, gene_name=gene_name, product=product
    )


def names(families) -> list[str | None]:
    return [family.name for family in families]


def test_family_names_come_from_gene_then_product_qualifiers() -> None:
    genes = [
        named(0, "sitA", "iron transporter"),
        named(1, None, "ABC transporter"),
        named(2),
    ]
    families = group_families(genes, [], FamilyGrouping())
    assert names(families) == ["sitA", "ABC transporter", None]


def test_paralog_names_are_unique_and_joined_in_genomic_order() -> None:
    genes = [named(0, "tnpA", "transposase"), named(1, "tnpB"), named(2, "tnpA")]
    links = [hit(genes[0], genes[1]), hit(genes[1], genes[2])]
    (family,) = group_families(genes, links, FamilyGrouping())
    assert family.name == "tnpA/tnpB"


def test_product_is_used_only_when_no_member_has_a_gene_name() -> None:
    genes = [named(0, None, "transposase"), named(1, "insB", "IS1 protein")]
    (family,) = group_families(genes, [hit(genes[0], genes[1])], FamilyGrouping())
    assert family.name == "insB"
    products = [named(0, None, "transposase"), named(1, None, "IS1 protein")]
    (family,) = group_families(products, [hit(products[0], products[1])], FamilyGrouping())
    assert family.name == "transposase/IS1 protein"


def test_families_are_named_with_grouping_disabled() -> None:
    genes = [named(0, "sitA"), named(1, "sitB")]
    families = group_families(genes, [hit(genes[0], genes[1])], FamilyGrouping(None))
    assert names(families) == ["sitA", "sitB"]
