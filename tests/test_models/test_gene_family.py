"""GeneFamily and ReferenceFamilies (spec 5.5)."""

import pytest

from homolocus import GeneFamily, ReferenceFamilies

A, B, C = "c|1-9|+", "c|11-19|+", "c|21-29|-"


def test_representative_is_the_first_member() -> None:
    family = GeneFamily("F001", (A, B))
    assert family.representative_id == A
    assert family.gene_ids == (A, B)


@pytest.mark.parametrize(
    ("family_id", "gene_ids", "message"),
    [("", (A,), "family ID"), ("F001", (), "at least one gene"), ("F001", (A, A), "duplicate")],
)
def test_illegal_families_are_rejected(
    family_id: str, gene_ids: tuple[str, ...], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        GeneFamily(family_id, gene_ids)


def test_family_of_finds_each_member() -> None:
    first, second = GeneFamily("F001", (A, C)), GeneFamily("F002", (B,))
    families = ReferenceFamilies((first, second))
    assert families.family_of(A) is first
    assert families.family_of(C) is first
    assert families.family_of(B) is second


def test_unknown_gene_raises_key_error() -> None:
    with pytest.raises(KeyError):
        ReferenceFamilies((GeneFamily("F001", (A,)),)).family_of(B)


def test_length_and_iteration_follow_the_families() -> None:
    families = ReferenceFamilies((GeneFamily("F001", (A, B)), GeneFamily("F002", (C,))))
    assert len(families) == 2
    assert [family.id for family in families] == ["F001", "F002"]


def test_no_families_is_allowed() -> None:
    assert len(ReferenceFamilies(())) == 0


def test_gene_in_two_families_is_rejected() -> None:
    with pytest.raises(ValueError, match="more than one family"):
        ReferenceFamilies((GeneFamily("F001", (A, B)), GeneFamily("F002", (B,))))


def test_duplicate_family_ids_are_rejected() -> None:
    with pytest.raises(ValueError, match="duplicate family ID"):
        ReferenceFamilies((GeneFamily("F001", (A,)), GeneFamily("F001", (B,))))


def test_a_family_may_carry_a_name_but_never_an_empty_one() -> None:
    assert GeneFamily("F001", (A,)).name is None
    assert GeneFamily("F001", (A,), name="sitA").name == "sitA"
    with pytest.raises(ValueError, match="name must not be empty"):
        GeneFamily("F001", (A,), name="")
