"""Gene model (spec 4.1, 5.1)."""

from dataclasses import FrozenInstanceError
from typing import Any

import pytest

from homolocus import Gene


def make_gene(**overrides: Any) -> Gene:
    fields: dict[str, Any] = {
        "contig_id": "NC_000913.3",
        "start": 1_900_000,
        "end": 1_900_009,
        "strand": 1,
        "nucleotide_sequence": "ATGAAATAA",
        "protein_sequence": "MK",
    }
    fields.update(overrides)
    return Gene(**fields)


def test_forward_gene_id_is_canonical_one_based() -> None:
    assert make_gene().id == "NC_000913.3|1900001-1900009|+"


def test_reverse_gene_id_uses_minus() -> None:
    assert make_gene(strand=-1).id == "NC_000913.3|1900001-1900009|-"


def test_id_ignores_locus_tag() -> None:
    assert make_gene(locus_tag="b1234").id == make_gene(locus_tag="other").id


def test_length_is_end_minus_start() -> None:
    assert make_gene().length == 9


def test_metadata_defaults() -> None:
    gene = make_gene()
    assert gene.gene_name is None
    assert gene.product is None
    assert gene.locus_tag is None
    assert gene.protein_id is None
    assert gene.partial is False
    assert gene.qualifiers == {}


def test_gene_without_protein_is_allowed() -> None:
    assert make_gene(protein_sequence=None).protein_sequence is None


def test_gene_is_immutable() -> None:
    gene = make_gene()
    with pytest.raises(FrozenInstanceError):
        gene.start = 0  # type: ignore[misc]


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"start": -1, "end": 8}, "interval"),
        ({"start": 9, "end": 9, "nucleotide_sequence": ""}, "interval"),
        ({"start": 10, "end": 1}, "interval"),
        ({"strand": 0}, "strand"),
        ({"strand": 2}, "strand"),
        ({"nucleotide_sequence": "ATGAAA"}, "nucleotide sequence"),
        ({"nucleotide_sequence": "ATGAAATAAT"}, "nucleotide sequence"),
        ({"protein_sequence": ""}, "protein sequence"),
        ({"contig_id": ""}, "contig"),
    ],
)
def test_illegal_genes_are_rejected(overrides: dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        make_gene(**overrides)
