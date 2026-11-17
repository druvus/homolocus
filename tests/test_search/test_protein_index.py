"""DIAMOND tokens, FASTA input and the mapping file (spec 4.1, 6.5)."""

import pytest

from homolocus import Gene
from homolocus.search.protein_index import ProteinIndex


def gene(start: int, protein: str | None = "MKV", contig_id: str = "c1") -> Gene:
    return Gene(contig_id, start, start + 9, 1, "ATGAAAGTA", protein)


REFERENCE = [gene(0, "MKA"), gene(20, None), gene(40, "MKC")]
TARGET = [gene(0, "MTA", contig_id="t1"), gene(30, "MTB", contig_id="t1")]


def test_tokens_are_assigned_per_side_in_input_order() -> None:
    index = ProteinIndex.build(REFERENCE, TARGET)
    assert index.mapping_rows() == (
        ("r000001", "c1|1-9|+"),
        ("r000002", "c1|41-49|+"),
        ("t000001", "t1|1-9|+"),
        ("t000002", "t1|31-39|+"),
    )
    assert index.gene_id_for("r000002") == "c1|41-49|+"
    assert index.gene_id_for("t000001") == "t1|1-9|+"


def test_genes_without_protein_are_not_searchable() -> None:
    index = ProteinIndex.build(REFERENCE, TARGET)
    assert (index.n_reference, index.n_target) == (2, 2)
    assert "c1|21-29|+" not in index.reference_fasta()


def test_fasta_uses_tokens_as_headers() -> None:
    index = ProteinIndex.build(REFERENCE, TARGET)
    assert index.reference_fasta() == ">r000001\nMKA\n>r000002\nMKC\n"
    assert index.target_fasta() == ">t000001\nMTA\n>t000002\nMTB\n"


def test_the_same_gene_may_be_on_both_sides() -> None:
    shared = gene(0, "MKA")
    index = ProteinIndex.build([shared], [shared])
    assert index.gene_id_for("r000001") == index.gene_id_for("t000001") == shared.id


def test_unknown_token_raises_key_error() -> None:
    with pytest.raises(KeyError):
        ProteinIndex.build(REFERENCE, TARGET).gene_id_for("t000009")


def test_empty_sides() -> None:
    index = ProteinIndex.build([], [gene(0, None)])
    assert (index.n_reference, index.n_target) == (0, 0)
    assert index.reference_fasta() == index.target_fasta() == ""
    assert index.mapping_rows() == ()


def test_mapping_file_has_header_and_one_row_per_token(tmp_path) -> None:
    path = tmp_path / "protein_ids.tsv"
    ProteinIndex.build(REFERENCE, TARGET).write_mapping(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "token\tgene_id"
    assert lines[1:] == [
        "r000001\tc1|1-9|+",
        "r000002\tc1|41-49|+",
        "t000001\tt1|1-9|+",
        "t000002\tt1|31-39|+",
    ]
