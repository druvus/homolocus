"""Contig model (spec 4, 5.2)."""

import pytest

from homolocus import Contig, Gene, Strand

SEQUENCE = "ACGT" * 25  # 100 bp


def gene(start: int, end: int, strand: Strand = 1, contig_id: str = "c1") -> Gene:
    return Gene(
        contig_id=contig_id,
        start=start,
        end=end,
        strand=strand,
        nucleotide_sequence="A" * (end - start),
        protein_sequence="M",
    )


def test_contig_without_genes_is_kept() -> None:
    contig = Contig(id="c1", sequence=SEQUENCE)
    assert contig.genes == ()
    assert contig.n_genes == 0
    assert contig.length == 100
    assert contig.description is None


def test_genes_are_sorted_by_start_irrespective_of_strand() -> None:
    genes = (gene(60, 90, -1), gene(0, 30, 1), gene(30, 45, -1))
    contig = Contig(id="c1", sequence=SEQUENCE, genes=genes)
    assert [g.start for g in contig.genes] == [0, 30, 60]
    assert contig.n_genes == 3


def test_ties_on_start_are_ordered_by_end_then_strand() -> None:
    genes = (gene(0, 30, 1), gene(0, 30, -1), gene(0, 15, 1))
    contig = Contig(id="c1", sequence=SEQUENCE, genes=genes)
    assert [g.id for g in contig.genes] == ["c1|1-15|+", "c1|1-30|-", "c1|1-30|+"]


def test_gene_may_end_at_contig_end() -> None:
    contig = Contig(id="c1", sequence=SEQUENCE, genes=(gene(70, 100),))
    assert contig.n_genes == 1


@pytest.mark.parametrize(
    ("contig_id", "sequence", "genes", "message"),
    [
        ("", SEQUENCE, (), "contig ID"),
        ("c1", "", (), "sequence"),
        ("c1", SEQUENCE, (gene(0, 30, contig_id="c2"),), "belongs to contig"),
        ("c1", SEQUENCE, (gene(80, 101),), "extends past"),
        ("c1", SEQUENCE, (gene(0, 30), gene(0, 30)), "duplicate"),
    ],
)
def test_illegal_contigs_are_rejected(
    contig_id: str, sequence: str, genes: tuple[Gene, ...], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        Contig(id=contig_id, sequence=sequence, genes=genes)


def test_duplicate_gene_error_names_the_gene() -> None:
    with pytest.raises(ValueError, match=r"contig c1 has duplicate gene IDs: c1\|1-30\|\+"):
        Contig(id="c1", sequence=SEQUENCE, genes=(gene(0, 30), gene(30, 45), gene(0, 30)))
