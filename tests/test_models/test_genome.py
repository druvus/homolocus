"""Genome model (spec 2.3, 5.3)."""

import pytest

from homolocus import Contig, Gene, Genome, SkippedFeature, SkipReason


def gene(contig_id: str, start: int, end: int) -> Gene:
    return Gene(
        contig_id=contig_id,
        start=start,
        end=end,
        strand=1,
        nucleotide_sequence="A" * (end - start),
        protein_sequence="M",
    )


def contig(contig_id: str, *genes: Gene) -> Contig:
    return Contig(id=contig_id, sequence="A" * 100, genes=genes)


def make_genome() -> Genome:
    return Genome(
        contigs={
            "c2": contig("c2", gene("c2", 50, 80), gene("c2", 0, 30)),
            "empty": contig("empty"),
            "c1": contig("c1", gene("c1", 10, 40)),
        },
        name="target",
        source="target.fna",
    )


def test_all_contigs_are_preserved_in_input_order() -> None:
    assert list(make_genome().contigs) == ["c2", "empty", "c1"]


def test_genes_iterates_contigs_in_order_and_genes_by_start() -> None:
    ids = [g.id for g in make_genome().genes()]
    assert ids == ["c2|1-30|+", "c2|51-80|+", "c1|11-40|+"]


def test_gene_lookup_by_canonical_id() -> None:
    found = make_genome().gene("c2|51-80|+")
    assert (found.contig_id, found.start, found.end) == ("c2", 50, 80)


def test_unknown_gene_id_raises_key_error() -> None:
    with pytest.raises(KeyError):
        make_genome().gene("c3|1-30|+")


def test_contigs_cannot_be_changed_through_the_input_mapping() -> None:
    contigs = {"c1": contig("c1")}
    genome = Genome(contigs=contigs)
    contigs["c2"] = contig("c2")
    assert list(genome.contigs) == ["c1"]


def test_metadata_defaults() -> None:
    genome = Genome(contigs={"c1": contig("c1")})
    assert genome.name is None
    assert genome.source is None


def test_key_must_match_contig_id() -> None:
    with pytest.raises(ValueError, match="key"):
        Genome(contigs={"c2": contig("c1")})


def test_genome_needs_a_contig() -> None:
    with pytest.raises(ValueError, match="at least one contig"):
        Genome(contigs={})


def test_genome_carries_skipped_features() -> None:
    skipped = (SkippedFeature("c1", 0, 30, None, SkipReason.PSEUDOGENE),)
    genome = Genome(contigs={"c1": contig("c1")}, skipped_features=skipped)
    assert genome.skipped_features == skipped
    assert Genome(contigs={"c1": contig("c1")}).skipped_features == ()


def test_from_file_reads_a_genome(tmp_path) -> None:
    path = tmp_path / "target.fna"
    path.write_text(">c1 first\nACGTACGTAC\n", encoding="utf-8")
    genome = Genome.from_file(path)
    assert list(genome.contigs) == ["c1"]
    assert genome.source == str(path)
    assert Genome.from_file(str(path), annotation="pyrodigal") == genome
