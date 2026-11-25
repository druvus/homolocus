"""Printing a model never prints sequences or raw search output (spec 5)."""

from dataclasses import replace

from homolocus import Contig, Gene, Genome, GenomeSetResults
from sample_world import GENOME, REFERENCE, results

LARGE = 5_000_000


def large_genome() -> Genome:
    gene = Gene("c1", 3000, 3900, 1, "A" * 900, "M" * 300, qualifiers={"note": ["x" * 10_000]})
    contig = Contig("c1", "A" * LARGE, (gene,))
    return Genome({"c1": contig}, name="strain", display_name="E. coli K-12")


def test_genes_leave_sequences_and_qualifiers_out() -> None:
    (gene,) = large_genome().genes()
    text = repr(gene)
    assert text.startswith("Gene(contig_id='c1', start=3000, end=3900, strand=1, gene_name=None")
    for content in ("AAAA", "MMMM", "xxxx"):
        assert content not in text


def test_a_contig_prints_its_id_length_and_gene_count() -> None:
    assert repr(large_genome().contigs["c1"]) == "Contig(id='c1', length=5000000, genes=1)"


def test_a_genome_prints_its_name_label_contigs_and_genes() -> None:
    assert repr(large_genome()) == (
        "Genome(name='strain', label='E. coli K-12', contigs=1, genes=1)"
    )


def test_a_reference_prints_its_contig_one_based_span_and_genes() -> None:
    assert repr(REFERENCE) == "ReferenceRegion(contig='ref', span='1-4500', genes=4)"


def test_results_print_genome_families_candidates_and_best_score() -> None:
    assert repr(results()) == (
        "RegionResults(genome='target', families=4, candidates=2, best_score=0.827)"
    )
    assert repr(results(candidates=())) == (
        "RegionResults(genome='target', families=4, candidates=0, best_score=None)"
    )


def test_a_genome_set_prints_its_genome_and_candidate_counts() -> None:
    other = results(genome=replace(GENOME, name="other"), candidates=())
    found = GenomeSetResults((results(), other))
    assert repr(found) == "GenomeSetResults(genomes=2, candidates=2)"


def test_results_of_a_large_genome_print_in_one_short_line() -> None:
    text = repr(results(genome=large_genome(), candidates=(), raw_main="x" * LARGE))
    assert "\n" not in text
    assert len(text) < 200
