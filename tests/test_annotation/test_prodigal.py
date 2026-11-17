"""Gene prediction with Pyrodigal (spec 2.2)."""

import pyrodigal
import pytest
from Bio.Seq import Seq

from homolocus import Gene
from homolocus.annotation.prodigal import predict_genes

TABLE_11_BINS = pyrodigal.MetagenomicBins(
    [b for b in pyrodigal.METAGENOMIC_BINS if b.training_info.translation_table == 11]
)


def signature(genes) -> list[tuple[int, int, int]]:
    return [(g.start, g.end, g.strand) for g in genes]


def oracle(finder: pyrodigal.GeneFinder, sequence: str) -> list[tuple[int, int, int]]:
    return [(g.begin - 1, g.end, g.strand) for g in finder.find_genes(sequence)]


def planted_contig(dna):
    forward, reverse = dna.orf(300), dna.orf(250)
    sequence = (
        dna.random(1000)
        + forward
        + dna.random(1000)
        + str(Seq(reverse).reverse_complement())
        + dna.random(1000)
    )
    forward_end = 1000 + len(forward)
    reverse_start = forward_end + 1000
    return sequence, forward_end, reverse_start


def test_planted_forward_and_reverse_orfs_are_predicted(dna) -> None:
    sequence, forward_end, reverse_start = planted_contig(dna)
    genes = predict_genes({"c1": sequence})["c1"]
    forward = [g for g in genes if g.strand == 1 and g.end == forward_end]
    reverse = [g for g in genes if g.strand == -1 and g.start == reverse_start]
    assert len(forward) == len(reverse) == 1
    assert forward[0].nucleotide_sequence == sequence[forward[0].start : forward[0].end]
    reverse_slice = sequence[reverse[0].start : reverse[0].end]
    assert reverse[0].nucleotide_sequence == str(Seq(reverse_slice).reverse_complement())


def test_proteins_have_no_stop_and_match_the_coding_sequence(dna) -> None:
    sequence, _, _ = planted_contig(dna)
    for gene in predict_genes({"c1": sequence})["c1"]:
        assert gene.protein_sequence is not None
        assert "*" not in gene.protein_sequence
        assert gene.protein_sequence[0] == "M" or gene.partial
        literal = str(Seq(gene.nucleotide_sequence).translate(table=11)).removesuffix("*")
        assert gene.protein_sequence[1:] == literal[1:]


def test_genes_are_plain_predictions_on_their_contig(dna) -> None:
    sequence, _, _ = planted_contig(dna)
    gene = predict_genes({"c1": sequence})["c1"][0]
    assert isinstance(gene, Gene)
    assert gene.contig_id == "c1"
    assert (gene.locus_tag, gene.gene_name, gene.product, gene.protein_id) == (None,) * 4


def test_orf_running_off_the_contig_end_is_partial(dna) -> None:
    sequence = dna.random(1000) + dna.orf(300)[:-33]
    genes = predict_genes({"c1": sequence})["c1"]
    edge = [g for g in genes if g.end == len(sequence)]
    assert len(edge) == 1
    assert edge[0].partial


def test_complete_genes_are_not_partial(dna) -> None:
    sequence, forward_end, _ = planted_contig(dna)
    genes = predict_genes({"c1": sequence})["c1"]
    assert not [g for g in genes if g.end == forward_end][0].partial


def test_short_contigs_use_table_11_meta_mode(dna) -> None:
    sequence, _, _ = planted_contig(dna)
    meta = pyrodigal.GeneFinder(meta=True, metagenomic_bins=TABLE_11_BINS)
    assert signature(predict_genes({"c1": sequence})["c1"]) == oracle(meta, sequence)


SMALL_TRAINING_SET = pytest.mark.filterwarnings("ignore:sequence should be at least")


@SMALL_TRAINING_SET
def test_long_contigs_use_single_mode_trained_on_all_contigs(dna) -> None:
    long_contig = "".join(planted_contig(dna)[0] for _ in range(5))  # > 20 kb
    short_contig = planted_contig(dna)[0]
    predicted = predict_genes({"long": long_contig, "short": short_contig})
    trained = pyrodigal.GeneFinder()
    trained.train(long_contig, short_contig)
    meta = pyrodigal.GeneFinder(meta=True, metagenomic_bins=TABLE_11_BINS)
    assert len(long_contig) >= pyrodigal.MIN_SINGLE_GENOME
    assert signature(predicted["long"]) == oracle(trained, long_contig)
    assert signature(predicted["short"]) == oracle(meta, short_contig)


@SMALL_TRAINING_SET
def test_contig_of_exactly_the_training_minimum_uses_single_mode(dna) -> None:
    contig = (planted_contig(dna)[0] * 5)[: pyrodigal.MIN_SINGLE_GENOME]
    trained = pyrodigal.GeneFinder()
    trained.train(contig)
    assert signature(predict_genes({"c1": contig})["c1"]) == oracle(trained, contig)


def test_contig_without_genes_is_kept_and_order_preserved(dna) -> None:
    sequence, _, _ = planted_contig(dna)
    predicted = predict_genes({"tiny": "ACGTACGTAC", "c1": sequence})
    assert list(predicted) == ["tiny", "c1"]
    assert predicted["tiny"] == ()
    assert predicted["c1"]
