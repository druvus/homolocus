"""Target genome input: format detection and annotation modes (spec 2.2, 2.3)."""

import gzip
import warnings
from dataclasses import replace
from pathlib import Path

import pytest
from Bio.SeqFeature import FeatureLocation, SeqFeature

from homolocus import Genome
from homolocus.annotation.genbank import AnnotationWarning, read_genbank
from homolocus.annotation.genome_file import detect_format, read_genome
from homolocus.annotation.prodigal import predict_genes


def write_text(tmp_path: Path, text: str, name: str = "input.txt") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def gene_ids(genome: Genome) -> list[str]:
    return [g.id for g in genome.genes()]


def predicted_ids(contigs: dict[str, str]) -> list[str]:
    return [g.id for genes in predict_genes(contigs).values() for g in genes]


def contig_sequence(dna) -> str:
    return dna.random(500) + dna.orf(300) + dna.random(500)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("LOCUS       c1  12 bp    DNA\n", "genbank"),
        (">c1\nACGT\n", "fasta"),
        ("\n  \n>c1\nACGT\n", "fasta"),
    ],
)
def test_format_is_detected_from_the_first_non_empty_line(tmp_path, text, expected) -> None:
    assert detect_format(write_text(tmp_path, text, name="no_extension")) == expected


def test_unknown_format_is_rejected(tmp_path) -> None:
    with pytest.raises(ValueError, match="format"):
        detect_format(write_text(tmp_path, "##gff-version 3\n"))


def test_empty_file_is_rejected(tmp_path) -> None:
    with pytest.raises(ValueError, match="empty"):
        detect_format(write_text(tmp_path, "\n\n"))


def test_fasta_genes_are_predicted(tmp_path, dna) -> None:
    sequence = contig_sequence(dna)
    path = write_text(tmp_path, f">c1 plasmid one\n{sequence}\n>tiny\nACGTACGTAC\n", "x.fna")
    genome = read_genome(path)
    assert list(genome.contigs) == ["c1", "tiny"]
    assert gene_ids(genome) == predicted_ids({"c1": sequence, "tiny": "ACGTACGTAC"})
    assert gene_ids(genome)
    assert genome.contigs["tiny"].n_genes == 0
    assert genome.contigs["c1"].description == "plasmid one"
    assert genome.contigs["tiny"].description is None
    assert genome.source == str(path)
    assert genome.skipped_features == ()


def test_duplicate_fasta_ids_are_rejected(tmp_path) -> None:
    path = write_text(tmp_path, ">c1\nACGT\n>c1\nACGT\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_genome(path)


def test_genbank_with_cds_uses_supplied_annotation(write_genbank, make_record, dna) -> None:
    sequence = contig_sequence(dna)
    cds = SeqFeature(FeatureLocation(500, 1403, strand=1), type="CDS")
    path = write_genbank([make_record("c1", sequence, [cds])])
    genome = read_genome(path)
    assert genome == read_genbank(path)
    assert gene_ids(genome) == ["c1|501-1403|+"]


def test_genbank_without_cds_is_predicted(write_genbank, make_record, dna) -> None:
    sequence = contig_sequence(dna)
    gene_feature = SeqFeature(FeatureLocation(500, 1403, strand=1), type="gene")
    path = write_genbank([make_record("c1", sequence, [gene_feature])])
    genome = read_genome(path)
    assert gene_ids(genome) == predicted_ids({"c1": sequence})
    assert genome.contigs["c1"].description == "synthetic"


def test_pyrodigal_mode_ignores_supplied_cds(write_genbank, make_record, dna) -> None:
    sequence = contig_sequence(dna)
    features = [
        SeqFeature(FeatureLocation(0, 9, strand=1), type="CDS"),
        SeqFeature(FeatureLocation(9, 18, strand=1), type="CDS", qualifiers={"pseudo": [""]}),
    ]
    path = write_genbank([make_record("c1", sequence, features)])
    with warnings.catch_warnings():
        warnings.simplefilter("error", AnnotationWarning)
        genome = read_genome(path, annotation="pyrodigal")
    assert gene_ids(genome) == predicted_ids({"c1": sequence})
    assert genome.skipped_features == ()


def test_pyrodigal_mode_on_fasta(tmp_path, dna) -> None:
    sequence = contig_sequence(dna)
    path = write_text(tmp_path, f">c1\n{sequence}\n")
    assert gene_ids(read_genome(path, annotation="pyrodigal")) == predicted_ids({"c1": sequence})


def test_duplicate_genbank_ids_are_rejected_when_predicting(write_genbank, make_record) -> None:
    path = write_genbank([make_record("c1", "ACGT" * 10, []), make_record("c1", "ACGT" * 10, [])])
    with pytest.raises(ValueError, match="duplicate"):
        read_genome(path)


@pytest.mark.parametrize("annotation", ["genbank", "", "AUTO"])
def test_unknown_annotation_mode_is_rejected(tmp_path, annotation) -> None:
    path = write_text(tmp_path, ">c1\nACGT\n")
    with pytest.raises(ValueError, match="annotation"):
        read_genome(path, annotation=annotation)


def test_genome_is_named_after_its_file(tmp_path, write_genbank, make_record, dna) -> None:
    fasta = write_text(tmp_path, f">c1\n{contig_sequence(dna)}\n", "strain_A.fna")
    assert read_genome(fasta).name == "strain_A"
    cds = SeqFeature(FeatureLocation(500, 1403, strand=1), type="CDS")
    gbk = write_genbank([make_record("c1", contig_sequence(dna), [cds])])
    assert read_genome(gbk).name == read_genbank(gbk).name == "input"
    assert Genome.from_file(gbk, annotation="pyrodigal").name == "input"


def gzip_copy(path: Path, name: str) -> Path:
    """Write a gzip-compressed copy of a file next to it."""
    target = path.with_name(name)
    target.write_bytes(gzip.compress(path.read_bytes()))
    return target


def test_gzip_fasta_gives_the_same_genome_as_the_plain_file(tmp_path, dna) -> None:
    plain = write_text(tmp_path, f">c1\n{contig_sequence(dna)}\n", "strain_A.fna")
    compressed = gzip_copy(plain, "strain_A.fna.gz")
    genome = read_genome(compressed)
    assert genome.name == "strain_A"
    assert genome.source == str(compressed)
    assert replace(genome, source=str(plain)) == read_genome(plain)


def test_gzip_genbank_gives_the_same_genome_as_the_plain_file(
    write_genbank, make_record, dna
) -> None:
    cds = SeqFeature(FeatureLocation(500, 1403, strand=1), type="CDS")
    plain = write_genbank([make_record("c1", contig_sequence(dna), [cds])])
    compressed = gzip_copy(plain, "input.gbk.gz")
    for read in (read_genome, read_genbank):
        genome = read(compressed)
        assert genome.name == "input"
        assert replace(genome, source=str(plain)) == read(plain)


def test_compression_is_detected_from_the_content_not_the_extension(tmp_path) -> None:
    plain = write_text(tmp_path, ">c1\nACGT\n", "plain.txt")
    assert detect_format(gzip_copy(plain, "looks_plain.fna")) == "fasta"


@pytest.mark.parametrize(
    ("file_name", "name"),
    [("strain.fna.gz", "strain"), ("strain.FNA.GZ", "strain"), ("strain.gz", "strain")],
)
def test_a_trailing_gz_is_removed_from_the_genome_name(tmp_path, dna, file_name, name) -> None:
    plain = write_text(tmp_path, f">c1\n{contig_sequence(dna)}\n", "plain.fna")
    assert read_genome(gzip_copy(plain, file_name)).name == name


@pytest.mark.parametrize(
    "damage",
    [lambda data: data[:-20], lambda data: data[:2] + b"not a gzip header"],
    ids=["truncated", "bad-header"],
)
def test_a_damaged_gzip_file_is_rejected_naming_the_file(tmp_path, dna, damage) -> None:
    plain = write_text(tmp_path, f">c1\n{contig_sequence(dna)}\n", "plain.fna")
    compressed = gzip_copy(plain, "cut.fna.gz")
    compressed.write_bytes(damage(compressed.read_bytes()))
    with pytest.raises(ValueError, match=r"cut\.fna\.gz: damaged gzip file"):
        read_genome(compressed)
