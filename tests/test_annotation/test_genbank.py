"""GenBank parsing into a Genome (spec 2.2, 3.1, 3.2)."""

import warnings
from pathlib import Path

import pytest
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    FeatureLocation,
    SeqFeature,
)
from Bio.SeqRecord import SeqRecord

from homolocus import Genome, SkippedFeature, SkipReason
from homolocus.annotation.genbank import (
    AnnotationWarning,
    parse_record,
    read_genbank,
    record_contig,
)

# Forward CDS at 0..12 (ATG AAA CCC TAA -> MKP), spacer, reverse CDS at 15..27.
# Reverse strand of 15..27 is TTAGGGTTTCAT; its reverse complement is ATGAAACCCTAA.
FORWARD = "ATGAAACCCTAA"
REVERSE = "TTAGGGTTTCAT"
SEQUENCE = FORWARD + "GGG" + REVERSE + "GGG"


def cds(start: int, end: int, strand: int = 1, **qualifiers: str) -> SeqFeature:
    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type="CDS",
        qualifiers={key: [value] for key, value in qualifiers.items()},
    )


def read(write_genbank, make_record, *features: SeqFeature, sequence: str = SEQUENCE) -> Genome:
    path = write_genbank([make_record("c1", sequence, features)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AnnotationWarning)
        return read_genbank(path)


def only_gene(genome: Genome):
    genes = list(genome.genes())
    assert len(genes) == 1
    return genes[0]


def test_supplied_translation_is_used(write_genbank, make_record) -> None:
    genome = read(write_genbank, make_record, cds(0, 12, translation="MKQ"))
    assert only_gene(genome).protein_sequence == "MKQ"


def test_missing_translation_is_derived(write_genbank, make_record) -> None:
    genome = read(write_genbank, make_record, cds(0, 12))
    assert only_gene(genome).protein_sequence == "MKP"


def test_invalid_supplied_translation_falls_back_to_derived(write_genbank, make_record) -> None:
    genome = read(write_genbank, make_record, cds(0, 12, translation="MK1"))
    assert only_gene(genome).protein_sequence == "MKP"


def test_reverse_strand_cds(write_genbank, make_record) -> None:
    gene = only_gene(read(write_genbank, make_record, cds(15, 27, strand=-1)))
    assert gene.id == "c1|16-27|-"
    assert gene.nucleotide_sequence == "ATGAAACCCTAA"
    assert gene.protein_sequence == "MKP"


def test_transl_table_is_respected(write_genbank, make_record) -> None:
    sequence = "ATGTGATAA" + "GGG"
    genome = read(write_genbank, make_record, cds(0, 9, transl_table="4"), sequence=sequence)
    assert only_gene(genome).protein_sequence == "MW"


def test_codon_start_is_respected(write_genbank, make_record) -> None:
    sequence = "GAAACCCTAA" + "GGG"
    location = FeatureLocation(BeforePosition(0), 10, strand=1)
    feature = SeqFeature(location, type="CDS", qualifiers={"codon_start": ["2"]})
    gene = only_gene(read(write_genbank, make_record, feature, sequence=sequence))
    assert gene.protein_sequence == "KP"


def test_unknown_translation_table_gives_untranslatable_gene(write_genbank, make_record) -> None:
    gene = only_gene(read(write_genbank, make_record, cds(0, 12, transl_table="99")))
    assert gene.protein_sequence is None


@pytest.mark.parametrize(
    ("start", "end", "partial"),
    [
        (BeforePosition(0), 12, True),
        (0, AfterPosition(12), True),
        (0, 12, False),
    ],
)
def test_partial_cds_is_flagged(write_genbank, make_record, start, end, partial) -> None:
    feature = SeqFeature(FeatureLocation(start, end, strand=1), type="CDS")
    assert only_gene(read(write_genbank, make_record, feature)).partial is partial


def test_metadata_and_qualifiers(write_genbank, make_record) -> None:
    feature = cds(0, 12, gene="abcA", product="ABC protein", locus_tag="b0001", protein_id="P1")
    gene = only_gene(read(write_genbank, make_record, feature))
    assert (gene.gene_name, gene.product, gene.locus_tag, gene.protein_id) == (
        "abcA",
        "ABC protein",
        "b0001",
        "P1",
    )
    assert gene.qualifiers["locus_tag"] == ["b0001"]


@pytest.mark.parametrize("qualifier", ["pseudo", "pseudogene"])
def test_pseudogenes_are_excluded_and_recorded(write_genbank, make_record, qualifier) -> None:
    feature = cds(0, 12, locus_tag="b0001", **{qualifier: "unknown"})
    genome = read(write_genbank, make_record, feature, cds(15, 27, strand=-1))
    assert [g.id for g in genome.genes()] == ["c1|16-27|-"]
    assert genome.skipped_features == (SkippedFeature("c1", 0, 12, "b0001", SkipReason.PSEUDOGENE),)


@pytest.mark.parametrize("translation", [None, "MK*P"])
def test_internal_stop_is_excluded_and_recorded(write_genbank, make_record, translation) -> None:
    sequence = "ATGTAACCCTAA" + "GGG"
    qualifiers = {} if translation is None else {"translation": translation}
    genome = read(write_genbank, make_record, cds(0, 12, **qualifiers), sequence=sequence)
    assert list(genome.genes()) == []
    assert genome.skipped_features[0].reason is SkipReason.INTERNAL_STOP


def test_compound_cds_is_skipped_and_recorded(write_genbank, make_record) -> None:
    location = CompoundLocation([FeatureLocation(0, 6, 1), FeatureLocation(9, 12, 1)])
    genome = read(write_genbank, make_record, SeqFeature(location, type="CDS"))
    assert list(genome.genes()) == []
    assert genome.skipped_features == (
        SkippedFeature("c1", 0, 12, None, SkipReason.COMPOUND_LOCATION),
    )


def test_non_cds_features_are_ignored(write_genbank, make_record) -> None:
    other = SeqFeature(FeatureLocation(0, 12, strand=1), type="gene")
    genome = read(write_genbank, make_record, other)
    assert list(genome.genes()) == []
    assert genome.skipped_features == ()


def test_all_records_become_contigs(write_genbank, make_record) -> None:
    path = write_genbank(
        [
            make_record("c1", SEQUENCE, [cds(0, 12)]),
            make_record("empty", "ACGTACGT", []),
        ]
    )
    genome = read_genbank(path)
    assert list(genome.contigs) == ["c1", "empty"]
    assert genome.contigs["empty"].n_genes == 0
    assert genome.contigs["c1"].description == "synthetic"
    assert genome.source == str(path)


def test_duplicate_record_ids_are_rejected(write_genbank, make_record) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, []), make_record("c1", SEQUENCE, [])])
    with pytest.raises(ValueError, match="duplicate"):
        read_genbank(path)


def test_one_summary_warning_per_contig(write_genbank, make_record) -> None:
    features = [
        cds(0, 12, pseudo=""),
        cds(15, 27, strand=-1, transl_table="99"),
        SeqFeature(FeatureLocation(0, 11, strand=1), type="CDS"),
    ]
    path = write_genbank([make_record("c1", SEQUENCE, features)])
    with pytest.warns(AnnotationWarning) as record:
        read_genbank(path)
    assert len(record) == 1
    message = str(record[0].message)
    assert "c1" in message
    assert "pseudogene: 1" in message
    assert "1 untranslatable" in message
    assert "1 trimmed" in message
    assert message.startswith("contig c1: pseudogene: 1;")


def test_clean_contig_gives_no_warning(write_genbank, make_record) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, [cds(0, 12)])])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        read_genbank(path)


def test_path_may_be_a_string(write_genbank, make_record) -> None:
    path: Path = write_genbank([make_record("c1", SEQUENCE, [cds(0, 12)])])
    assert read_genbank(str(path)).source == str(path)


def test_cds_without_location_is_counted_in_the_warning(make_record) -> None:
    record = make_record("c1", SEQUENCE, [SeqFeature(None, type="CDS"), cds(0, 12)])
    with pytest.warns(AnnotationWarning, match="1 with an unparseable location"):
        parsed = parse_record(record)
    assert [g.id for g in parsed.genes] == ["c1|1-12|+"]


def test_record_without_id_is_rejected(make_record) -> None:
    record = make_record("c1", SEQUENCE, [])
    record.id = ""
    with pytest.raises(ValueError, match="without an ID"):
        parse_record(record)


ALT_START = "GTGAAATAA"  # translates to MK when 5'-complete, VK when 5'-partial
ALT_START_REVERSE = "TTATTTCAC"  # reverse complement of ALT_START


@pytest.mark.parametrize(
    ("start", "end", "protein"),
    [(0, 9, "MK"), (BeforePosition(0), 9, "VK"), (0, AfterPosition(9), "MK")],
)
def test_forward_five_prime_end_decides_alternative_start(
    write_genbank, make_record, start, end, protein
) -> None:
    feature = SeqFeature(FeatureLocation(start, end, strand=1), type="CDS")
    genome = read(write_genbank, make_record, feature, sequence=ALT_START + "GGG")
    assert only_gene(genome).protein_sequence == protein


@pytest.mark.parametrize(
    ("start", "end", "protein"),
    [(0, 9, "MK"), (0, AfterPosition(9), "VK"), (BeforePosition(0), 9, "MK")],
)
def test_reverse_five_prime_end_decides_alternative_start(
    write_genbank, make_record, start, end, protein
) -> None:
    feature = SeqFeature(FeatureLocation(start, end, strand=-1), type="CDS")
    genome = read(write_genbank, make_record, feature, sequence=ALT_START_REVERSE + "GGG")
    assert only_gene(genome).protein_sequence == protein


@pytest.mark.parametrize(
    ("qualifier", "value"),
    [("codon_start", "x"), ("codon_start", "4"), ("transl_table", "x"), ("transl_table", "99")],
)
def test_invalid_numeric_qualifier_gives_untranslatable_gene(
    write_genbank, make_record, qualifier, value
) -> None:
    genome = read(write_genbank, make_record, cds(0, 12, **{qualifier: value}))
    assert only_gene(genome).protein_sequence is None


def test_region_keeps_only_fully_contained_cds(make_record) -> None:
    features = [cds(0, 12), cds(15, 27, strand=-1), cds(3, 9)]
    record = make_record("c1", SEQUENCE, features)
    assert {g.id for g in parse_record(record, region=(0, 12)).genes} == {"c1|1-12|+", "c1|4-9|+"}
    assert {g.id for g in parse_record(record, region=(1, 27)).genes} == {"c1|4-9|+", "c1|16-27|-"}
    assert parse_record(record, region=(0, 11)).genes == [parse_record(record).genes[2]]


def test_region_filters_skipped_features_too(make_record) -> None:
    record = make_record("c1", SEQUENCE, [cds(0, 12, pseudo=""), cds(15, 27, strand=-1)])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        parsed = parse_record(record, region=(15, 27))
    assert parsed.skipped == []


def test_record_without_sequence_is_rejected() -> None:
    record = SeqRecord(None, id="c1", name="c1", description="no sequence")
    with pytest.raises(ValueError, match="no sequence"):
        record_contig(record, [])
