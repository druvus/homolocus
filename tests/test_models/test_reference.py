"""ReferenceRegion model and GenBank input modes (spec 2.1, 5.4)."""

import warnings

import pytest
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

from homolocus import Gene, ReferenceRegion, SkippedFeature, SkipReason, Strand
from homolocus.annotation.genbank import AnnotationWarning

# 40 bp: CDS a at 1..12, CDS b (reverse) at 16..27, CDS c at 31..39 (1-based).
SEQUENCE = "ATGAAACCCTAA" + "GGG" + "TTAGGGTTTCAT" + "GGG" + "ATGAAATAA" + "G"


def cds(start: int, end: int, strand: int = 1, **qualifiers: str) -> SeqFeature:
    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type="CDS",
        qualifiers={key: [value] for key, value in qualifiers.items()},
    )


FEATURES = [cds(0, 12, locus_tag="a"), cds(15, 27, -1, locus_tag="b"), cds(30, 39, locus_tag="c")]


def gene(start: int, end: int, strand: Strand = 1, contig_id: str = "c1") -> Gene:
    return Gene(contig_id, start, end, strand, "A" * (end - start), "M")


def region(**overrides) -> ReferenceRegion:
    fields = {"contig_id": "c1", "start": 10, "end": 40, "sequence": "A" * 30}
    fields.update(overrides)
    return ReferenceRegion(**fields)


def tags(reference: ReferenceRegion) -> list[str | None]:
    return [g.locus_tag for g in reference.genes]


# Construction


def test_genes_are_sorted_and_keep_contig_coordinates() -> None:
    reference = region(genes=(gene(30, 40, -1), gene(10, 20)))
    assert [g.start for g in reference.genes] == [10, 30]
    assert reference.name is None
    assert reference.source_genome is None
    assert reference.skipped_features == ()


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"contig_id": ""}, "contig ID"),
        ({"start": 40, "end": 40, "sequence": ""}, "interval"),
        ({"sequence": "A" * 29}, "sequence length"),
        ({"genes": (gene(10, 20, contig_id="c2"),)}, "belongs to contig"),
        ({"genes": (gene(5, 20),)}, "extends past"),
        ({"genes": (gene(30, 41),)}, "extends past"),
        ({"genes": (gene(10, 20), gene(10, 20))}, "duplicate"),
        (
            {"skipped_features": (SkippedFeature("c1", 0, 20, None, SkipReason.PSEUDOGENE),)},
            "skipped feature",
        ),
    ],
)
def test_illegal_regions_are_rejected(overrides, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        region(**overrides)


# from_genbank


def test_coordinates_are_one_based_inclusive_and_genes_not_rebased(
    write_genbank, make_record
) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, FEATURES)])
    reference = ReferenceRegion.from_genbank(path, contig="c1", start=16, end=27)
    assert (reference.contig_id, reference.start, reference.end) == ("c1", 15, 27)
    assert reference.sequence == SEQUENCE[15:27]
    assert [g.id for g in reference.genes] == ["c1|16-27|-"]
    assert reference.source_genome == str(path)


def test_only_fully_contained_cds_are_included(write_genbank, make_record) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, FEATURES)])
    assert tags(ReferenceRegion.from_genbank(path, contig="c1", start=1, end=27)) == ["a", "b"]
    assert tags(ReferenceRegion.from_genbank(path, contig="c1", start=2, end=39)) == ["b", "c"]
    assert tags(ReferenceRegion.from_genbank(path, contig="c1", start=1, end=26)) == ["a"]


def test_whole_record_mode(write_genbank, make_record) -> None:
    path = write_genbank([make_record("pai", SEQUENCE, FEATURES)])
    reference = ReferenceRegion.from_genbank(path)
    assert (reference.contig_id, reference.start, reference.end) == ("pai", 0, 40)
    assert tags(reference) == ["a", "b", "c"]


def test_contig_is_selected_among_several_records(write_genbank, make_record) -> None:
    path = write_genbank(
        [make_record("c1", SEQUENCE, FEATURES[:1]), make_record("c2", SEQUENCE, FEATURES[1:])]
    )
    assert tags(ReferenceRegion.from_genbank(path, contig="c2")) == ["b", "c"]


def test_contig_is_required_when_the_file_has_several_records(write_genbank, make_record) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, []), make_record("c2", SEQUENCE, [])])
    with pytest.raises(ValueError, match="2 records"):
        ReferenceRegion.from_genbank(path)


def test_unknown_contig_is_rejected(write_genbank, make_record) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, [])])
    with pytest.raises(ValueError, match="'c9' not found"):
        ReferenceRegion.from_genbank(path, contig="c9")


def test_duplicate_selected_contig_is_rejected(write_genbank, make_record) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, []), make_record("c1", SEQUENCE, [])])
    with pytest.raises(ValueError, match="duplicate"):
        ReferenceRegion.from_genbank(path, contig="c1")


@pytest.mark.parametrize(("start", "end"), [(1, None), (None, 10)])
def test_start_and_end_go_together(write_genbank, make_record, start, end) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, [])])
    with pytest.raises(ValueError, match="together"):
        ReferenceRegion.from_genbank(path, start=start, end=end)


@pytest.mark.parametrize(("start", "end"), [(0, 10), (11, 10), (1, 41)])
def test_region_must_lie_within_the_record(write_genbank, make_record, start, end) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, [])])
    with pytest.raises(ValueError, match="interval|record"):
        ReferenceRegion.from_genbank(path, start=start, end=end)


def test_region_may_end_at_record_end(write_genbank, make_record) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, FEATURES)])
    assert ReferenceRegion.from_genbank(path, start=28, end=40).end == 40


def test_gene_without_protein_is_kept(write_genbank, make_record) -> None:
    path = write_genbank([make_record("c1", SEQUENCE, [cds(0, 12, transl_table="99")])])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AnnotationWarning)
        reference = ReferenceRegion.from_genbank(path)
    assert reference.genes[0].protein_sequence is None


def test_only_skipped_features_inside_the_region_are_kept(write_genbank, make_record) -> None:
    compound = SeqFeature(
        CompoundLocation([FeatureLocation(30, 33, 1), FeatureLocation(36, 39, 1)]), type="CDS"
    )
    features = [cds(0, 12, pseudo=""), compound]
    path = write_genbank([make_record("c1", SEQUENCE, features)])
    with pytest.warns(AnnotationWarning, match="compound_location: 1"):
        reference = ReferenceRegion.from_genbank(path, start=20, end=40)
    assert reference.skipped_features == (
        SkippedFeature("c1", 30, 39, None, SkipReason.COMPOUND_LOCATION),
    )


def test_sequence_length_error_names_both_lengths() -> None:
    with pytest.raises(ValueError, match="length 29 does not match region length 30"):
        region(sequence="A" * 29)


def test_skipped_feature_may_span_the_whole_region() -> None:
    skipped = (SkippedFeature("c1", 10, 40, None, SkipReason.PSEUDOGENE),)
    assert region(skipped_features=skipped).skipped_features == skipped


@pytest.mark.parametrize(
    "skipped",
    [
        SkippedFeature("c2", 10, 20, None, SkipReason.PSEUDOGENE),
        SkippedFeature("c1", 30, 41, None, SkipReason.PSEUDOGENE),
    ],
)
def test_skipped_feature_on_another_contig_or_past_the_end_is_rejected(skipped) -> None:
    with pytest.raises(ValueError, match="skipped feature"):
        region(skipped_features=(skipped,))
