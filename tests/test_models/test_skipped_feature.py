"""SkippedFeature model (spec 5.10)."""

import pytest

from homolocus import SkippedFeature, SkipReason


def test_skipped_feature_records_location_and_reason() -> None:
    feature = SkippedFeature("c1", 10, 40, "b0001", SkipReason.PSEUDOGENE)
    assert (feature.contig_id, feature.start, feature.end) == ("c1", 10, 40)
    assert feature.locus_tag == "b0001"
    assert feature.reason == "pseudogene"


def test_reason_values_are_stable_strings() -> None:
    assert [reason.value for reason in SkipReason] == [
        "pseudogene",
        "internal_stop",
        "compound_location",
    ]


@pytest.mark.parametrize(
    ("contig_id", "start", "end", "message"),
    [("", 0, 10, "contig"), ("c1", 10, 10, "interval"), ("c1", -1, 10, "interval")],
)
def test_illegal_skipped_features_are_rejected(
    contig_id: str, start: int, end: int, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        SkippedFeature(contig_id, start, end, None, SkipReason.COMPOUND_LOCATION)
