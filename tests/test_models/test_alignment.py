"""Alignment model: one DIAMOND HSP (spec 5.6)."""

from typing import Any

import pytest

from homolocus import Alignment

REF = "ref|1-300|+"
TARGET = "tgt|1-330|-"


def make_alignment(**overrides: Any) -> Alignment:
    fields: dict[str, Any] = {
        "query_id": REF,
        "subject_id": TARGET,
        "identity": 45.5,
        "alignment_length": 60,
        "query_start": 0,
        "query_end": 50,
        "subject_start": 10,
        "subject_end": 65,
        "query_length": 100,
        "subject_length": 110,
        "evalue": 1e-20,
        "bit_score": 120.0,
    }
    fields.update(overrides)
    return Alignment(**fields)


def test_alignment_keeps_its_evidence() -> None:
    hsp = make_alignment()
    assert (hsp.query_id, hsp.subject_id) == (REF, TARGET)
    assert hsp.query_interval == (0, 50)
    assert hsp.subject_interval == (10, 65)
    assert (hsp.identity, hsp.evalue, hsp.bit_score) == (45.5, 1e-20, 120.0)


def test_interval_may_cover_the_whole_protein() -> None:
    hsp = make_alignment(query_end=100, subject_start=0, subject_end=110, alignment_length=110)
    assert (hsp.query_interval, hsp.subject_interval) == ((0, 100), (0, 110))


@pytest.mark.parametrize("identity", [0.0, 100.0])
def test_identity_bounds_are_allowed(identity: float) -> None:
    assert make_alignment(identity=identity).identity == identity


def test_zero_evalue_and_bit_score_are_allowed() -> None:
    hsp = make_alignment(evalue=0.0, bit_score=0.0)
    assert (hsp.evalue, hsp.bit_score) == (0.0, 0.0)


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"query_id": ""}, "query ID"),
        ({"subject_id": ""}, "subject ID"),
        ({"identity": -0.1}, "identity"),
        ({"identity": 100.1}, "identity"),
        ({"query_length": 0}, "query length 0 must be positive"),
        ({"subject_length": 0}, "subject length 0 must be positive"),
        ({"query_start": 50, "query_end": 50}, "interval"),
        ({"query_start": -1}, "interval"),
        ({"subject_start": 70, "subject_end": 65}, "interval"),
        ({"query_end": 101}, "query interval"),
        ({"subject_end": 111}, "subject interval"),
        ({"alignment_length": 54}, "alignment length"),
        ({"alignment_length": 49, "subject_end": 59}, "alignment length"),
        ({"evalue": -1e-9}, "e-value"),
        ({"bit_score": -1.0}, "bit score"),
    ],
)
def test_illegal_alignments_are_rejected(overrides: dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        make_alignment(**overrides)


def test_span_is_measured_from_the_interval_start() -> None:
    hsp = make_alignment(query_start=20, query_end=60, alignment_length=55)
    assert hsp.query_interval == (20, 60)


def test_one_residue_proteins_are_allowed() -> None:
    hsp = make_alignment(
        query_end=1, query_length=1, subject_start=0, subject_end=1, subject_length=1
    )
    assert (hsp.query_length, hsp.subject_length) == (1, 1)
