"""Conversion between internal and user-facing coordinates (spec 4)."""

import pytest
from hypothesis import given
from hypothesis import strategies as st

from homolocus.coordinates import strand_char, to_one_based, to_zero_based


def test_zero_based_half_open_becomes_one_based_inclusive() -> None:
    assert to_one_based(1_900_000, 1_901_200) == (1_900_001, 1_901_200)


def test_one_based_inclusive_becomes_zero_based_half_open() -> None:
    assert to_zero_based(1_900_001, 1_901_200) == (1_900_000, 1_901_200)


def test_single_base_interval() -> None:
    assert to_one_based(0, 1) == (1, 1)
    assert to_zero_based(1, 1) == (0, 1)


@pytest.mark.parametrize(("start", "end"), [(5, 5), (6, 5), (-1, 5)])
def test_to_one_based_rejects_empty_or_negative_intervals(start: int, end: int) -> None:
    with pytest.raises(ValueError, match="interval"):
        to_one_based(start, end)


@pytest.mark.parametrize(("start", "end"), [(6, 5), (0, 5), (-2, 5)])
def test_to_zero_based_rejects_empty_or_non_positive_intervals(start: int, end: int) -> None:
    with pytest.raises(ValueError, match="interval"):
        to_zero_based(start, end)


@given(start=st.integers(min_value=0, max_value=10**9), length=st.integers(1, 10**6))
def test_round_trip_preserves_interval(start: int, length: int) -> None:
    end = start + length
    assert to_zero_based(*to_one_based(start, end)) == (start, end)


@pytest.mark.parametrize(("strand", "char"), [(1, "+"), (-1, "-")])
def test_strand_char(strand: int, char: str) -> None:
    assert strand_char(strand) == char


@pytest.mark.parametrize("strand", [0, 2, -2])
def test_strand_char_rejects_invalid_strand(strand: int) -> None:
    with pytest.raises(ValueError, match="strand"):
        strand_char(strand)
