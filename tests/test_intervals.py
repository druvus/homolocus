"""Interval unions used for coverage (spec 5.7, 7.7)."""

import pytest
from hypothesis import given
from hypothesis import strategies as st

from homolocus.intervals import merge_intervals, union_length


@pytest.mark.parametrize(
    ("intervals", "merged"),
    [
        ([], ()),
        ([(0, 10)], ((0, 10),)),
        ([(0, 10), (20, 30)], ((0, 10), (20, 30))),
        ([(0, 10), (5, 15)], ((0, 15),)),
        ([(0, 30), (5, 10)], ((0, 30),)),
        ([(0, 10), (10, 20)], ((0, 20),)),
        ([(20, 30), (0, 10), (8, 12)], ((0, 12), (20, 30))),
    ],
)
def test_merge_intervals(intervals, merged) -> None:
    assert merge_intervals(intervals) == merged


def test_union_length_counts_overlaps_once() -> None:
    assert union_length([(0, 10), (5, 15), (20, 25)]) == 20


def test_union_length_of_nothing_is_zero() -> None:
    assert union_length([]) == 0


@pytest.mark.parametrize("bad", [(5, 5), (6, 5), (-1, 5)])
def test_invalid_intervals_are_rejected(bad: tuple[int, int]) -> None:
    with pytest.raises(ValueError, match="interval"):
        merge_intervals([(0, 10), bad])


intervals_strategy = st.lists(
    st.tuples(st.integers(0, 200), st.integers(1, 50)).map(lambda p: (p[0], p[0] + p[1])),
    max_size=20,
)


@given(intervals_strategy)
def test_union_length_matches_brute_force(intervals: list[tuple[int, int]]) -> None:
    covered = {position for start, end in intervals for position in range(start, end)}
    assert union_length(intervals) == len(covered)


@given(intervals_strategy)
def test_union_length_is_bounded(intervals: list[tuple[int, int]]) -> None:
    lengths = [end - start for start, end in intervals]
    assert max(lengths, default=0) <= union_length(intervals) <= sum(lengths)


@given(intervals_strategy)
def test_merged_intervals_are_sorted_and_disjoint(intervals: list[tuple[int, int]]) -> None:
    merged = merge_intervals(intervals)
    assert all(left[1] < right[0] for left, right in zip(merged, merged[1:], strict=False))
