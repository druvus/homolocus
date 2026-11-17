"""Unions of 0-based half-open intervals, used for coverage (spec 5.7, 7.7)."""

from collections.abc import Iterable

from homolocus.coordinates import to_one_based

Interval = tuple[int, int]


def merge_intervals(intervals: Iterable[Interval]) -> tuple[Interval, ...]:
    """Sorted union of intervals; overlapping and touching intervals are joined."""
    merged: list[Interval] = []
    for start, end in sorted(intervals):
        to_one_based(start, end)
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


def union_length(intervals: Iterable[Interval]) -> int:
    """Number of positions covered by at least one interval."""
    return sum(end - start for start, end in merge_intervals(intervals))
