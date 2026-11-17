"""Conversion between internal and user-facing coordinates (spec 4).

Internally, intervals are 0-based and half-open: ``[start, end)``.
Everything a user types or reads is 1-based and inclusive. This module is
the only place where one is converted into the other.
"""

_STRAND_CHARS = {1: "+", -1: "-"}


def to_one_based(start: int, end: int) -> tuple[int, int]:
    """Convert a 0-based half-open interval to a 1-based inclusive one."""
    if start < 0 or start >= end:
        raise ValueError(f"invalid 0-based half-open interval [{start}, {end})")
    return start + 1, end


def to_zero_based(start: int, end: int) -> tuple[int, int]:
    """Convert a 1-based inclusive interval to a 0-based half-open one."""
    if start < 1 or start > end:
        raise ValueError(f"invalid 1-based inclusive interval {start}..{end}")
    return start - 1, end


def strand_char(strand: int) -> str:
    """Return ``+`` for strand +1 and ``-`` for strand -1."""
    try:
        return _STRAND_CHARS[strand]
    except KeyError:
        raise ValueError(f"invalid strand {strand!r}; expected 1 or -1") from None
