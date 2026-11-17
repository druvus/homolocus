"""Score components for a candidate (spec 8). Pure functions.

Each represented family contributes one pairing: its best reference member
and that member's best copy in the candidate (spec 7.7). Order and
orientation are reversal-invariant.
"""

from bisect import bisect_right
from collections.abc import Sequence
from dataclasses import dataclass

from homolocus.models.cluster import ScoreComponents
from homolocus.models.gene import Gene, genomic_order


@dataclass(frozen=True)
class FamilyPairing:
    """A represented family's best member, its best copy, and their evidence."""

    member: Gene
    copy: Gene
    identity: float
    coverage: float


def components(
    pairings: Sequence[FamilyPairing], *, total_families: int, candidate_span: int
) -> ScoreComponents:
    """The five components for one candidate (spec 8).

    With fewer than two represented families, order, orientation and spacing
    are 0, and similarity is 0 when none is represented.
    """
    ordered = sorted(pairings, key=lambda pairing: genomic_order(pairing.member))
    structured = len(ordered) >= 2
    return ScoreComponents(
        content=len(ordered) / total_families,
        similarity=_similarity(ordered),
        order=_order(ordered) if structured else 0.0,
        orientation=_orientation(ordered) if structured else 0.0,
        spacing=_spacing(ordered, candidate_span) if structured else 0.0,
    )


def _similarity(pairings: Sequence[FamilyPairing]) -> float:
    if not pairings:
        return 0.0
    values = [pairing.identity / 100 * pairing.coverage for pairing in pairings]
    return sum(values) / len(values)


def _order(pairings: Sequence[FamilyPairing]) -> float:
    positions = [genomic_order(pairing.copy) for pairing in pairings]
    increasing = _longest_non_decreasing(positions)
    decreasing = _longest_non_decreasing([tuple(-v for v in key) for key in positions])
    return max(increasing, decreasing) / len(pairings)


def _longest_non_decreasing(keys: Sequence[tuple[int, ...]]) -> int:
    tails: list[tuple[int, ...]] = []
    for key in keys:
        slot = bisect_right(tails, key)
        if slot == len(tails):
            tails.append(key)
        else:
            tails[slot] = key
    return len(tails)


def _orientation(pairings: Sequence[FamilyPairing]) -> float:
    neighbours = list(zip(pairings, pairings[1:], strict=False))
    agreeing = sum(
        (left.member.strand == right.member.strand) == (left.copy.strand == right.copy.strand)
        for left, right in neighbours
    )
    return agreeing / len(neighbours)


def _spacing(pairings: Sequence[FamilyPairing], candidate_span: int) -> float:
    reference_span = max(p.member.end for p in pairings) - min(p.member.start for p in pairings)
    return min(candidate_span, reference_span) / max(candidate_span, reference_span)
