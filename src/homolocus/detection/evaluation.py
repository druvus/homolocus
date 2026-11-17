"""Candidate evaluation: coverage, family evidence, filtering and ranking (spec 7.7, 7.9, 8)."""

from collections.abc import Iterable, Mapping

from homolocus.intervals import union_length
from homolocus.models.cluster import FamilyEvidence
from homolocus.models.gene_family import GeneFamily
from homolocus.models.gene_match import GeneMatch


def family_evidence(
    family: GeneFamily,
    matches: Iterable[GeneMatch],
    *,
    target_order: Mapping[str, int],
    min_reference_gene_coverage: float,
) -> FamilyEvidence:
    """Support for one family from a candidate's passing evidence matches (spec 7.7).

    ``matches`` are the candidate's GeneMatches that pass the pair-level
    thresholds; ``target_order`` gives each evidence gene's genomic position.
    A member's coverage is the union of aligned reference positions over all
    its matches. The best member has the highest coverage (ties: genomic
    order); the best copy is its match with the highest bit score (ties:
    target order). The family is represented when that coverage reaches
    ``min_reference_gene_coverage``.
    """
    passing = tuple(matches)
    by_member = {
        member: [m for m in passing if m.reference_gene_id == member] for member in family.gene_ids
    }
    supported = [member for member in family.gene_ids if by_member[member]]
    copies = sorted(
        {m.target_gene_id for ms in by_member.values() for m in ms}, key=target_order.__getitem__
    )
    if not supported:
        return FamilyEvidence(family.id, 0.0, False, None, None, ())
    coverage = {member: _coverage(by_member[member]) for member in supported}
    best_member = max(supported, key=lambda member: coverage[member])
    best_copy = min(
        by_member[best_member],
        key=lambda m: (-m.bit_score, target_order[m.target_gene_id]),
    ).target_gene_id
    return FamilyEvidence(
        family_id=family.id,
        coverage=coverage[best_member],
        represented=coverage[best_member] >= min_reference_gene_coverage,
        best_member_id=best_member,
        best_copy_id=best_copy,
        copies=tuple(copies),
    )


def _coverage(matches: list[GeneMatch]) -> float:
    intervals = [hsp.query_interval for match in matches for hsp in match.alignments]
    return union_length(intervals) / matches[0].query_length
