"""Candidate evaluation: coverage, family evidence, filtering and ranking (spec 7.7, 7.9, 8)."""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from homolocus.config import DiamondSearch, RegionDetection, ScoringWeights
from homolocus.detection.candidate_generation import EventWindow, padded_interval, windows
from homolocus.detection.flags import fusion_events, split_events, touches_contig_edge
from homolocus.detection.matching import matching_events
from homolocus.detection.scoring import FamilyPairing, components, is_reversed
from homolocus.intervals import union_length
from homolocus.models.cluster import Cluster, FamilyEvidence
from homolocus.models.contig import Contig
from homolocus.models.gene import Gene
from homolocus.models.gene_family import GeneFamily, ReferenceFamilies
from homolocus.models.gene_match import GeneMatch
from homolocus.models.genome import Genome
from homolocus.models.reference import ReferenceRegion

SCORE_DECIMALS = 12
"""Scores are compared at this precision when ranking, so floating-point noise
never decides a tie that spec 7.9 resolves by reference coverage."""


@dataclass(frozen=True)
class DetectionContext:
    """Everything detection needs besides the target genome and its matches.

    Construction checks that ``min_matches`` does not exceed the number of
    families and that every family member is a reference gene (spec 7.5).
    """

    reference: ReferenceRegion
    families: ReferenceFamilies
    search: DiamondSearch
    detection: RegionDetection
    weights: ScoringWeights

    def __post_init__(self) -> None:
        check_min_matches(self.families, self.detection.min_matches)
        reference_ids = {gene.id for gene in self.reference.genes}
        for family in self.families:
            for gene_id in family.gene_ids:
                if gene_id not in reference_ids:
                    raise ValueError(f"family member {gene_id} is not a gene of the reference")


def check_min_matches(families: ReferenceFamilies, min_matches: int) -> None:
    """Raise ValueError when no candidate could represent ``min_matches`` families (spec 7.5).

    A reference without searchable genes has no families at all, which is
    reported as such rather than as a min_matches problem.
    """
    if not len(families):
        raise ValueError(
            "the reference region has no searchable protein-coding genes; "
            "check the region coordinates and its CDS annotation"
        )
    if min_matches > len(families):
        raise ValueError(
            f"min_matches {min_matches} exceeds the {len(families)} reference families; "
            "lower it with --min-matches (RegionDetection.min_matches)"
        )


def detect(
    context: DetectionContext, genome: Genome, matches: Iterable[GeneMatch]
) -> tuple[Cluster, ...]:
    """All passing candidates in the genome, ranked and deduplicated (spec 7, 8).

    Each contig is processed on its own: matching events, windows,
    evaluation, then the min_matches and min_reference_coverage filters
    (both must hold). Candidates are ranked by score, reference coverage,
    contig order, start, end and event IDs; a candidate overlapping a kept,
    higher-ranked one on the same contig by more than max_candidate_overlap
    of the shorter span is dropped.
    """
    all_matches = tuple(matches)
    found = [
        cluster
        for contig in genome.contigs.values()
        for cluster in _contig_candidates(context, contig, all_matches)
        if _passes(cluster, context.detection)
    ]
    contig_rank = {contig_id: rank for rank, contig_id in enumerate(genome.contigs)}
    ranked = sorted(
        found,
        key=lambda c: (
            -round(c.score, SCORE_DECIMALS),
            -c.reference_coverage,
            contig_rank[c.contig_id],
            c.start,
            c.end,
            c.event_ids,
        ),
    )
    return _deduplicate(ranked, context.detection.max_candidate_overlap)


def _contig_candidates(
    context: DetectionContext, contig: Contig, matches: Sequence[GeneMatch]
) -> list[Cluster]:
    gene_ids = {gene.id for gene in contig.genes}
    on_contig = [match for match in matches if match.target_gene_id in gene_ids]
    events = matching_events(
        contig,
        on_contig,
        max_evalue=context.search.evalue,
        min_identity=context.search.min_identity,
        min_event_coverage=context.detection.min_event_coverage,
    )
    return [
        evaluate_window(context, window, contig, on_contig)
        for window in windows(contig, events, context.detection)
    ]


def _passes(cluster: Cluster, detection: RegionDetection) -> bool:
    enough = len(cluster.represented_families) >= detection.min_matches
    return enough and cluster.reference_coverage >= detection.min_reference_coverage


def _deduplicate(ranked: Sequence[Cluster], max_overlap: float) -> tuple[Cluster, ...]:
    kept: list[Cluster] = []
    for cluster in ranked:
        if not any(_overlap(cluster, other) > max_overlap for other in kept):
            kept.append(cluster)
    return tuple(kept)


def _overlap(first: Cluster, second: Cluster) -> float:
    """Shared span as a fraction of the shorter span; 0 on different contigs."""
    if first.contig_id != second.contig_id:
        return 0.0
    shared = min(first.end, second.end) - max(first.start, second.start)
    return max(shared, 0) / min(first.span, second.span)


def evaluate_window(
    context: DetectionContext,
    window: EventWindow,
    contig: Contig,
    matches: Iterable[GeneMatch],
) -> Cluster:
    """Evaluate one candidate window into a Cluster (spec 7.7, 8, 10).

    Evidence genes lie fully inside the unpadded span; all their matches are
    kept, and those passing the pair-level thresholds count for coverage.
    """
    detection = context.detection
    display = padded_interval(
        window, contig, padding_genes=detection.padding_genes, padding_bp=detection.padding_bp
    )
    genes = tuple(g for g in contig.genes if display[0] <= g.start and g.end <= display[1])
    evidence = [g for g in genes if window.start <= g.start and g.end <= window.end]
    evidence_ids = {gene.id for gene in evidence}
    kept = tuple(match for match in matches if match.target_gene_id in evidence_ids)
    passing = [
        m
        for m in kept
        if m.passes(max_evalue=context.search.evalue, min_identity=context.search.min_identity)
    ]
    target_order = {gene.id: rank for rank, gene in enumerate(evidence)}
    families = tuple(
        family_evidence(
            family,
            passing,
            target_order=target_order,
            min_reference_gene_coverage=detection.min_reference_gene_coverage,
        )
        for family in context.families
    )
    pairings = _pairings(context.reference, families, evidence, passing)
    return Cluster(
        contig_id=contig.id,
        start=window.start,
        end=window.end,
        display_start=display[0],
        display_end=display[1],
        genes=genes,
        event_ids=window.event_ids,
        matches=kept,
        families=families,
        components=components(pairings, total_families=len(families), candidate_span=window.span),
        weights=context.weights,
        touches_contig_edge=touches_contig_edge(
            window, contig, edge_margin_bp=detection.edge_margin_bp
        ),
        splits=split_events(contig, passing, detection),
        fusions=fusion_events(contig, passing, context.families, detection),
        reversed=is_reversed(pairings),
    )


def _pairings(
    reference: ReferenceRegion,
    families: Sequence[FamilyEvidence],
    evidence: Sequence[Gene],
    passing: Sequence[GeneMatch],
) -> list[FamilyPairing]:
    members = {gene.id: gene for gene in reference.genes}
    copies = {gene.id: gene for gene in evidence}
    by_pair = {(m.reference_gene_id, m.target_gene_id): m for m in passing}
    pairings = []
    for family in families:
        if family.represented and family.best_member_id and family.best_copy_id:
            match = by_pair[family.best_member_id, family.best_copy_id]
            pairings.append(
                FamilyPairing(
                    member=members[family.best_member_id],
                    copy=copies[family.best_copy_id],
                    identity=match.identity,
                    coverage=family.coverage,
                )
            )
    return pairings


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
        return FamilyEvidence(family.id, 0.0, False, None, None, (), family.name)
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
        family_name=family.name,
    )


def _coverage(matches: list[GeneMatch]) -> float:
    intervals = [hsp.query_interval for match in matches for hsp in match.alignments]
    return union_length(intervals) / matches[0].query_length
