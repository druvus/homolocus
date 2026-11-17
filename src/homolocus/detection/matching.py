"""Matching events: target genes eligible to seed candidates (spec 7.2)."""

from collections.abc import Iterable

from homolocus.models.contig import Contig
from homolocus.models.gene import Gene
from homolocus.models.gene_match import GeneMatch


def matching_events(
    contig: Contig,
    matches: Iterable[GeneMatch],
    *,
    max_evalue: float,
    min_identity: float,
    min_event_coverage: float,
) -> tuple[Gene, ...]:
    """Genes of ``contig`` with at least one qualifying match, in genomic order.

    A match qualifies when it passes the pair-level thresholds and covers at
    least ``min_event_coverage`` of the reference or of the target protein.
    Each duplicated copy is its own event.
    """
    qualifying = {
        match.target_gene_id
        for match in matches
        if match.passes(max_evalue=max_evalue, min_identity=min_identity)
        and _covers(match, min_event_coverage)
    }
    return tuple(gene for gene in contig.genes if gene.id in qualifying)


def _covers(match: GeneMatch, min_event_coverage: float) -> bool:
    return (
        match.query_coverage >= min_event_coverage or match.subject_coverage >= min_event_coverage
    )
