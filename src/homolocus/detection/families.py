"""Group paralogous reference genes into families (spec 7.6)."""

from collections.abc import Iterable, Sequence

from homolocus.config import FamilyGrouping
from homolocus.models.alignment import Alignment
from homolocus.models.gene import Gene, genomic_order
from homolocus.models.gene_family import GeneFamily, ReferenceFamilies
from homolocus.models.gene_match import GeneMatch, aggregate


def group_families(
    reference_genes: Iterable[Gene],
    self_alignments: Iterable[Alignment],
    config: FamilyGrouping,
) -> ReferenceFamilies:
    """Families of searchable reference genes from the self-search alignments.

    Two genes are linked when either direction of their aggregated evidence
    reaches the identity threshold and both coverage thresholds. Families are
    connected components, numbered F001, F002, ... in genomic order of their
    first member. Genes without a protein are left out.
    """
    genes = sorted(
        (g for g in reference_genes if g.protein_sequence is not None), key=genomic_order
    )
    parent = {gene.id: gene.id for gene in genes}
    if config.family_min_identity is not None:
        for match in aggregate(self_alignments):
            if _links(match, config.family_min_identity, config.family_min_coverage, parent):
                _union(parent, match.reference_gene_id, match.target_gene_id)
    components: dict[str, list[str]] = {}
    for gene in genes:
        components.setdefault(_find(parent, gene.id), []).append(gene.id)
    by_id = {gene.id: gene for gene in genes}
    families = (
        GeneFamily(f"F{number:03d}", tuple(members), family_name([by_id[m] for m in members]))
        for number, members in enumerate(components.values(), start=1)
    )
    return ReferenceFamilies(tuple(families))


def family_name(members: Sequence[Gene]) -> str | None:
    """Name of a family from its members' annotation (spec 7.6).

    The ``/gene`` names, unique and in member order, joined with ``/``; when no
    member has one, the ``/product`` descriptions the same way; else None.
    """
    for values in ([m.gene_name for m in members], [m.product for m in members]):
        unique = list(dict.fromkeys(value for value in values if value))
        if unique:
            return "/".join(unique)
    return None


def _links(match: GeneMatch, identity: float, coverage: float, parent: dict[str, str]) -> bool:
    known = match.reference_gene_id in parent and match.target_gene_id in parent
    distinct = match.reference_gene_id != match.target_gene_id
    covered = match.query_coverage >= coverage and match.subject_coverage >= coverage
    return known and distinct and covered and match.identity >= identity


def _find(parent: dict[str, str], gene_id: str) -> str:
    root = gene_id
    while parent[root] != root:
        root = parent[root]
    while parent[gene_id] != root:
        parent[gene_id], gene_id = root, parent[gene_id]
    return root


def _union(parent: dict[str, str], first: str, second: str) -> None:
    parent[_find(parent, second)] = _find(parent, first)
