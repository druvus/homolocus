"""Candidate result models (spec 5.8, 8, 8.1).

A Cluster is the result of evaluating one candidate window. It keeps the
evidence it was built from; the score and coverage metrics are derived, so
they cannot disagree with that evidence.
"""

from dataclasses import dataclass, fields
from typing import Literal

from homolocus.config import ScoringWeights
from homolocus.models.gene import Gene
from homolocus.models.gene_match import GeneMatch

Role = Literal["matched", "additional"]


@dataclass(frozen=True)
class ScoreComponents:
    """The five score components, each in [0, 1] (spec 8)."""

    content: float
    similarity: float
    order: float
    orientation: float
    spacing: float

    def __post_init__(self) -> None:
        for field in fields(self):
            value = getattr(self, field.name)
            if not 0 <= value <= 1:
                raise ValueError(f"score component {field.name} {value} is not in [0, 1]")


@dataclass(frozen=True)
class FamilyEvidence:
    """How one reference family is supported in a candidate (spec 7.4, 7.7).

    ``copies`` are the evidence genes with a passing match to any member;
    ``best_member_id`` and ``best_copy_id`` form the pairing used for scoring.
    """

    family_id: str
    coverage: float
    represented: bool
    best_member_id: str | None
    best_copy_id: str | None
    copies: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.family_id:
            raise ValueError("family ID must not be empty")
        if not 0 <= self.coverage <= 1:
            raise ValueError(f"family {self.family_id} coverage {self.coverage} is not a fraction")
        if len(set(self.copies)) != len(self.copies):
            raise ValueError(f"family {self.family_id} lists duplicate copies")
        self._check_pairing()

    def _check_pairing(self) -> None:
        if self.best_copy_id is not None and self.best_copy_id not in self.copies:
            raise ValueError(f"family {self.family_id} best copy is not one of its copies")
        if self.best_copy_id is not None and self.best_member_id is None:
            raise ValueError(f"family {self.family_id} has a best copy but no best member")
        if self.represented and self.best_copy_id is None:
            raise ValueError(f"represented family {self.family_id} needs a best copy")
        if self.best_member_id is None and self.coverage != 0:
            raise ValueError(f"family {self.family_id} has coverage without a best member")

    @property
    def copy_count(self) -> int:
        """Number of duplicated homologue copies in the candidate."""
        return len(self.copies)


@dataclass(frozen=True)
class Cluster:
    """One evaluated candidate region on one target contig (spec 5.8).

    ``start``/``end`` is the unpadded span (0-based half-open) and
    ``display_start``/``display_end`` the padded interval. ``genes`` are the
    target genes in the padded interval; evidence genes lie fully inside the
    unpadded span. ``matches`` holds every GeneMatch of an evidence gene.
    """

    contig_id: str
    start: int
    end: int
    display_start: int
    display_end: int
    genes: tuple[Gene, ...]
    event_ids: tuple[str, ...]
    matches: tuple[GeneMatch, ...]
    families: tuple[FamilyEvidence, ...]
    components: ScoreComponents
    weights: ScoringWeights
    touches_contig_edge: bool

    def __post_init__(self) -> None:
        if not 0 <= self.start < self.end:
            raise ValueError(f"invalid candidate interval [{self.start}, {self.end})")
        if not self.display_start <= self.start or not self.end <= self.display_end:
            raise ValueError("display interval must contain the candidate span")
        if any(gene.contig_id != self.contig_id for gene in self.genes):
            raise ValueError(f"cluster genes must lie on contig {self.contig_id}")
        self._check_evidence()
        self._check_families()

    def _check_evidence(self) -> None:
        evidence = {gene.id for gene in self.evidence_genes}
        if not self.event_ids:
            raise ValueError("a cluster needs at least one event")
        if not set(self.event_ids) <= evidence:
            raise ValueError("every event must be an evidence gene of the cluster")
        if any(match.target_gene_id not in evidence for match in self.matches):
            raise ValueError("every match must belong to an evidence gene of the cluster")

    def _check_families(self) -> None:
        if not self.families:
            raise ValueError("a cluster needs at least one family")
        if len({family.family_id for family in self.families}) != len(self.families):
            raise ValueError("duplicate family in cluster")
        evidence = {gene.id for gene in self.evidence_genes}
        if any(not set(family.copies) <= evidence for family in self.families):
            raise ValueError("every family copy must be an evidence gene of the cluster")

    @property
    def span(self) -> int:
        """Unpadded length in base pairs."""
        return self.end - self.start

    @property
    def evidence_genes(self) -> tuple[Gene, ...]:
        """Genes lying fully inside the unpadded span (spec 7.7)."""
        return tuple(g for g in self.genes if self.start <= g.start and g.end <= self.end)

    def role(self, gene_id: str) -> Role:
        """``matched`` for matching events, ``additional`` for other cluster genes."""
        if gene_id not in {gene.id for gene in self.genes}:
            raise KeyError(gene_id)
        return "matched" if gene_id in self.event_ids else "additional"

    @property
    def represented_families(self) -> tuple[FamilyEvidence, ...]:
        """Families represented in the candidate (spec 7.7)."""
        return tuple(family for family in self.families if family.represented)

    @property
    def missing_families(self) -> tuple[FamilyEvidence, ...]:
        """Families not represented in the candidate (spec 8.2)."""
        return tuple(family for family in self.families if not family.represented)

    @property
    def reference_coverage(self) -> float:
        """Represented families divided by all families (spec 8.1)."""
        return len(self.represented_families) / len(self.families)

    @property
    def mean_protein_coverage(self) -> float:
        """Mean 7.7 coverage of the represented families; 0 when none (spec 8.1)."""
        represented = self.represented_families
        if not represented:
            return 0.0
        return sum(family.coverage for family in represented) / len(represented)

    @property
    def contains_partial_genes(self) -> bool:
        """Whether any evidence gene is partial."""
        return any(gene.partial for gene in self.evidence_genes)

    @property
    def score(self) -> float:
        """Weighted sum of the components with normalized weights (spec 8)."""
        weights, parts = self.weights.normalized(), self.components
        return (
            weights.content * parts.content
            + weights.similarity * parts.similarity
            + weights.order * parts.order
            + weights.orientation * parts.orientation
            + weights.spacing * parts.spacing
        )
