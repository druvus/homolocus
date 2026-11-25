"""Run configuration (spec 6.2). Values are provisional and not validated on real data."""

from dataclasses import dataclass
from typing import Literal, get_args

Sensitivity = Literal[
    "fast",
    "default",
    "mid-sensitive",
    "sensitive",
    "more-sensitive",
    "very-sensitive",
    "ultra-sensitive",
]
SENSITIVITIES: tuple[Sensitivity, ...] = get_args(Sensitivity)


@dataclass(frozen=True)
class DiamondSearch:
    """DIAMOND search settings and pair-level thresholds (spec 6.2).

    ``evalue`` and ``min_identity`` (percent) are pair-level thresholds applied
    to aggregated GeneMatch evidence. ``hsp_evalue`` is the looser per-HSP
    cutoff passed to DIAMOND. ``max_target_seqs`` 0 means unlimited and
    ``max_hsps`` 0 reports all HSPs per pair. ``threads`` None leaves DIAMOND's
    default. ``timeout`` is how many seconds each DIAMOND command may run before
    it is stopped; None sets no limit.
    """

    evalue: float = 1e-5
    min_identity: float = 30.0
    sensitivity: Sensitivity = "sensitive"
    max_target_seqs: int = 100
    max_hsps: int = 0
    hsp_evalue: float = 1e-3
    threads: int | None = None
    timeout: float | None = 3600.0

    def __post_init__(self) -> None:
        if self.evalue < 0:
            raise ValueError(f"evalue {self.evalue} must not be negative")
        if self.hsp_evalue < self.evalue:
            raise ValueError(f"hsp_evalue {self.hsp_evalue} is stricter than evalue {self.evalue}")
        if not 0 <= self.min_identity <= 100:
            raise ValueError(f"min_identity {self.min_identity} is not a percentage")
        if self.sensitivity not in SENSITIVITIES:
            raise ValueError(
                f"unknown sensitivity {self.sensitivity!r}; expected one of {SENSITIVITIES}"
            )
        self._check_counts()

    def _check_counts(self) -> None:
        if self.max_target_seqs < 0:
            raise ValueError(f"max_target_seqs {self.max_target_seqs} must not be negative")
        if self.max_hsps < 0:
            raise ValueError(f"max_hsps {self.max_hsps} must not be negative")
        if self.threads is not None and self.threads < 1:
            raise ValueError(f"threads {self.threads} must be at least 1")
        if self.timeout is not None and self.timeout <= 0:
            raise ValueError(f"timeout {self.timeout} must be positive")


@dataclass(frozen=True)
class FamilyGrouping:
    """Paralog family thresholds for the reference self-search (spec 7.6).

    ``family_min_identity`` is a percentage; None disables grouping so every
    gene is its own family. ``family_min_coverage`` is a fraction that both
    proteins must reach.
    """

    family_min_identity: float | None = 90.0
    family_min_coverage: float = 0.8

    def __post_init__(self) -> None:
        identity = self.family_min_identity
        if identity is not None and not 0 <= identity <= 100:
            raise ValueError(f"family_min_identity {identity} is not a percentage")
        if not 0 <= self.family_min_coverage <= 1:
            raise ValueError(f"family_min_coverage {self.family_min_coverage} is not a fraction")

    @property
    def enabled(self) -> bool:
        """Whether paralogs are grouped at all."""
        return self.family_min_identity is not None


_DETECTION_FRACTIONS = (
    "min_reference_coverage",
    "min_reference_gene_coverage",
    "min_event_coverage",
    "max_candidate_overlap",
    "min_split_contribution",
    "max_fragment_overlap",
    "min_fusion_part_coverage",
)
_DETECTION_NON_NEGATIVE = (
    "max_intervening_genes",
    "max_inter_match_bp_gap",
    "padding_genes",
    "padding_bp",
    "edge_margin_bp",
)


@dataclass(frozen=True)
class RegionDetection:
    """Candidate generation and detection thresholds (spec 7.5, 10).

    Coverages, ``max_candidate_overlap`` and the split and fusion settings
    (spec 9.4) are fractions. ``min_matches`` counts distinct represented
    families and must be at least 2. Distances are in base pairs or genes;
    ``max_region_bp`` applies to the unpadded span.
    """

    min_matches: int = 3
    min_reference_coverage: float = 0.3
    min_reference_gene_coverage: float = 0.5
    min_event_coverage: float = 0.5
    max_intervening_genes: int = 10
    max_inter_match_bp_gap: int = 20_000
    max_region_bp: int = 50_000
    max_candidate_overlap: float = 0.5
    padding_genes: int = 0
    padding_bp: int = 0
    edge_margin_bp: int = 500
    min_split_contribution: float = 0.2
    max_fragment_overlap: float = 0.2
    min_fusion_part_coverage: float = 0.2

    def __post_init__(self) -> None:
        if self.min_matches < 2:
            raise ValueError(f"min_matches {self.min_matches} must be at least 2")
        if self.max_region_bp < 1:
            raise ValueError(f"max_region_bp {self.max_region_bp} must be at least 1")
        for name in _DETECTION_FRACTIONS:
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} {getattr(self, name)} is not a fraction")
        for name in _DETECTION_NON_NEGATIVE:
            if getattr(self, name) < 0:
                raise ValueError(f"{name} {getattr(self, name)} must not be negative")


_COMPONENTS = ("content", "similarity", "order", "orientation", "spacing")


@dataclass(frozen=True)
class ScoringWeights:
    """Weights of the five score components (spec 8). Unvalidated starting values.

    Weights must be non-negative with a positive sum; they are normalized to
    sum to 1 when the score is computed.
    """

    content: float = 0.40
    similarity: float = 0.30
    order: float = 0.15
    orientation: float = 0.10
    spacing: float = 0.05

    def __post_init__(self) -> None:
        for name in _COMPONENTS:
            if getattr(self, name) < 0:
                raise ValueError(f"weight {name} {getattr(self, name)} must not be negative")
        if self._total() <= 0:
            raise ValueError("at least one weight must be positive")

    def normalized(self) -> "ScoringWeights":
        """The same weights scaled to sum to 1."""
        total = self._total()
        return ScoringWeights(**{name: getattr(self, name) / total for name in _COMPONENTS})

    def _total(self) -> float:
        return float(sum(getattr(self, name) for name in _COMPONENTS))


@dataclass(frozen=True)
class PlotOptions:
    """What the written overview plots show (spec 11.1, 11.4).

    ``candidates`` is the most candidate lines drawn per genome, best-ranked
    first; 0 draws all. An overview longer than ``lines_per_page`` lines is
    split into pages, each starting with the reference line.
    """

    candidates: int = 3
    lines_per_page: int = 30

    def __post_init__(self) -> None:
        if self.candidates < 0:
            raise ValueError(f"candidates {self.candidates} must be 0 (all) or more")
        if self.lines_per_page < 2:
            raise ValueError(f"lines_per_page {self.lines_per_page} must be at least 2")


DEFAULT_PLOTS = PlotOptions()
"""The plot options used when none are given."""
