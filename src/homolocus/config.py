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
    default.
    """

    evalue: float = 1e-5
    min_identity: float = 30.0
    sensitivity: Sensitivity = "sensitive"
    max_target_seqs: int = 100
    max_hsps: int = 0
    hsp_evalue: float = 1e-3
    threads: int | None = None

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
