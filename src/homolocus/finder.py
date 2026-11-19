"""RegionFinder: run the whole pipeline for one reference region and one genome (spec 12)."""

import tempfile
import warnings
from collections.abc import Callable, Iterable
from dataclasses import dataclass, replace
from functools import partial
from pathlib import Path
from typing import TypeVar

from homolocus.config import DiamondSearch, FamilyGrouping, RegionDetection, ScoringWeights
from homolocus.detection.evaluation import DetectionContext, detect
from homolocus.detection.families import group_families
from homolocus.models.alignment import Alignment
from homolocus.models.gene_family import ReferenceFamilies
from homolocus.models.gene_match import aggregate
from homolocus.models.genome import Genome
from homolocus.models.reference import ReferenceRegion
from homolocus.results import (
    GenomeSetResults,
    RegionResults,
    check_genome_labels,
    check_genome_names,
)
from homolocus.search.diamond import diamond_version, search
from homolocus.search.protein_index import ProteinIndex


@dataclass(frozen=True)
class RegionFinder:
    """Find regions homologous to ``reference`` in target genomes (spec 12).

    Settings default to the provisional values of the specification.
    """

    reference: ReferenceRegion
    search: DiamondSearch = DiamondSearch()
    family_grouping: FamilyGrouping = FamilyGrouping()
    detection: RegionDetection = RegionDetection()
    scoring: ScoringWeights = ScoringWeights()

    def find(self, genome: Genome) -> RegionResults:
        """Search, group families, detect and rank candidates in ``genome``.

        DIAMOND runs in a temporary directory; nothing is written elsewhere.
        Warnings raised during the run are kept in the results and re-issued.
        """
        shared, early = _recorded(self._shared)
        results, late = _recorded(lambda: self._run(genome, shared))
        return replace(results, warnings=early + late)

    def find_all(self, genomes: Iterable[Genome]) -> GenomeSetResults:
        """Find regions in several genomes, one search each (spec 2.2, 12).

        The self-search and families are computed once and shared, so each
        genome's results equal ``find`` on that genome alone. Genome names
        are checked before any search.
        """
        targets = tuple(genomes)
        check_genome_names(genome.name for genome in targets)
        check_genome_labels(genome.label for genome in targets)
        shared, early = _recorded(self._shared)
        runs = [_recorded(partial(self._run, genome, shared)) for genome in targets]
        return GenomeSetResults(
            tuple(replace(results, warnings=early + late) for results, late in runs)
        )

    def _shared(self) -> "_Shared":
        """What every genome shares: the DIAMOND version, families and self-search."""
        version = diamond_version()
        index = ProteinIndex.build(self.reference.genes, ())
        with tempfile.TemporaryDirectory(prefix="homolocus-") as workdir:
            self_hits, raw_self = self._self_search(index, Path(workdir) / "reference_self.tsv")
        families = group_families(self.reference.genes, self_hits, self.family_grouping)
        return _Shared(version, families, raw_self)

    def _run(self, genome: Genome, shared: "_Shared") -> RegionResults:
        index = ProteinIndex.build(self.reference.genes, genome.genes())
        with tempfile.TemporaryDirectory(prefix="homolocus-") as workdir:
            main_path = Path(workdir) / "diamond.tsv"
            matches = aggregate(search(index, self.search, main_path))
            raw_main = _read_unchanged(main_path)
        context = DetectionContext(
            self.reference, shared.families, self.search, self.detection, self.scoring
        )
        return RegionResults(
            reference=self.reference,
            genome=genome,
            families=shared.families,
            matches=matches,
            candidates=detect(context, genome, matches),
            search=self.search,
            family_grouping=self.family_grouping,
            detection=self.detection,
            scoring=self.scoring,
            diamond_version=shared.version,
            raw_main=raw_main,
            raw_self=shared.raw_self,
            protein_index=index,
        )

    def _self_search(self, index: ProteinIndex, path: Path) -> tuple[tuple[Alignment, ...], str]:
        if not self.family_grouping.enabled:
            return (), ""
        hits = search(index, self.search, path, self_search=True)
        return hits, _read_unchanged(path)


@dataclass(frozen=True)
class _Shared:
    version: str
    families: ReferenceFamilies
    raw_self: str


T = TypeVar("T")


def _recorded(action: Callable[[], T]) -> tuple[T, tuple[str, ...]]:
    """Run ``action``, keep the warnings it raises as text, and re-issue them."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        value = action()
    for warning in caught:
        warnings.warn_explicit(warning.message, warning.category, warning.filename, warning.lineno)
    return value, tuple(str(w.message) for w in caught)


def _read_unchanged(path: Path) -> str:
    """DIAMOND's output as text, without newline translation (works on Python 3.10+)."""
    return path.read_bytes().decode("utf-8")
