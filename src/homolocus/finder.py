"""RegionFinder: run the whole pipeline for one reference region and one genome (spec 12)."""

import tempfile
import warnings
from dataclasses import dataclass, replace
from pathlib import Path

from homolocus.config import DiamondSearch, FamilyGrouping, RegionDetection, ScoringWeights
from homolocus.detection.evaluation import DetectionContext, detect
from homolocus.detection.families import group_families
from homolocus.models.alignment import Alignment
from homolocus.models.gene_match import aggregate
from homolocus.models.genome import Genome
from homolocus.models.reference import ReferenceRegion
from homolocus.results import RegionResults
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
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            results = self._run(genome)
        for warning in caught:
            warnings.warn_explicit(
                warning.message, warning.category, warning.filename, warning.lineno
            )
        return replace(results, warnings=tuple(str(w.message) for w in caught))

    def _run(self, genome: Genome) -> RegionResults:
        index = ProteinIndex.build(self.reference.genes, genome.genes())
        version = diamond_version()
        with tempfile.TemporaryDirectory(prefix="homolocus-") as workdir:
            self_hits, raw_self = self._self_search(index, Path(workdir) / "reference_self.tsv")
            main_path = Path(workdir) / "diamond.tsv"
            matches = aggregate(search(index, self.search, main_path))
            raw_main = _read_unchanged(main_path)
        families = group_families(self.reference.genes, self_hits, self.family_grouping)
        context = DetectionContext(
            self.reference, families, self.search, self.detection, self.scoring
        )
        return RegionResults(
            reference=self.reference,
            genome=genome,
            families=families,
            matches=matches,
            candidates=detect(context, genome, matches),
            search=self.search,
            family_grouping=self.family_grouping,
            detection=self.detection,
            scoring=self.scoring,
            diamond_version=version,
            raw_main=raw_main,
            raw_self=raw_self,
            protein_index=index,
        )

    def _self_search(self, index: ProteinIndex, path: Path) -> tuple[tuple[Alignment, ...], str]:
        if not self.family_grouping.enabled:
            return (), ""
        hits = search(index, self.search, path, self_search=True)
        return hits, _read_unchanged(path)


def _read_unchanged(path: Path) -> str:
    """DIAMOND's output as text, without newline translation (works on Python 3.10+)."""
    return path.read_bytes().decode("utf-8")
