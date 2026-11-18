"""Results of a RegionFinder run: candidates, summary and gene rows, data frames (spec 12, 14).

Row layouts live here once and feed both the data frames and the TSV files.
User-facing coordinates are 1-based inclusive (spec 4).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from homolocus.config import DiamondSearch, FamilyGrouping, RegionDetection, ScoringWeights
from homolocus.coordinates import strand_char, to_one_based
from homolocus.models.cluster import Cluster, FamilyEvidence
from homolocus.models.gene import Gene
from homolocus.models.gene_family import ReferenceFamilies
from homolocus.models.gene_match import GeneMatch
from homolocus.models.genome import Genome
from homolocus.models.reference import ReferenceRegion
from homolocus.search.protein_index import ProteinIndex

if TYPE_CHECKING:
    from matplotlib.figure import Figure

Row = dict[str, Any]

SUMMARY_COLUMNS = (
    "rank",
    "candidate_id",
    "contig",
    "start",
    "end",
    "display_start",
    "display_end",
    "span_bp",
    "score",
    "content",
    "similarity",
    "order",
    "orientation",
    "spacing",
    "reference_coverage",
    "mean_protein_coverage",
    "n_target_genes",
    "n_matched_genes",
    "n_missing_families",
    "flags",
)
GENE_COLUMNS = (
    "target_gene_id",
    "locus_tag",
    "product",
    "start",
    "end",
    "strand",
    "role",
    "reference_gene_id",
    "family_id",
    "family_name",
    "identity",
    "query_coverage",
    "subject_coverage",
    "bit_score",
    "is_best_copy",
    "flags",
)


@dataclass(frozen=True)
class RegionResults:
    """Everything a RegionFinder run produced (spec 12).

    ``candidates`` are ranked; ``raw_main`` and ``raw_self`` are DIAMOND's
    tabular output, unchanged, and ``protein_index`` maps their tokens.
    """

    reference: ReferenceRegion
    genome: Genome
    families: ReferenceFamilies
    matches: tuple[GeneMatch, ...]
    candidates: tuple[Cluster, ...]
    search: DiamondSearch
    family_grouping: FamilyGrouping
    detection: RegionDetection
    scoring: ScoringWeights
    diamond_version: str
    raw_main: str
    raw_self: str
    protein_index: ProteinIndex
    warnings: tuple[str, ...] = ()

    @property
    def best_candidate(self) -> Cluster | None:
        """The top-ranked candidate, or None when nothing passed."""
        return self.candidates[0] if self.candidates else None

    @staticmethod
    def candidate_id(rank: int) -> str:
        """ID of the candidate at a 1-based rank: ``candidate_001``, ..."""
        return f"candidate_{rank:03d}"

    def to_dataframe(self) -> pd.DataFrame:
        """One row per candidate, the ``summary.tsv`` columns (spec 12, 14)."""
        rows = [_summary_row(rank, cluster) for rank, cluster in enumerate(self.candidates, 1)]
        frame: pd.DataFrame = pd.DataFrame(rows, columns=list(SUMMARY_COLUMNS))
        return frame

    def genes_dataframe(self) -> pd.DataFrame:
        """Gene rows of every candidate, with ``candidate_id`` first (spec 12, 14)."""
        rows = [
            {"candidate_id": self.candidate_id(rank), **row}
            for rank, cluster in enumerate(self.candidates, 1)
            for row in self.gene_rows(cluster)
        ]
        frame: pd.DataFrame = pd.DataFrame(rows, columns=["candidate_id", *GENE_COLUMNS])
        return frame

    def gene_rows(self, cluster: Cluster) -> list[Row]:
        """Rows of one candidate TSV: per gene, one per GeneMatch or one empty row."""
        families = {family.family_id: family for family in cluster.families}
        rows = []
        for gene in cluster.genes:
            matches = sorted(
                (m for m in cluster.matches if m.target_gene_id == gene.id),
                key=lambda m: m.reference_gene_id,
            )
            base = _gene_columns(cluster, gene)
            if not matches:
                rows.append({**base, **_EMPTY_MATCH})
            for match in matches:
                family = families[self.families.family_of(match.reference_gene_id).id]
                rows.append({**base, **_match_columns(match, family)})
        return rows

    def write(self, output_dir: str | Path) -> None:
        """Write every output file of spec 14 to ``output_dir``."""
        from homolocus.output import write_results  # noqa: PLC0415

        write_results(self, output_dir)

    def plot(self) -> "Figure":
        """Overview figure: reference line and one line per candidate (spec 11.1)."""
        from homolocus.visualization.region_plot import overview_figure  # noqa: PLC0415

        return overview_figure(self)

    def plot_candidate(self, rank: int) -> "Figure":
        """Detailed figure for the candidate at a 1-based rank (spec 11.2)."""
        from homolocus.visualization.region_plot import candidate_figure  # noqa: PLC0415

        return candidate_figure(self, rank)


_EMPTY_MATCH: Row = {
    "reference_gene_id": None,
    "family_id": None,
    "family_name": None,
    "identity": float("nan"),
    "query_coverage": float("nan"),
    "subject_coverage": float("nan"),
    "bit_score": float("nan"),
    "is_best_copy": False,
}


def _summary_row(rank: int, cluster: Cluster) -> Row:
    start, end = to_one_based(cluster.start, cluster.end)
    display_start, display_end = to_one_based(cluster.display_start, cluster.display_end)
    parts = cluster.components
    return {
        "rank": rank,
        "candidate_id": RegionResults.candidate_id(rank),
        "contig": cluster.contig_id,
        "start": start,
        "end": end,
        "display_start": display_start,
        "display_end": display_end,
        "span_bp": cluster.span,
        "score": cluster.score,
        "content": parts.content,
        "similarity": parts.similarity,
        "order": parts.order,
        "orientation": parts.orientation,
        "spacing": parts.spacing,
        "reference_coverage": cluster.reference_coverage,
        "mean_protein_coverage": cluster.mean_protein_coverage,
        "n_target_genes": len(cluster.evidence_genes),
        "n_matched_genes": len(cluster.event_ids),
        "n_missing_families": len(cluster.missing_families),
        "flags": candidate_flags(cluster),
    }


def candidate_flags(cluster: Cluster) -> str:
    """Set candidate flags joined with ``;``: split, fusion, partial, edge (spec 14)."""
    return _flags(
        split=cluster.potential_split,
        fusion=cluster.potential_fusion,
        partial=cluster.contains_partial_genes,
        edge=cluster.touches_contig_edge,
    )


def _gene_columns(cluster: Cluster, gene: Gene) -> Row:
    start, end = to_one_based(gene.start, gene.end)
    fragments = {gene_id for split in cluster.splits for gene_id in split.fragment_ids}
    fused = {fusion.target_gene_id for fusion in cluster.fusions}
    return {
        "target_gene_id": gene.id,
        "locus_tag": gene.locus_tag,
        "product": gene.product,
        "start": start,
        "end": end,
        "strand": strand_char(gene.strand),
        "role": cluster.role(gene.id),
        "flags": _flags(
            partial=gene.partial,
            split_fragment=gene.id in fragments,
            fused=gene.id in fused,
        ),
    }


def _match_columns(match: GeneMatch, family: FamilyEvidence) -> Row:
    pairing = (family.best_member_id, family.best_copy_id)
    return {
        "reference_gene_id": match.reference_gene_id,
        "family_id": family.family_id,
        "family_name": family.family_name,
        "identity": match.identity,
        "query_coverage": match.query_coverage,
        "subject_coverage": match.subject_coverage,
        "bit_score": match.bit_score,
        "is_best_copy": family.represented
        and pairing == (match.reference_gene_id, match.target_gene_id),
    }


def _flags(**flags: bool) -> str:
    """Names of the set flags, in the given order, joined with ``;``."""
    return ";".join(name for name, is_set in flags.items() if is_set)
