"""Write a run's output directory (spec 14). The I/O layer for results.

Every file is deterministic: identical results give byte-identical files.
Managed files are summary.tsv, candidates/candidate_*.tsv, alignments/*,
analysis.json and plots/ (regions.png, candidate_*.png); stale candidate files
are removed, other files are left alone. Several genomes get one such
directory each, plus a combined summary.tsv and plots/regions.png.
"""

import json
import logging
from collections.abc import Iterable, Sequence
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from homolocus.config import DEFAULT_PLOTS, PlotOptions
from homolocus.coordinates import to_one_based
from homolocus.models.gene import Gene
from homolocus.models.skipped_feature import SkippedFeature
from homolocus.results import GENE_COLUMNS, GenomeSetResults, RegionResults
from homolocus.visualization.region_plot import (
    candidate_figure,
    genome_set_pages,
    overview_pages,
)

if TYPE_CHECKING:
    from matplotlib.figure import Figure

FLOAT_FORMAT = "%.6g"
"""Floats in TSV files are written with 6 significant digits."""

PLOT_DPI = 150
"""Resolution of the PNG plots (spec 11.4)."""

logger = logging.getLogger(__name__)


def write_results(
    results: RegionResults, output_dir: str | Path, plots: PlotOptions | None = DEFAULT_PLOTS
) -> None:
    """Write summary, candidate tables, raw alignments, run metadata and plots.

    ``plots=None`` writes no plots and removes those of an earlier run.
    """
    root = Path(output_dir)
    logger.info("writing results to %s", root)
    candidates, alignments = root / "candidates", root / "alignments"
    for directory in (candidates, alignments):
        directory.mkdir(parents=True, exist_ok=True)
    _write_tsv(results.to_dataframe(), root / "summary.tsv")
    _write_candidates(results, candidates)
    (alignments / "diamond.tsv").write_text(results.raw_main, encoding="utf-8", newline="")
    (alignments / "reference_self.tsv").write_text(results.raw_self, encoding="utf-8", newline="")
    results.protein_index.write_mapping(alignments / "protein_ids.tsv")
    metadata = json.dumps(analysis(results, plots), indent=2, sort_keys=True)
    (root / "analysis.json").write_text(metadata + "\n", encoding="utf-8")
    if plots is None:
        _remove_plots(root / "plots")
    else:
        _write_plots(results, root / "plots", plots)


def write_genome_set(
    results: GenomeSetResults, output_dir: str | Path, plots: PlotOptions | None = DEFAULT_PLOTS
) -> None:
    """Write each genome to ``output_dir/<genome>/``, then the combined files (spec 14)."""
    root = Path(output_dir)
    logger.info("writing results to %s", root)
    for name, genome_results in results.items():
        write_results(genome_results, root / name, plots)
    _write_tsv(results.to_dataframe(), root / "summary.tsv")  # root exists: a set has genomes
    if plots is None:
        _remove_plots(root / "plots")
    else:
        _write_pages(genome_set_pages(results, plots), root / "plots")


def _write_candidates(results: RegionResults, directory: Path) -> None:
    written = set()
    for rank, cluster in enumerate(results.candidates, start=1):
        path = directory / f"{results.candidate_id(rank)}.tsv"
        _write_tsv(pd.DataFrame(results.gene_rows(cluster), columns=list(GENE_COLUMNS)), path)
        written.add(path.name)
    _remove_stale(directory, "candidate_*.tsv", written)


def _write_plots(results: RegionResults, directory: Path, plots: PlotOptions) -> None:
    _write_pages(overview_pages(results, plots), directory)
    written = set()
    for rank in range(1, len(results.candidates) + 1):
        path = directory / f"{results.candidate_id(rank)}.png"
        _save_png(candidate_figure(results, rank), path)
        written.add(path.name)
    _remove_stale(directory, "candidate_*.png", written)


def _write_pages(pages: Sequence["Figure"], directory: Path) -> None:
    """Page 1 as ``regions.png``, later pages as ``regions_002.png``, ... (spec 11.4)."""
    directory.mkdir(parents=True, exist_ok=True)
    written = set()
    for number, page in enumerate(pages, start=1):
        name = "regions.png" if number == 1 else f"regions_{number:03d}.png"
        _save_png(page, directory / name)
        written.add(name)
    _remove_stale(directory, "regions_*.png", written)


def _remove_plots(directory: Path) -> None:
    """Remove the plot files an earlier run wrote; other files stay (spec 11.4, 14)."""
    for pattern in ("regions.png", "regions_*.png", "candidate_*.png"):
        _remove_stale(directory, pattern, set())


def _save_png(figure: "Figure", path: Path) -> None:
    """PNG without the Software stamp, so identical results give identical bytes."""
    figure.savefig(path, format="png", dpi=PLOT_DPI, metadata={"Software": None})


def _remove_stale(directory: Path, pattern: str, written: set[str]) -> None:
    for stale in directory.glob(pattern):
        if stale.name not in written:
            stale.unlink()


def _write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n", float_format=FLOAT_FORMAT)


def analysis(results: RegionResults, plots: PlotOptions | None = DEFAULT_PLOTS) -> dict[str, Any]:
    """Run metadata for ``analysis.json``: no timestamp, coordinates 1-based (spec 14)."""
    reference, genome = results.reference, results.genome
    reference_start, reference_end = to_one_based(reference.start, reference.end)
    return {
        "homolocus_version": version("homolocus"),
        "tool_versions": {
            "diamond": results.diamond_version,
            "pyrodigal": version("pyrodigal"),
            "biopython": version("biopython"),
        },
        "parameters": {
            "search": asdict(results.search),
            "family_grouping": asdict(results.family_grouping),
            "detection": asdict(results.detection),
            "scoring": asdict(results.scoring),
            "plots": None if plots is None else asdict(plots),
        },
        "reference": {
            "source": reference.source_genome,
            "contig": reference.contig_id,
            "start": reference_start,
            "end": reference_end,
            "n_genes": len(reference.genes),
        },
        "target": {
            "source": genome.source,
            "name": genome.name,
            "display_name": genome.display_name,
            "n_contigs": len(genome.contigs),
            "n_genes": sum(1 for _ in genome.genes()),
        },
        "families": [
            {"id": f.id, "name": f.name, "gene_ids": list(f.gene_ids)} for f in results.families
        ],
        "reference_genes_without_protein": _without_protein(reference.genes),
        "untranslatable_genes": {
            "reference": _without_protein(reference.genes),
            "target": _without_protein(genome.genes()),
        },
        "skipped_features": {
            "reference": [_skipped(f) for f in reference.skipped_features],
            "target": [_skipped(f) for f in genome.skipped_features],
        },
        "warnings": [
            *(f"reference: {warning}" for warning in reference.parse_warnings),
            *(f"target: {warning}" for warning in genome.parse_warnings),
            *results.warnings,
        ],
        "n_candidates": len(results.candidates),
    }


def _without_protein(genes: Iterable[Gene]) -> list[str]:
    return [gene.id for gene in genes if gene.protein_sequence is None]


def _skipped(feature: SkippedFeature) -> dict[str, Any]:
    start, end = to_one_based(feature.start, feature.end)
    return {
        "contig": feature.contig_id,
        "start": start,
        "end": end,
        "locus_tag": feature.locus_tag,
        "reason": feature.reason.value,
    }
