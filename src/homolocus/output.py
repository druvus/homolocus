"""Write a run's output directory (spec 14). The I/O layer for results.

Every file is deterministic: identical results give byte-identical files.
Managed files are summary.tsv, candidates/candidate_*.tsv, alignments/* and
analysis.json; stale candidate files are removed, other files are left alone.
"""

import json
from collections.abc import Iterable
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path
from typing import Any

import pandas as pd

from homolocus.coordinates import to_one_based
from homolocus.models.gene import Gene
from homolocus.models.skipped_feature import SkippedFeature
from homolocus.results import GENE_COLUMNS, RegionResults

FLOAT_FORMAT = "%.6g"
"""Floats in TSV files are written with 6 significant digits."""


def write_results(results: RegionResults, output_dir: str | Path) -> None:
    """Write summary, candidate tables, raw alignments and run metadata."""
    root = Path(output_dir)
    candidates, alignments = root / "candidates", root / "alignments"
    for directory in (candidates, alignments):
        directory.mkdir(parents=True, exist_ok=True)
    _write_tsv(results.to_dataframe(), root / "summary.tsv")
    _write_candidates(results, candidates)
    (alignments / "diamond.tsv").write_text(results.raw_main, encoding="utf-8", newline="")
    (alignments / "reference_self.tsv").write_text(results.raw_self, encoding="utf-8", newline="")
    results.protein_index.write_mapping(alignments / "protein_ids.tsv")
    metadata = json.dumps(analysis(results), indent=2, sort_keys=True)
    (root / "analysis.json").write_text(metadata + "\n", encoding="utf-8")


def _write_candidates(results: RegionResults, directory: Path) -> None:
    written = set()
    for rank, cluster in enumerate(results.candidates, start=1):
        path = directory / f"{results.candidate_id(rank)}.tsv"
        _write_tsv(pd.DataFrame(results.gene_rows(cluster), columns=list(GENE_COLUMNS)), path)
        written.add(path.name)
    for stale in directory.glob("candidate_*.tsv"):
        if stale.name not in written:
            stale.unlink()


def _write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n", float_format=FLOAT_FORMAT)


def analysis(results: RegionResults) -> dict[str, Any]:
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
            "n_contigs": len(genome.contigs),
            "n_genes": sum(1 for _ in genome.genes()),
        },
        "families": [{"id": f.id, "gene_ids": list(f.gene_ids)} for f in results.families],
        "reference_genes_without_protein": _without_protein(reference.genes),
        "untranslatable_genes": {
            "reference": _without_protein(reference.genes),
            "target": _without_protein(genome.genes()),
        },
        "skipped_features": {
            "reference": [_skipped(f) for f in reference.skipped_features],
            "target": [_skipped(f) for f in genome.skipped_features],
        },
        "warnings": list(results.warnings),
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
