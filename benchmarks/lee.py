"""Does homolocus find the LEE where it belongs, and nowhere else? (spec 17, 18)

The locus of enterocyte effacement is a well-studied pathogenicity island.
It is present in enterohaemorrhagic and enteropathogenic Escherichia coli
and absent from commensal K-12 and from uropathogenic CFT073, which makes
it a reference region with known answers.

The reference region is derived from the Sakai annotation rather than from
hard-coded coordinates: the densest cluster of CDS carrying LEE gene names
defines it, and the script reports what it used. Target genomes are read as
FASTA, so their genes are predicted by Pyrodigal and no result depends on
how each submission was annotated.

This is a report, not a test. It prints what it found; it asserts nothing
and changes no default. Run it with ``python benchmarks/lee.py``; it needs
network access and DIAMOND, and downloads about 45 MB into benchmarks/data.
"""

import argparse
import logging
import time
import urllib.request
from dataclasses import dataclass, replace
from importlib.metadata import version
from pathlib import Path

from homolocus import Genome, ReferenceRegion, RegionFinder
from homolocus.models.cluster import Cluster
from homolocus.results import GenomeSetResults, RegionResults, candidate_flags
from homolocus.search.diamond import diamond_version

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
USER_AGENT = "homolocus-benchmark"
"""Identifies the client to NCBI. It carries no personal information."""

REQUEST_PAUSE_SECONDS = 1.0
"""NCBI asks for no more than three requests a second without an API key."""

LEE_GENE_NAMES = frozenset(
    {"ler", "eae", "tir", "map", "grlA", "grlR", "etgA", "mpc", "escN", "escV"}
)
LEE_GENE_PREFIXES = ("esc", "esp", "ces", "sep")
"""The LEE's gene names: the secretion apparatus, its chaperones and effectors."""

CLUSTER_GAP_BP = 15_000
"""Named genes further apart than this belong to different loci."""


@dataclass(frozen=True)
class Target:
    """One target genome: its accession, a readable label and what we expect."""

    accession: str
    label: str
    carries_lee: bool


TARGETS = (
    Target("NC_002695.2", "O157:H7 Sakai", carries_lee=True),
    Target("AE005174.2", "O157:H7 EDL933", carries_lee=True),
    Target("FM180568.1", "O127:H6 E2348/69", carries_lee=True),
    Target("NC_000913.3", "K-12 MG1655", carries_lee=False),
    Target("AE014075.1", "CFT073", carries_lee=False),
)
REFERENCE_ACCESSION = "NC_002695.2"
"""Sakai provides the reference region, and is also searched as a control."""

logger = logging.getLogger("benchmarks.lee")


def download(accession: str, directory: Path, *, annotated: bool) -> Path:
    """Fetch one record from NCBI into ``directory``, or reuse the local copy."""
    suffix, rettype = ("gbk", "gbwithparts") if annotated else ("fna", "fasta")
    path = directory / f"{accession}.{suffix}"
    if path.exists():
        logger.info("have %s", path.name)
        return path
    url = f"{EFETCH}?db=nuccore&id={accession}&rettype={rettype}&retmode=text"
    logger.info("downloading %s", path.name)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=300) as response:
        path.write_bytes(response.read())
    time.sleep(REQUEST_PAUSE_SECONDS)
    return path


def lee_genes(genome: Genome) -> list[tuple[int, int, str]]:
    """Start, end and name of every CDS whose gene name belongs to the LEE."""
    return sorted(
        (gene.start, gene.end, gene.gene_name)
        for gene in genome.genes()
        if gene.gene_name and _is_lee_name(gene.gene_name)
    )


def _is_lee_name(name: str) -> bool:
    return name in LEE_GENE_NAMES or name.startswith(LEE_GENE_PREFIXES)


def densest_cluster(
    genes: list[tuple[int, int, str]], gap: int = CLUSTER_GAP_BP
) -> list[tuple[int, int, str]]:
    """The largest group of genes with no gap wider than ``gap`` between them.

    LEE gene names also occur on prophages elsewhere in the genome, so the
    island is the cluster holding most of them, not the outermost pair.
    """
    if not genes:
        raise ValueError("no LEE gene names in the reference annotation")
    clusters: list[list[tuple[int, int, str]]] = [[genes[0]]]
    for gene in genes[1:]:
        if gene[0] - clusters[-1][-1][1] > gap:
            clusters.append([gene])
        else:
            clusters[-1].append(gene)
    return max(clusters, key=len)


def reference_region(path: Path) -> tuple[ReferenceRegion, list[str]]:
    """The LEE of the annotated reference, and the gene names that defined it."""
    annotated = Genome.from_file(path)
    (contig_id,) = annotated.contigs
    cluster = densest_cluster(lee_genes(annotated))
    start, end = cluster[0][0] + 1, cluster[-1][1]  # 1-based inclusive
    logger.info(
        "LEE in %s: %s:%d-%d, %d named genes", path.name, contig_id, start, end, len(cluster)
    )
    region = ReferenceRegion.from_genbank(path, contig=contig_id, start=start, end=end)
    return region, [name for _, _, name in cluster]


def read_targets(directory: Path) -> list[Genome]:
    """Every target genome as FASTA, so Pyrodigal predicts its genes."""
    genomes = []
    for target in TARGETS:
        path = download(target.accession, directory, annotated=False)
        genomes.append(replace(Genome.from_file(path), display_name=target.label))
    return genomes


def report(results: GenomeSetResults, region: ReferenceRegion, names: list[str]) -> str:
    """The findings as Markdown: one section of context, then one row per genome."""
    lines = [
        "# LEE benchmark",
        "",
        "Generated by `benchmarks/lee.py`. It reports what homolocus found with its",
        "default settings; it asserts nothing. See section 18 of the specification.",
        "",
        f"- homolocus {version('homolocus')}, DIAMOND {next(iter(results)).diamond_version}",
        f"- Reference: {REFERENCE_ACCESSION} {region.contig_id}:{region.start + 1}-{region.end}",
        f"  ({region.end - region.start} bp, {len(region.genes)} CDS,"
        f" {len(results.families)} families)",
        f"- Named LEE genes defining it: {', '.join(names)}",
        "",
        "## Candidates",
        "",
        "| Genome | LEE expected | Candidates | Best score | Reference coverage"
        " | Missing families | Flags |",
        "|---|---|---|---|---|---|---|",
    ]
    lines += [_row(target, results[target.label]) for target in TARGETS]
    lines += [
        "",
        "## Every candidate, with its score components",
        "",
        "| Genome | Rank | Location | Score | content | similarity | order"
        " | orientation | spacing |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    lines += [
        _candidate_row(label, rank, candidate)
        for label, found in results.labelled_items()
        for rank, candidate in enumerate(found.candidates, start=1)
    ]
    lines += [
        "",
        "## Families not represented",
        "",
        "With, for each one, the best match its reference genes had anywhere in the",
        "genome. A family can be missing because nothing matched; because no match",
        "reached the coverage a matching event needs (spec 7.2, either side); or",
        "because no copy covered enough of the reference gene to count as represented",
        "(spec 7.7, the reference side only).",
        "",
    ]
    lines += [_missing_line(label, found) for label, found in results.labelled_items()]
    return "\n".join(lines) + "\n"


def _candidate_row(label: str, rank: int, candidate: Cluster) -> str:
    parts = candidate.components
    return (
        f"| {label} | #{rank} | {candidate.contig_id}:{candidate.start + 1}-{candidate.end} |"
        f" {candidate.score:.3f} | {parts.content:.3f} | {parts.similarity:.3f} |"
        f" {parts.order:.3f} | {parts.orientation:.3f} | {parts.spacing:.3f} |"
    )


def _missing_line(label: str, found: RegionResults) -> str:
    """Which families the genome's best candidate did not represent, and why."""
    best = found.best_candidate
    if best is None:
        return f"- **{label}**: no candidate."
    if not best.missing_families:
        return f"- **{label}**: none."
    parts = [
        f"{family.family_id} {family.family_name or '(unnamed)'}"
        f" ({_best_match(found, family.family_id)})"
        for family in best.missing_families
    ]
    return f"- **{label}**: {'; '.join(parts)}."


def _best_match(found: RegionResults, family_id: str) -> str:
    """The strongest match of a family's reference genes, as identity and coverage."""
    gene_ids = {
        gene_id
        for family in found.families
        if family.id == family_id
        for gene_id in family.gene_ids
    }
    matches = [match for match in found.matches if match.reference_gene_id in gene_ids]
    if not matches:
        return "no match"
    best = max(matches, key=lambda match: match.bit_score)
    return (
        f"best match {best.identity:.0f} % identity,"
        f" coverage {best.query_coverage:.2f} of the reference gene"
        f" and {best.subject_coverage:.2f} of the target gene"
    )


def _row(target: Target, found: RegionResults) -> str:
    best = found.best_candidate
    if best is None:
        return f"| {target.label} | {_yes(target.carries_lee)} | 0 | - | - | - | - |"
    return (
        f"| {target.label} | {_yes(target.carries_lee)} | {len(found.candidates)} |"
        f" {best.score:.3f} | {best.reference_coverage:.2f} |"
        f" {len(best.missing_families)} | {candidate_flags(best) or '-'} |"
    )


def _yes(carries: bool) -> str:
    return "yes" if carries else "no"


def main() -> int:
    """Download, search and write the report. Returns a process exit code."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    here = Path(__file__).parent
    parser.add_argument("--data", type=Path, default=here / "data", help="downloaded genomes")
    parser.add_argument("--out", type=Path, default=here / "results" / "lee.md", help="report")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args.data.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    logger.info("DIAMOND %s", diamond_version())
    region, names = reference_region(download(REFERENCE_ACCESSION, args.data, annotated=True))
    results = RegionFinder(region).find_all(read_targets(args.data))
    args.out.write_text(report(results, region, names), encoding="utf-8")
    logger.info("report written to %s", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
