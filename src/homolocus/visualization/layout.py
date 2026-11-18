"""What to draw: tracks of gene glyphs built from results (spec 11). Pure, no Matplotlib.

Families take the eight validated palette hues in reference order; family 9
reuses hue 1. Every coloured gene carries its family ID as a label, so colour
is never the only cue. Coordinates are 1-based inclusive.
"""

from collections.abc import Collection, Mapping
from dataclasses import dataclass, replace
from typing import Literal

from homolocus.coordinates import to_one_based
from homolocus.models.cluster import Cluster
from homolocus.models.gene import Gene, Strand
from homolocus.models.gene_family import ReferenceFamilies
from homolocus.results import GenomeSetResults, RegionResults, candidate_flags

PALETTE = (
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#6250d6",
    "#e34948",
)
"""Validated categorical hues (blue, orange, aqua, yellow, magenta, green, violet, red)."""

NEUTRAL_GREY = "#d9d8d2"
"""Additional target genes and reference genes without a protein."""

MISSING_GREY = "#9a998f"
"""Reference genes of families missing from a candidate (detail plot)."""

LABEL_NAME_LENGTH = 24
"""Longest family name shown in a plot label; longer names end in an ellipsis."""

GlyphRole = Literal["reference", "missing", "unsearchable", "matched", "additional"]


@dataclass(frozen=True)
class GeneGlyph:
    """One gene arrow: 1-based span, strand, colour, an optional label and its family."""

    gene_id: str
    start: int
    end: int
    strand: Strand
    colour: str
    label: str | None
    role: GlyphRole
    outlined: bool = False
    family_id: str | None = None


@dataclass(frozen=True)
class Track:
    """One line of a plot: a contig interval (1-based) and its gene glyphs.

    A ``reversed`` track is drawn on a descending axis, so it reads in the
    reference direction (spec 11.1).
    """

    title: str
    contig_id: str
    start: int
    end: int
    glyphs: tuple[GeneGlyph, ...]
    note: str = ""
    reversed: bool = False


def family_colours(families: ReferenceFamilies) -> dict[str, str]:
    """Family ID to hue, in family order; hues are reused past eight families."""
    return {family.id: PALETTE[n % len(PALETTE)] for n, family in enumerate(families)}


def family_label(family_id: str, name: str | None) -> str:
    """Plot label of a family: ``F001 sitA``, or the ID alone when unnamed (spec 11)."""
    if name is None:
        return family_id
    if len(name) > LABEL_NAME_LENGTH:
        name = name[: LABEL_NAME_LENGTH - 1] + "…"
    return f"{family_id} {name}"


def family_labels(families: ReferenceFamilies) -> dict[str, str]:
    """Family ID to plot label, in family order."""
    return {family.id: family_label(family.id, family.name) for family in families}


def reference_track(results: RegionResults, *, missing: Collection[str] = ()) -> Track:
    """The reference region, coloured by family; ``missing`` families in mid grey."""
    reference = results.reference
    colours, labels = family_colours(results.families), family_labels(results.families)
    glyphs = []
    for gene in reference.genes:
        if gene.protein_sequence is None:
            glyphs.append(_glyph(gene, NEUTRAL_GREY, "unsearchable"))
            continue
        family_id = results.families.family_of(gene.id).id
        if family_id in missing:
            glyphs.append(_glyph(gene, MISSING_GREY, "missing", (family_id, labels[family_id])))
        else:
            family = (family_id, labels[family_id])
            glyphs.append(_glyph(gene, colours[family_id], "reference", family))
    start, end = to_one_based(reference.start, reference.end)
    title = f"reference · {reference.contig_id}:{start}-{end}"
    return Track(title, reference.contig_id, start, end, tuple(glyphs))


def candidate_track(results: RegionResults, rank: int, *, detail: bool = False) -> Track:
    """A candidate's display interval; matched genes coloured by their best family."""
    cluster = results.candidates[rank - 1]
    best_family = _best_families(results, cluster)
    colours, labels = family_colours(results.families), family_labels(results.families)
    best_copies = {family.best_copy_id for family in cluster.represented_families}
    glyphs = []
    for gene in cluster.genes:
        if gene.id in cluster.event_ids:
            family_id = best_family[gene.id]
            outlined = detail and gene.id in best_copies
            family = (family_id, labels[family_id])
            glyphs.append(_glyph(gene, colours[family_id], "matched", family, outlined))
        else:
            glyphs.append(_glyph(gene, NEUTRAL_GREY, "additional"))
    start, end = to_one_based(cluster.display_start, cluster.display_end)
    return Track(
        title=_candidate_title(results, rank, cluster),
        contig_id=cluster.contig_id,
        start=start,
        end=end,
        glyphs=tuple(glyphs),
        note=_missing_note(cluster),
        reversed=cluster.reversed,
    )


def overview_tracks(results: RegionResults) -> tuple[Track, ...]:
    """The reference line followed by one line per candidate (spec 11.1)."""
    reference = reference_track(results)
    if not results.candidates:
        return (
            Track(
                reference.title,
                reference.contig_id,
                reference.start,
                reference.end,
                reference.glyphs,
                note="no candidates",
            ),
        )
    candidates = (candidate_track(results, rank) for rank in range(1, len(results.candidates) + 1))
    return (reference, *candidates)


def genome_set_tracks(results: GenomeSetResults) -> tuple[Track, ...]:
    """The reference line, then every candidate of every genome (spec 11).

    Genomes keep their input order and candidates their rank. Each candidate
    title starts with its genome; genomes without candidates are named in a
    note on the reference line.
    """
    empty = [name for name, found in results.items() if not found.candidates]
    reference = replace(reference_track(next(iter(results))), note=_no_candidates_note(empty))
    candidates = (
        replace(track, title=f"{name} · {track.title}")
        for name, found in results.items()
        for track in overview_tracks(found)[1:]
        if found.candidates
    )
    return (reference, *candidates)


def _no_candidates_note(names: list[str]) -> str:
    return f"no candidates: {', '.join(names)}" if names else ""


def detail_tracks(results: RegionResults, rank: int) -> tuple[Track, Track]:
    """Reference with missing families greyed, and the candidate with best copies outlined."""
    cluster = results.candidates[rank - 1]
    missing = {family.family_id for family in cluster.missing_families}
    return (
        reference_track(results, missing=missing),
        candidate_track(results, rank, detail=True),
    )


def _glyph(
    gene: Gene,
    colour: str,
    role: GlyphRole,
    family: tuple[str, str] | None = None,
    outlined: bool = False,
) -> GeneGlyph:
    """A glyph; ``family`` is the (ID, label) pair of a coloured gene."""
    family_id, label = family if family else (None, None)
    start, end = to_one_based(gene.start, gene.end)
    return GeneGlyph(gene.id, start, end, gene.strand, colour, label, role, outlined, family_id)


def _best_families(results: RegionResults, cluster: Cluster) -> Mapping[str, str]:
    """For each matched gene, the family of its best-bit-score passing match (spec 7.6).

    Ties go to the family that comes first in reference order.
    """
    search = results.search
    rank = {family.id: n for n, family in enumerate(results.families)}
    by_gene: dict[str, list[tuple[float, int, str]]] = {}
    for match in cluster.matches:
        if match.passes(max_evalue=search.evalue, min_identity=search.min_identity):
            family_id = results.families.family_of(match.reference_gene_id).id
            option = (match.bit_score, -rank[family_id], family_id)
            by_gene.setdefault(match.target_gene_id, []).append(option)
    return {gene_id: max(options)[2] for gene_id, options in by_gene.items()}


def _candidate_title(results: RegionResults, rank: int, cluster: Cluster) -> str:
    start, end = to_one_based(cluster.start, cluster.end)
    parts = [
        f"#{rank} {results.candidate_id(rank)}",
        f"{cluster.contig_id}:{start}-{end}",
        f"score {cluster.score:.3f}",
        f"coverage {cluster.reference_coverage:.2f}",
    ]
    flags = candidate_flags(cluster)
    marks = [flags] if flags else []
    return " · ".join([*parts, *marks, *(["reversed"] if cluster.reversed else [])])


def _missing_note(cluster: Cluster) -> str:
    missing = [family_label(f.family_id, f.family_name) for f in cluster.missing_families]
    return f"missing: {', '.join(missing)}" if missing else ""
