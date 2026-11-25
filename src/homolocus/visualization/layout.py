"""What to draw: tracks of gene glyphs built from results (spec 11). Pure, no Matplotlib.

Families take the eight validated palette hues in reference order; family 9
reuses hue 1. Every coloured gene carries its family ID as a label, so colour
is never the only cue. Coordinates are 1-based inclusive.
"""

from collections.abc import Collection, Iterable, Mapping, Sequence
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
    opacity: float = 1.0


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


def identity_opacity(identity: float, min_identity: float) -> float:
    """Opacity of a matched gene: 0.5 at ``min_identity``, 1 at 100 % identity (spec 11)."""
    if min_identity >= 100:
        return 1.0
    share = (identity - min_identity) / (100 - min_identity)
    return 0.5 + 0.5 * min(max(share, 0.0), 1.0)


def arrow_labels(families: ReferenceFamilies, genes: Iterable[Gene]) -> dict[str, str]:
    """Family ID to its arrow label: the members' ``/gene`` names, else the ID (spec 11).

    Names are unique, in genomic order, joined with ``/``. Products never label arrows.
    """
    by_id = {gene.id: gene for gene in genes}
    labels = {}
    for family in families:
        names = (by_id[member].gene_name for member in family.gene_ids)
        unique = list(dict.fromkeys(name for name in names if name))
        labels[family.id] = "/".join(unique) or family.id
    return labels


def reference_track(results: RegionResults, *, missing: Collection[str] = ()) -> Track:
    """The reference region, coloured by family; ``missing`` families in mid grey."""
    reference = results.reference
    colours = family_colours(results.families)
    labels = arrow_labels(results.families, results.reference.genes)
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
    threshold = results.search.min_identity
    colours = family_colours(results.families)
    labels = arrow_labels(results.families, results.reference.genes)
    best_copies = {family.best_copy_id for family in cluster.represented_families}
    glyphs = []
    for gene in cluster.genes:
        if gene.id in cluster.event_ids:
            family_id, identity = best_family[gene.id]
            outlined = detail and gene.id in best_copies
            family = (family_id, labels[family_id])
            glyph = _glyph(gene, colours[family_id], "matched", family, outlined)
            glyphs.append(replace(glyph, opacity=identity_opacity(identity, threshold)))
        else:
            glyphs.append(_glyph(gene, NEUTRAL_GREY, "additional"))
    start, end = to_one_based(cluster.display_start, cluster.display_end)
    return Track(
        title=_candidate_title(results, rank, cluster),
        contig_id=cluster.contig_id,
        start=start,
        end=end,
        glyphs=tuple(glyphs),
        reversed=cluster.reversed,
    )


def overview_tracks(results: RegionResults, candidates: int = 0) -> tuple[Track, ...]:
    """The reference line followed by one line per drawn candidate (spec 11.1).

    At most ``candidates`` lines are drawn, best-ranked first (0 draws all);
    the reference note says how many more are in summary.tsv.
    """
    reference = reference_track(results)
    if not results.candidates:
        return (replace(reference, note="no candidates"),)
    drawn = _drawn(results, candidates)
    hidden = len(results.candidates) - drawn
    lines = (candidate_track(results, rank) for rank in range(1, drawn + 1))
    return (replace(reference, note=_more_note(hidden)), *lines)


def genome_set_tracks(results: GenomeSetResults, candidates: int = 0) -> tuple[Track, ...]:
    """The reference line, then the drawn candidates of every genome (spec 11).

    Genomes keep their input order and candidates their rank, at most
    ``candidates`` per genome (0 draws all). Each candidate title starts with
    its genome; the reference note names genomes without candidates and
    those with more candidates than drawn.
    """
    genomes = list(results.labelled_items())
    empty = [label for label, found in genomes if not found.candidates]
    hidden = [
        (label, len(found.candidates) - _drawn(found, candidates)) for label, found in genomes
    ]
    note = _joined(_no_candidates_note(empty), _hidden_note(hidden))
    reference = replace(reference_track(next(iter(results))), note=note)
    lines = (
        replace(track, title=f"{name} · {track.title}")
        for name, found in genomes
        for track in overview_tracks(found, candidates)[1:]
        if found.candidates
    )
    return (reference, *lines)


def paginate(tracks: Sequence[Track], lines_per_page: int) -> tuple[tuple[Track, ...], ...]:
    """Split overview tracks into pages of at most ``lines_per_page`` lines (spec 11.1).

    Every page starts with the reference line (the first track); with more
    than one page, its note ends with ``page i/n``.
    """
    reference, *lines = tracks
    size = lines_per_page - 1
    chunks = [lines[at : at + size] for at in range(0, len(lines), size)] or [[]]
    if len(chunks) == 1:
        return ((reference, *chunks[0]),)
    return tuple(
        (replace(reference, note=_joined(reference.note, f"page {number}/{len(chunks)}")), *chunk)
        for number, chunk in enumerate(chunks, start=1)
    )


def _drawn(results: RegionResults, candidates: int) -> int:
    total = len(results.candidates)
    return total if candidates == 0 else min(total, candidates)


def _more_note(hidden: int) -> str:
    if not hidden:
        return ""
    noun = "candidate" if hidden == 1 else "candidates"
    return f"{hidden} more {noun} in summary.tsv"


def _hidden_note(hidden: list[tuple[str, int]]) -> str:
    parts = [f"{label} +{count}" for label, count in hidden if count]
    return f"more candidates in summary.tsv: {', '.join(parts)}" if parts else ""


def _no_candidates_note(names: list[str]) -> str:
    return f"no candidates: {', '.join(names)}" if names else ""


def _joined(*notes: str) -> str:
    return " · ".join(note for note in notes if note)


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


def _best_families(results: RegionResults, cluster: Cluster) -> Mapping[str, tuple[str, float]]:
    """For each matched gene, the family and identity of its best-bit-score passing match.

    Ties go to the family that comes first in reference order (spec 7.6, 11).
    """
    search = results.search
    rank = {family.id: n for n, family in enumerate(results.families)}
    by_gene: dict[str, list[tuple[float, int, str, float]]] = {}
    for match in cluster.matches:
        if match.passes(max_evalue=search.evalue, min_identity=search.min_identity):
            family_id = results.families.family_of(match.reference_gene_id).id
            option = (match.bit_score, -rank[family_id], family_id, match.identity)
            by_gene.setdefault(match.target_gene_id, []).append(option)
    return {gene_id: max(options)[2:] for gene_id, options in by_gene.items()}


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
