"""Draw tracks with Matplotlib (spec 11). Rendering only: no searches, no detection.

Uses ``matplotlib.figure.Figure`` directly, so nothing touches pyplot's global
state and figures are never shown. Text uses ink colours, never series hues.
"""

from collections.abc import Sequence

from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch, Polygon
from matplotlib.textpath import TextToPath

from homolocus.config import DEFAULT_PLOTS, PlotOptions
from homolocus.models.gene_family import ReferenceFamilies
from homolocus.results import GenomeSetResults, RegionResults
from homolocus.visualization.layout import (
    GeneGlyph,
    Track,
    detail_tracks,
    family_colours,
    family_labels,
    genome_set_tracks,
    identity_opacity,
    overview_tracks,
    paginate,
)

TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
AXIS_INK = "#c3c2b7"
SURFACE = "#fcfcfb"
COMPONENTS = ("content", "similarity", "order", "orientation", "spacing")

_ARROW_HEIGHT = 0.6
_HEAD_SHARE_OF_GENE = 0.3
_HEAD_SHARE_OF_AXIS = 0.015
_PAD_SHARE_OF_AXIS = 0.01
_LABEL_ROWS = (0.42, 0.98)
_TOP = 1.6
_LEGEND_FONT = 8
_LEGEND_EXTRA_EMS = 2.0 + 0.8 + 2.4  # handle, handle-to-text pad, column spacing and margin
POINTS_PER_INCH = 72


def arrow_polygon(glyph: GeneGlyph, *, axis_span: int) -> list[tuple[float, float]]:
    """Corner points of a gene arrow pointing along its strand, centred on y = 0.

    The head is at most 30 % of the gene's length and 1.5 % of the axis span.
    """
    length = glyph.end - glyph.start + 1
    head = min(_HEAD_SHARE_OF_GENE * length, _HEAD_SHARE_OF_AXIS * axis_span)
    half = _ARROW_HEIGHT / 2
    if glyph.strand == 1:
        shoulder, base, tip = glyph.end - head, glyph.start, glyph.end
    else:
        shoulder, base, tip = glyph.start + head, glyph.end, glyph.start
    return [(base, -half), (shoulder, -half), (tip, 0.0), (shoulder, half), (base, half)]


def overview_figure(results: RegionResults, options: PlotOptions = DEFAULT_PLOTS) -> Figure:
    """First page of the overview: reference line and one line per drawn candidate (11.1)."""
    return overview_pages(results, options, first_only=True)[0]


def overview_pages(
    results: RegionResults, options: PlotOptions = DEFAULT_PLOTS, *, first_only: bool = False
) -> tuple[Figure, ...]:
    """Every page of the overview, each starting with the reference line (spec 11.1)."""
    tracks = overview_tracks(results, options.candidates)
    pages = paginate(tracks, options.lines_per_page)
    return _page_figures(pages, results.families, results.search.min_identity, first_only)


def genome_set_figure(results: GenomeSetResults, options: PlotOptions = DEFAULT_PLOTS) -> Figure:
    """First page of the combined overview of every genome (spec 11)."""
    return genome_set_pages(results, options, first_only=True)[0]


def genome_set_pages(
    results: GenomeSetResults, options: PlotOptions = DEFAULT_PLOTS, *, first_only: bool = False
) -> tuple[Figure, ...]:
    """Every page of the combined overview, each starting with the reference line (11)."""
    threshold = next(iter(results)).search.min_identity  # one finder, one search for all
    pages = paginate(genome_set_tracks(results, options.candidates), options.lines_per_page)
    return _page_figures(pages, results.families, threshold, first_only)


def _page_figures(
    pages: Sequence[Sequence[Track]],
    families: ReferenceFamilies,
    min_identity: float,
    first_only: bool,
) -> tuple[Figure, ...]:
    chosen = pages[:1] if first_only else pages
    return tuple(_lines_figure(page, families, min_identity) for page in chosen)


def _lines_figure(
    tracks: Sequence[Track], families: ReferenceFamilies, min_identity: float
) -> Figure:
    figure = Figure(figsize=(11, 1.3 * len(tracks) + 0.9), layout="constrained")
    figure.set_facecolor(SURFACE)
    axes = figure.subplots(len(tracks), 1, squeeze=False)[:, 0]
    for axis, track in zip(axes, tracks, strict=True):
        _draw_track(axis, track)
    _family_legend(figure, families, tracks, min_identity)
    return figure


def candidate_figure(results: RegionResults, rank: int) -> Figure:
    """Reference, candidate with every copy, and the score components (spec 11.2)."""
    cluster = results.candidates[rank - 1]
    tracks = detail_tracks(results, rank)
    figure = Figure(figsize=(11, 5.2), layout="constrained")
    figure.set_facecolor(SURFACE)
    axes = figure.subplots(3, 1, gridspec_kw={"height_ratios": [1, 1, 1.4]})
    for axis, track in zip(axes[:2], tracks, strict=True):
        _draw_track(axis, track)
    parts = cluster.components
    values = [parts.content, parts.similarity, parts.order, parts.orientation, parts.spacing]
    _draw_components(axes[2], values)
    figure.suptitle(
        f"{results.candidate_id(rank)} · score {cluster.score:.3f}",
        color=TEXT_PRIMARY,
        fontsize=11,
        x=0.01,
        ha="left",
    )
    _family_legend(figure, results.families, tracks, results.search.min_identity)
    return figure


def _draw_track(axis: Axes, track: Track) -> None:
    span = track.end - track.start + 1
    for glyph in track.glyphs:
        axis.add_patch(
            Polygon(
                arrow_polygon(glyph, axis_span=span),
                closed=True,
                facecolor=to_rgba(glyph.colour, glyph.opacity),
                edgecolor=TEXT_PRIMARY if glyph.outlined else SURFACE,
                linewidth=1.8 if glyph.outlined else 0.6,
            )
        )
    pad = _PAD_SHARE_OF_AXIS * span
    low, high = track.start - pad, track.end + pad
    axis.set_xlim((high, low) if track.reversed else (low, high))
    axis.set_xticks([t for t in axis.get_xticks() if track.start <= t <= track.end])
    _draw_labels(axis, track)
    axis.set_ylim(-0.9, _TOP)
    axis.set_title(track.title, loc="left", fontsize=9, color=TEXT_PRIMARY)
    if track.note:
        axis.text(
            1.0,
            1.02,
            track.note,
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color=TEXT_SECONDARY,
        )
    direction = ", reversed" if track.reversed else ""
    _recessive_axes(axis, f"{track.contig_id} position (1-based{direction})")


def _draw_labels(axis: Axes, track: Track) -> None:
    """Gene labels above the arrows, alternating between two rows in drawn order."""
    labels = sorted(
        ((g.start, (g.start + g.end) / 2, g.label) for g in track.glyphs if g.label is not None),
        reverse=track.reversed,
    )
    for number, (_, middle, label) in enumerate(labels):
        axis.text(
            middle,
            _LABEL_ROWS[number % 2],
            label,
            ha="center",
            va="bottom",
            fontsize=7,
            color=TEXT_PRIMARY,
        )


def _draw_components(axis: Axes, values: Sequence[float]) -> None:
    positions = list(range(len(COMPONENTS)))
    axis.barh(positions, values, height=0.55, color=TEXT_SECONDARY)
    for position, value in zip(positions, values, strict=True):
        axis.text(
            value + 0.01, position, f"{value:.2f}", va="center", fontsize=8, color=TEXT_PRIMARY
        )
    axis.set_yticks(positions, labels=list(COMPONENTS), fontsize=8, color=TEXT_PRIMARY)
    axis.invert_yaxis()
    axis.set_xlim(0.0, 1.0)
    axis.set_facecolor(SURFACE)
    axis.tick_params(axis="x", colors=TEXT_SECONDARY, labelsize=8)
    axis.set_title("score components", loc="left", fontsize=9, color=TEXT_PRIMARY)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_color(AXIS_INK)


def _recessive_axes(axis: Axes, label: str) -> None:
    axis.set_facecolor(SURFACE)
    axis.set_yticks([])
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color(AXIS_INK)
    axis.tick_params(axis="x", colors=TEXT_SECONDARY, labelsize=8)
    axis.ticklabel_format(axis="x", style="plain", useOffset=False)
    axis.set_xlabel(label, fontsize=8, color=TEXT_SECONDARY)


def _family_legend(
    figure: Figure, families: ReferenceFamilies, tracks: Sequence[Track], min_identity: float
) -> None:
    shown = {g.family_id for track in tracks for g in track.glyphs if g.role != "missing"}
    colours, labels = family_colours(families), family_labels(families)
    handles = [
        Patch(facecolor=colour, edgecolor=SURFACE, label=labels[family_id])
        for family_id, colour in colours.items()
        if family_id in shown
    ]
    if any(g.role == "matched" for track in tracks for g in track.glyphs):
        handles += _opacity_key(min_identity)
    if handles:
        figure.legend(
            handles=handles,
            loc="outside lower center",
            ncols=legend_columns([str(h.get_label()) for h in handles], figure.get_figwidth()),
            frameon=False,
            fontsize=_LEGEND_FONT,
            labelcolor=TEXT_PRIMARY,
        )


def _opacity_key(min_identity: float) -> list[Patch]:
    """Three grey swatches: identity at the threshold, midway, and 100 % (spec 11)."""
    identities = (min_identity, (min_identity + 100) / 2, 100.0)
    return [
        Patch(
            facecolor=to_rgba(TEXT_SECONDARY, identity_opacity(identity, min_identity)),
            edgecolor=SURFACE,
            label=f"identity {identity:g}%",
        )
        for identity in identities
    ]


def legend_columns(labels: Sequence[str], width_inches: float) -> int:
    """How many legend columns of these labels fit across the figure (spec 11)."""
    font = FontProperties(size=_LEGEND_FONT)
    measure = TextToPath().get_text_width_height_descent
    widest = max(measure(label, font, ismath=False)[0] for label in labels)
    entry = widest + _LEGEND_EXTRA_EMS * _LEGEND_FONT
    fitting = int(width_inches * POINTS_PER_INCH // entry)
    return max(1, min(len(labels), fitting))
