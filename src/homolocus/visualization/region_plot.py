"""Draw tracks with Matplotlib (spec 11). Rendering only: no searches, no detection.

Uses ``matplotlib.figure.Figure`` directly, so nothing touches pyplot's global
state and figures are never shown. Text uses ink colours, never series hues.
"""

from collections.abc import Sequence

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Polygon

from homolocus.results import RegionResults
from homolocus.visualization.layout import (
    GeneGlyph,
    Track,
    detail_tracks,
    family_colours,
    family_labels,
    overview_tracks,
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


def overview_figure(results: RegionResults) -> Figure:
    """Reference line and one line per candidate (spec 11.1)."""
    tracks = overview_tracks(results)
    figure = Figure(figsize=(11, 1.3 * len(tracks) + 0.9), layout="constrained")
    figure.set_facecolor(SURFACE)
    axes = figure.subplots(len(tracks), 1, squeeze=False)[:, 0]
    for axis, track in zip(axes, tracks, strict=True):
        _draw_track(axis, track)
    _family_legend(figure, results, tracks)
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
    _family_legend(figure, results, tracks)
    return figure


def _draw_track(axis: Axes, track: Track) -> None:
    span = track.end - track.start + 1
    for glyph in track.glyphs:
        axis.add_patch(
            Polygon(
                arrow_polygon(glyph, axis_span=span),
                closed=True,
                facecolor=glyph.colour,
                edgecolor=TEXT_PRIMARY if glyph.outlined else SURFACE,
                linewidth=1.8 if glyph.outlined else 0.6,
            )
        )
        if glyph.label:
            middle = (glyph.start + glyph.end) / 2
            axis.text(
                middle, 0.48, glyph.label, ha="center", va="bottom", fontsize=7, color=TEXT_PRIMARY
            )
    pad = _PAD_SHARE_OF_AXIS * span
    axis.set_xlim(track.start - pad, track.end + pad)
    axis.set_xticks([t for t in axis.get_xticks() if track.start <= t <= track.end])
    axis.set_ylim(-0.9, 1.0)
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
    _recessive_axes(axis, f"{track.contig_id} position (1-based)")


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


def _family_legend(figure: Figure, results: RegionResults, tracks: Sequence[Track]) -> None:
    shown = {g.family_id for track in tracks for g in track.glyphs if g.role != "missing"}
    colours, labels = family_colours(results.families), family_labels(results.families)
    handles = [
        Patch(facecolor=colour, edgecolor=SURFACE, label=labels[family_id])
        for family_id, colour in colours.items()
        if family_id in shown
    ]
    if handles:
        figure.legend(
            handles=handles,
            loc="outside lower center",
            ncols=min(len(handles), 10),
            frameon=False,
            fontsize=8,
            labelcolor=TEXT_PRIMARY,
        )
