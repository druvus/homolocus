"""Rendering plots with Matplotlib (spec 11)."""

from dataclasses import replace

import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.patches import Polygon, Rectangle

from homolocus import DiamondSearch, Gene, GeneFamily, ReferenceFamilies, ReferenceRegion
from homolocus.visualization.layout import GeneGlyph, detail_tracks, overview_tracks
from homolocus.visualization.region_plot import (
    SURFACE,
    arrow_polygon,
    candidate_figure,
    legend_columns,
    overview_figure,
)
from sample_world import results

KEY = ["identity 30%", "identity 65%", "identity 100%"]  # the opacity key, default threshold
F003_LABEL = "F003 iron ABC transporter pe…"  # the sample F003 /product, cut at 24 characters


def glyph(start: int, end: int, strand: int) -> GeneGlyph:
    return GeneGlyph("g", start, end, strand, "#2a78d6", "F001", "matched")


def test_forward_arrow_points_right_with_a_capped_head() -> None:
    points = arrow_polygon(glyph(1001, 2000, 1), axis_span=10_000)
    xs = [x for x, _ in points]
    tip = max(points, key=lambda p: p[0])
    assert tip == (2000, 0.0)
    assert min(xs) == 1001
    assert sorted(set(xs))[-2] == pytest.approx(2000 - 150)  # head = 1.5 % of the axis


def test_reverse_arrow_points_left() -> None:
    points = arrow_polygon(glyph(1001, 2000, -1), axis_span=100_000)
    tip = min(points, key=lambda p: p[0])
    assert tip == (1001, 0.0)
    assert max(x for x, _ in points) == 2000


def test_short_gene_head_is_at_most_thirty_percent_of_its_length() -> None:
    points = arrow_polygon(glyph(1, 100, 1), axis_span=100_000)
    shoulder = sorted({x for x, _ in points})[-2]
    assert shoulder == pytest.approx(100 - 0.3 * 100)


def test_arrow_is_centred_on_zero_with_height_point_six() -> None:
    ys = [y for _, y in arrow_polygon(glyph(1, 100, 1), axis_span=1000)]
    assert ys == pytest.approx([-0.3, -0.3, 0.0, 0.3, 0.3])


def test_track_arrows_use_the_track_span_and_labels_sit_mid_gene() -> None:
    found = results()
    track = overview_tracks(found)[1]
    axis = overview_figure(found).axes[1]
    span = track.end - track.start + 1
    polygons = [p for p in axis.patches if isinstance(p, Polygon)]
    for patch, item in zip(polygons, track.glyphs, strict=True):
        assert patch.get_closed()
        expected = [v for point in arrow_polygon(item, axis_span=span) for v in point]
        drawn = patch.get_xy()[: len(expected) // 2].ravel().tolist()
        assert drawn == pytest.approx(expected)
    first = track.glyphs[0]
    label = next(t for t in axis.texts if t.get_text() == first.label)
    assert label.get_position()[0] == pytest.approx((first.start + first.end) / 2)


def test_axes_are_recessive_with_plain_ticks_and_a_frameless_legend() -> None:
    figure = candidate_figure(results(), 1)
    for axis in figure.axes:
        assert not axis.spines["top"].get_visible()
        assert not axis.spines["right"].get_visible()
    assert not figure.axes[0].spines["left"].get_visible()
    assert figure.axes[0].xaxis.get_major_formatter().get_useOffset() is False
    assert not figure.legends[0].get_frame_on()
    assert {axis.get_facecolor() for axis in figure.axes} == {to_rgba(SURFACE)}


def test_component_values_are_printed_just_right_of_their_bars() -> None:
    panel = candidate_figure(results(), 1).axes[2]
    bars = [p for p in panel.patches if isinstance(p, Rectangle)]
    for bar, text in zip(bars, panel.texts, strict=True):
        assert text.get_position()[0] == pytest.approx(bar.get_width() + 0.01)


def test_overview_has_one_axis_per_track_with_arrows_titles_and_notes() -> None:
    found = results()
    tracks = overview_tracks(found)
    figure = overview_figure(found)
    assert isinstance(figure, Figure)
    assert figure.canvas.manager is None
    assert len(figure.axes) == len(tracks) == 3
    for axis, track in zip(figure.axes, tracks, strict=True):
        assert axis.get_title(loc="left") == track.title
        assert len([p for p in axis.patches if isinstance(p, Polygon)]) == len(track.glyphs)
        pad = 0.01 * (track.end - track.start + 1)
        assert axis.get_xlim() == pytest.approx((track.start - pad, track.end + pad))
        assert all(track.start <= tick <= track.end for tick in axis.get_xticks())
    texts = [t.get_text() for t in figure.axes[1].texts]
    assert not [t for t in texts if t.startswith("missing")]
    assert "sitA" in texts


def test_overview_legend_lists_the_families() -> None:
    legend = overview_figure(results()).legends[0]
    assert [t.get_text() for t in legend.get_texts()] == [
        "F001 sitA",
        "F002",
        F003_LABEL,
        "F004",
        *KEY,
    ]


def test_detail_legend_leaves_out_families_drawn_grey() -> None:
    legend = candidate_figure(results(), 1).legends[0]
    assert [t.get_text() for t in legend.get_texts()] == [
        "F001 sitA",
        "F002",
        F003_LABEL,
        *KEY,
    ]


def test_best_copies_have_a_thicker_outline() -> None:
    figure = candidate_figure(results(), 1)
    reference, candidate = detail_tracks(results(), 1)
    widths = [p.get_linewidth() for p in figure.axes[1].patches if isinstance(p, Polygon)]
    outlined = [g.outlined for g in candidate.glyphs]
    assert all((w > 1) == flag for w, flag in zip(widths, outlined, strict=True))
    assert reference.glyphs[3].role == "missing"


def test_candidate_figure_has_a_component_panel() -> None:
    found = results()
    cluster = found.candidates[0]
    figure = candidate_figure(found, 1)
    assert len(figure.axes) == 3
    panel = figure.axes[2]
    bars = [p for p in panel.patches if isinstance(p, Rectangle)]
    parts = cluster.components
    expected = [parts.content, parts.similarity, parts.order, parts.orientation, parts.spacing]
    assert [bar.get_width() for bar in bars] == pytest.approx(expected)
    labels = [t.get_text() for t in panel.get_yticklabels()]
    assert labels == ["content", "similarity", "order", "orientation", "spacing"]
    assert panel.get_xlim() == (0.0, 1.0)
    assert f"{parts.order:.2f}" in [t.get_text() for t in panel.texts]
    assert f"score {cluster.score:.3f}" in figure.get_suptitle()


def test_results_return_figures_without_showing_them() -> None:
    found = results()
    assert len(found.plot().axes) == 3
    assert len(found.plot_candidate(2).axes) == 3
    assert found.plot().canvas.manager is None


def test_overview_without_candidates_shows_the_reference_and_a_note() -> None:
    figure = overview_figure(results(candidates=()))
    assert len(figure.axes) == 1
    assert "no candidates" in [t.get_text() for t in figure.axes[0].texts]


def test_reversed_lines_run_in_the_reference_direction_on_a_descending_axis() -> None:
    found = results()
    flipped = replace(found.candidates[1], reversed=True)
    mixed = replace(found, candidates=(found.candidates[0], flipped))
    reference, forward, backward = overview_figure(mixed).axes
    track = overview_tracks(mixed)[2]
    left, right = backward.get_xlim()
    assert left > right
    pad = 0.01 * (track.end - track.start + 1)
    assert (left, right) == pytest.approx((track.end + pad, track.start - pad))
    assert all(track.start <= tick <= track.end for tick in backward.get_xticks())
    assert backward.get_xlabel() == "t1 position (1-based, reversed)"
    for axis in (reference, forward):
        assert axis.get_xlim()[0] < axis.get_xlim()[1]
        assert "reversed" not in axis.get_xlabel()
    detail = candidate_figure(mixed, 2).axes[1]
    assert detail.get_xlim()[0] > detail.get_xlim()[1]


def label_texts(axis):
    return [t for t in axis.texts if t.get_text() and not t.get_text().startswith("no candidates")]


def test_neighbouring_labels_alternate_between_two_rows_in_drawn_order() -> None:
    found = results()
    flipped = replace(found.candidates[1], reversed=True)
    mixed = replace(found, candidates=(found.candidates[0], flipped))
    for axis in overview_figure(mixed).axes:
        drawn = sorted(label_texts(axis), key=lambda t: t.get_position()[0])
        if axis.xaxis_inverted():
            drawn.reverse()
        heights = [t.get_position()[1] for t in drawn]
        assert len(set(heights)) == min(2, len(heights))
        assert all(a != b for a, b in zip(heights, heights[1:], strict=False))
        assert heights[0] == min(heights)
        assert all(h < axis.get_ylim()[1] for h in heights)


def many_families():
    """A reference of 18 genes, each its own family with a long product name."""
    genes = tuple(
        Gene(
            "ref",
            n * 1000,
            n * 1000 + 900,
            1,
            "A" * 900,
            "M",
            product=f"type VI secretion protein {n}",
        )
        for n in range(18)
    )
    families = ReferenceFamilies(
        tuple(GeneFamily(f"F{n + 1:03d}", (g.id,), g.product) for n, g in enumerate(genes))
    )
    reference = ReferenceRegion("ref", 0, 18_000, "A" * 18_000, genes)
    return results(reference=reference, families=families, candidates=())


def test_a_long_legend_wraps_into_rows_and_fits_the_figure() -> None:
    figure = overview_figure(many_families())
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    (legend,) = figure.legends
    box, frame = legend.get_window_extent(), figure.bbox
    assert frame.x0 <= box.x0
    assert box.x1 <= frame.x1
    texts = legend.get_texts()
    assert len(texts) == 18
    assert len({round(t.get_window_extent().y0) for t in texts}) > 1  # several rows
    assert len({round(t.get_window_extent().x0) for t in texts}) > 1  # several columns


def test_legend_columns_are_bounded_by_the_entries_and_at_least_one() -> None:
    assert legend_columns(["a", "b"], 11) == 2
    assert legend_columns(["x" * 400], 2) == 1
    assert 1 < legend_columns(["F001 anmK"] * 20, 11) < 20


def test_matched_arrows_fade_by_opacity_but_outlines_stay_solid() -> None:
    found = results()
    axis = candidate_figure(found, 1).axes[1]
    _, track = detail_tracks(found, 1)
    polygons = [p for p in axis.patches if isinstance(p, Polygon)]
    for patch, item in zip(polygons, track.glyphs, strict=True):
        assert patch.get_facecolor()[3] == pytest.approx(item.opacity)
        assert patch.get_edgecolor()[3] == 1.0
    assert {round(g.opacity, 3) for g in track.glyphs} == {0.857, 1.0}


def test_the_opacity_key_spans_half_to_full_and_needs_matched_genes() -> None:
    legend = overview_figure(results()).legends[0]
    swatches = legend.legend_handles[-3:]
    assert [round(h.get_facecolor()[3], 3) for h in swatches] == [0.5, 0.75, 1.0]
    strict = overview_figure(results(search=DiamondSearch(min_identity=50.0))).legends[0]
    assert [t.get_text() for t in strict.get_texts()][-3:] == [
        "identity 50%",
        "identity 75%",
        "identity 100%",
    ]
    reference_only = overview_figure(results(candidates=())).legends[0]
    assert not [t for t in reference_only.get_texts() if "identity" in t.get_text()]
