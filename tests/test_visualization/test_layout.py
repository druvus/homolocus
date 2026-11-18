"""Plot layout: tracks and gene glyphs built from results (spec 11)."""

from dataclasses import replace

from homolocus import Gene, GeneFamily, ReferenceFamilies
from homolocus.models.cluster import SplitEvent
from homolocus.visualization.layout import (
    MISSING_GREY,
    NEUTRAL_GREY,
    PALETTE,
    candidate_track,
    detail_tracks,
    family_colours,
    family_label,
    overview_tracks,
    reference_track,
)
from sample_world import FAMILIES, INSERT, REF_GENES, REFERENCE, A, B, C, results

F003_LABEL = "F003 iron ABC transporter pe…"  # the sample F003 /product, cut at 24 characters

BLUE, ORANGE, AQUA, YELLOW = PALETTE[:4]


def test_palette_is_the_validated_eight_in_fixed_order() -> None:
    assert PALETTE == (
        "#2a78d6",
        "#eb6834",
        "#1baf7a",
        "#eda100",
        "#e87ba4",
        "#008300",
        "#6250d6",
        "#e34948",
    )


def test_families_take_hues_in_order_and_reuse_them_past_eight() -> None:
    nine = ReferenceFamilies(tuple(GeneFamily(f"F{n:03d}", (f"g{n}",)) for n in range(1, 10)))
    colours = family_colours(nine)
    assert [colours[f"F{n:03d}"] for n in range(1, 9)] == list(PALETTE)
    assert colours["F009"] == PALETTE[0]
    assert family_colours(FAMILIES) == {"F001": BLUE, "F002": ORANGE, "F003": AQUA, "F004": YELLOW}


def test_reference_track_shows_every_gene_one_based_and_coloured_by_family() -> None:
    track = reference_track(results())
    assert (track.contig_id, track.start, track.end) == ("ref", 1, 4500)
    assert "ref" in track.title
    assert [(g.start, g.end, g.strand) for g in track.glyphs] == [
        (1, 900, 1),
        (1201, 2100, 1),
        (2401, 3300, 1),
        (3601, 4500, 1),
    ]
    assert [g.colour for g in track.glyphs] == [BLUE, ORANGE, AQUA, YELLOW]
    assert [g.label for g in track.glyphs] == [
        "F001 sitA",
        "F002",
        F003_LABEL,
        "F004",
    ]
    assert [g.family_id for g in track.glyphs] == ["F001", "F002", "F003", "F004"]
    assert {g.role for g in track.glyphs} == {"reference"}


def test_paralogs_share_a_colour() -> None:
    paralogs = ReferenceFamilies(
        (GeneFamily("F001", (REF_GENES[0].id, REF_GENES[2].id)),)
        + tuple(GeneFamily(f"F00{n}", (REF_GENES[n - 1].id,)) for n in (2, 4))
    )
    track = reference_track(results(families=paralogs, candidates=()))
    assert track.glyphs[0].colour == track.glyphs[2].colour == BLUE
    assert track.glyphs[2].label == "F001"


def test_reference_gene_without_protein_is_grey_and_unlabelled() -> None:
    lost = Gene("ref", 4000, 4300, -1, "A" * 300, None)
    reference = replace(REFERENCE, genes=(*REF_GENES, lost))
    glyph = reference_track(results(reference=reference)).glyphs[-1]
    assert (glyph.colour, glyph.label, glyph.role, glyph.strand) == (
        NEUTRAL_GREY,
        None,
        "unsearchable",
        -1,
    )


def test_candidate_track_colours_matched_genes_by_their_best_family() -> None:
    track = candidate_track(results(), 1)
    assert (track.contig_id, track.start, track.end) == ("t1", 10_001, 14_500)
    by_id = {g.gene_id: g for g in track.glyphs}
    assert (by_id[A.id].colour, by_id[A.id].label, by_id[A.id].role) == (
        BLUE,
        "F001 sitA",
        "matched",
    )
    assert (by_id[B.id].colour, by_id[B.id].label) == (ORANGE, "F002")
    assert (by_id[C.id].colour, by_id[C.id].label) == (AQUA, F003_LABEL)


def test_additional_genes_are_grey_unlabelled_and_keep_their_strand() -> None:
    glyph = {g.gene_id: g for g in candidate_track(results(), 1).glyphs}[INSERT.id]
    assert (glyph.colour, glyph.label, glyph.role) == (NEUTRAL_GREY, None, "additional")
    assert (glyph.start, glyph.end, glyph.strand) == (12_401, 13_300, -1)


def test_candidate_title_and_note_carry_score_coverage_flags_and_missing_families() -> None:
    found = results()
    cluster = found.candidates[0]
    track = candidate_track(found, 1)
    assert "candidate_001" in track.title
    assert "t1:10001-14500" in track.title
    assert f"score {cluster.score:.3f}" in track.title
    assert "coverage 0.75" in track.title
    assert track.note == "missing: F004"
    flagged = replace(
        cluster,
        touches_contig_edge=True,
        splits=(SplitEvent(REF_GENES[0].id, (A.id, INSERT.id), adjacent=True, same_strand=False),),
    )
    title = candidate_track(replace(found, candidates=(flagged,)), 1).title
    assert title.endswith("split;edge")


def test_candidate_without_missing_families_has_no_note() -> None:
    found = results()
    second = candidate_track(found, 2)
    assert second.note == "missing: F003 iron ABC transporter pe…, F004"
    complete = replace(found.candidates[0], families=found.candidates[0].represented_families)
    assert candidate_track(replace(found, candidates=(complete,)), 1).note == ""


def test_overview_is_reference_then_each_candidate() -> None:
    tracks = overview_tracks(results())
    assert [t.contig_id for t in tracks] == ["ref", "t1", "t1"]
    assert [t.start for t in tracks[1:]] == [10_001, 100_001]
    assert not any(g.outlined for t in tracks for g in t.glyphs)


def test_overview_without_candidates_says_so() -> None:
    (only,) = overview_tracks(results(candidates=()))
    assert only.note == "no candidates"


def test_detail_greys_missing_families_and_outlines_best_copies() -> None:
    reference, candidate = detail_tracks(results(), 1)
    assert [g.role for g in reference.glyphs] == ["reference"] * 3 + ["missing"]
    assert (reference.glyphs[3].colour, reference.glyphs[3].label) == (MISSING_GREY, "F004")
    outlined = {g.gene_id for g in candidate.glyphs if g.outlined}
    assert outlined == {A.id, B.id, C.id}
    assert candidate.note == "missing: F004"


def test_failing_matches_never_decide_a_colour() -> None:
    from homolocus import GeneMatch  # noqa: PLC0415
    from sample_world import MATCHES, match  # noqa: PLC0415

    strong_failing = match(REF_GENES[3], A, identity=10.0, bits=900.0)
    found = results(matches=(*MATCHES, strong_failing))
    glyph = {g.gene_id: g for g in candidate_track(found, 1).glyphs}[A.id]
    assert isinstance(strong_failing, GeneMatch)
    assert (glyph.colour, glyph.family_id) == (BLUE, "F001")


def test_tied_bit_scores_pick_the_earlier_family() -> None:
    from sample_world import MATCHES, match  # noqa: PLC0415

    tie = match(REF_GENES[3], A, bits=300.0)
    glyph = {g.gene_id: g for g in candidate_track(results(matches=(tie, *MATCHES)), 1).glyphs}[
        A.id
    ]
    assert glyph.family_id == "F001"


def test_family_labels_add_the_name_and_cut_long_ones() -> None:
    assert family_label("F001", None) == "F001"
    assert family_label("F001", "sitA") == "F001 sitA"
    assert family_label("F002", "x" * 24) == "F002 " + "x" * 24
    assert family_label("F002", "x" * 25) == "F002 " + "x" * 23 + "…"
