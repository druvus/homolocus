"""RegionResults: candidates, summary and gene rows, data frames (spec 12, 14)."""

from dataclasses import replace

import pandas as pd
import pytest

from homolocus import Gene, Genome, RegionResults
from homolocus.models.cluster import FusionEvent, SplitEvent
from sample_world import FAR_A, FAR_B, GENOME, INSERT, REF_GENES, REFERENCE, A, B, C, results


def test_candidates_are_ranked_and_the_best_is_first() -> None:
    found = results()
    assert [c.start for c in found.candidates] == [10_000, 100_000]
    assert found.best_candidate is found.candidates[0]
    assert found.warnings == ()


def test_best_candidate_is_none_without_candidates() -> None:
    assert results(candidates=()).best_candidate is None


@pytest.mark.parametrize(
    ("rank", "expected"), [(1, "candidate_001"), (12, "candidate_012"), (1000, "candidate_1000")]
)
def test_candidate_ids_follow_the_rank(rank: int, expected: str) -> None:
    assert RegionResults.candidate_id(rank) == expected


def test_summary_has_one_row_per_candidate_with_one_based_coordinates() -> None:
    frame = results().to_dataframe()
    assert list(frame.columns) == [
        "rank",
        "candidate_id",
        "contig",
        "start",
        "end",
        "display_start",
        "display_end",
        "span_bp",
        "direction",
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
        "missing_families",
        "flags",
    ]
    first = frame.iloc[0]
    assert (first["rank"], first["candidate_id"], first["contig"]) == (1, "candidate_001", "t1")
    assert (first["start"], first["end"], first["span_bp"]) == (10_001, 14_500, 4_500)
    assert (first["display_start"], first["display_end"]) == (10_001, 14_500)
    assert (first["n_target_genes"], first["n_matched_genes"], first["n_missing_families"]) == (
        4,
        3,
        1,
    )
    assert first["reference_coverage"] == pytest.approx(0.75)
    best = results().best_candidate
    assert first["score"] == pytest.approx(best.score)
    assert first["order"] == best.components.order
    assert first["flags"] == ""
    assert first["direction"] == "forward"
    assert list(frame["missing_families"]) == ["F004", "F003;F004"]


def test_summary_tsv_lists_missing_families_and_blank_without_any(tmp_path) -> None:
    found = results()
    found.write(tmp_path)
    lines = (tmp_path / "summary.tsv").read_text().splitlines()
    column = lines[0].split("\t").index("missing_families")
    assert [line.split("\t")[column] for line in lines[1:]] == ["F004", "F003;F004"]
    complete = replace(found.candidates[0], families=found.candidates[0].represented_families)
    frame = replace(found, candidates=(complete,)).to_dataframe()
    assert list(frame["missing_families"]) == [""]


def test_summary_says_which_candidates_run_opposite_to_the_reference(tmp_path) -> None:
    found = results()
    flipped = replace(found.candidates[1], reversed=True)
    mixed = replace(found, candidates=(found.candidates[0], flipped))
    assert list(mixed.to_dataframe()["direction"]) == ["forward", "reversed"]
    mixed.write(tmp_path)
    lines = (tmp_path / "summary.tsv").read_text().splitlines()
    column = lines[0].split("\t").index("direction")
    assert [line.split("\t")[column] for line in lines[1:]] == ["forward", "reversed"]


def test_summary_flags_are_joined_in_a_fixed_order() -> None:
    found = results()
    flagged = replace(
        found.candidates[0],
        touches_contig_edge=True,
        splits=(SplitEvent(REF_GENES[0].id, (A.id, B.id), adjacent=True, same_strand=False),),
        fusions=(FusionEvent(B.id, ("F002", "F004")),),
    )
    frame = replace(found, candidates=(flagged, *found.candidates[1:])).to_dataframe()
    assert frame.iloc[0]["flags"] == "split;fusion;edge"


def test_empty_results_give_empty_frames_with_columns() -> None:
    empty = results(candidates=())
    assert empty.to_dataframe().empty
    assert "flags" in empty.to_dataframe().columns
    assert empty.genes_dataframe().empty
    assert "candidate_id" in empty.genes_dataframe().columns


def test_gene_rows_cover_every_gene_and_match() -> None:
    frame = results().genes_dataframe()
    first = frame[frame["candidate_id"] == "candidate_001"]
    assert list(first.columns) == [
        "candidate_id",
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
    ]
    assert list(first["target_gene_id"]) == [A.id, B.id, B.id, INSERT.id, C.id]
    assert list(first["reference_gene_id"].fillna("")) == [
        REF_GENES[0].id,
        REF_GENES[1].id,
        REF_GENES[3].id,
        "",
        REF_GENES[2].id,
    ]
    assert list(first["role"]) == ["matched", "matched", "matched", "additional", "matched"]
    assert list(first["family_id"].fillna("")) == ["F001", "F002", "F004", "", "F003"]
    assert list(first["family_name"].fillna("")) == [
        "sitA",
        "",
        "",
        "",
        "iron ABC transporter permease SitC",
    ]
    assert list(first["is_best_copy"]) == [True, True, False, False, True]


def test_gene_row_values_are_one_based_with_metadata() -> None:
    frame = results().genes_dataframe()
    row_a = frame.iloc[0]
    assert (row_a["start"], row_a["end"], row_a["strand"]) == (10_001, 10_900, "+")
    assert (row_a["locus_tag"], row_a["product"]) == ("tA", "alpha")
    assert row_a["identity"] == 80.0
    assert (row_a["query_coverage"], row_a["subject_coverage"], row_a["bit_score"]) == (
        1.0,
        1.0,
        300.0,
    )
    assert frame.iloc[1]["strand"] == "+"
    assert frame.iloc[3]["strand"] == "-"
    weak = frame.iloc[2]
    assert (weak["query_coverage"], weak["identity"]) == (0.2, 50.0)
    no_match = frame.iloc[3]
    assert pd.isna(no_match["identity"])
    assert not no_match["is_best_copy"]


def test_gene_row_flags() -> None:
    found = results()
    partial_insert = replace(INSERT, partial=True)
    genes = tuple(partial_insert if g.id == INSERT.id else g for g in found.candidates[0].genes)
    flagged = replace(
        found.candidates[0],
        genes=genes,
        splits=(
            SplitEvent(REF_GENES[0].id, (A.id, partial_insert.id), adjacent=True, same_strand=True),
        ),
        fusions=(FusionEvent(B.id, ("F002", "F004")),),
    )
    frame = replace(found, candidates=(flagged,)).genes_dataframe()
    flags = dict(zip(frame["target_gene_id"], frame["flags"], strict=False))
    assert flags[A.id] == "split_fragment"
    assert flags[B.id] == "fused"
    assert flags[INSERT.id] == "partial;split_fragment"
    assert flags[C.id] == ""


def test_gene_rows_of_every_candidate_carry_its_id() -> None:
    frame = results().genes_dataframe()
    assert list(frame["candidate_id"].unique()) == ["candidate_001", "candidate_002"]
    assert list(frame[frame["candidate_id"] == "candidate_002"]["target_gene_id"]) == [
        FAR_A.id,
        FAR_B.id,
    ]


# Writing the output directory (spec 14)

import json  # noqa: E402
import re  # noqa: E402

from homolocus import SkippedFeature, SkipReason  # noqa: E402


def tree(root) -> dict[str, bytes]:
    return {
        str(path.relative_to(root)): path.read_bytes()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def test_write_creates_the_whole_tree(tmp_path) -> None:
    results().write(tmp_path / "out")
    assert sorted(tree(tmp_path / "out")) == [
        "alignments/diamond.tsv",
        "alignments/protein_ids.tsv",
        "alignments/reference_self.tsv",
        "analysis.json",
        "candidates/candidate_001.tsv",
        "candidates/candidate_002.tsv",
        "plots/candidate_001.png",
        "plots/candidate_002.png",
        "plots/regions.png",
        "summary.tsv",
    ]


def test_plots_are_png_at_150_dpi(tmp_path) -> None:
    results().write(tmp_path)
    for name in ("regions.png", "candidate_001.png"):
        data = (tmp_path / "plots" / name).read_bytes()
        assert data.startswith(b"\x89PNG\r\n\x1a\n")
        assert int.from_bytes(data[16:20], "big") == 11 * 150  # width: 11 in at 150 dpi
        assert b"Software" not in data


def test_raw_diamond_tables_are_written_byte_for_byte(tmp_path) -> None:
    raw = "r000001\tt000001\t80.5\r\nodd line without newline"
    results(raw_main=raw, raw_self="self\n").write(tmp_path)
    assert (tmp_path / "alignments" / "diamond.tsv").read_bytes() == raw.encode()
    assert (tmp_path / "alignments" / "reference_self.tsv").read_text() == "self\n"
    mapping = (tmp_path / "alignments" / "protein_ids.tsv").read_text().splitlines()
    assert mapping[:2] == ["token\tgene_id", f"r000001\t{REF_GENES[0].id}"]


def test_tsv_files_match_the_data_frames(tmp_path) -> None:
    found = results()
    found.write(tmp_path)
    summary = pd.read_csv(tmp_path / "summary.tsv", sep="\t", keep_default_na=False)
    assert list(summary.columns) == list(found.to_dataframe().columns)
    assert list(summary["candidate_id"]) == ["candidate_001", "candidate_002"]
    assert list(summary["start"]) == [10_001, 100_001]
    genes = pd.read_csv(
        tmp_path / "candidates" / "candidate_001.tsv", sep="\t", keep_default_na=False
    )
    assert "candidate_id" not in genes.columns
    expected = found.genes_dataframe()
    expected = expected[expected["candidate_id"] == "candidate_001"]
    assert list(genes["target_gene_id"]) == list(expected["target_gene_id"])
    assert list(genes["reference_gene_id"]) == list(expected["reference_gene_id"].fillna(""))


def test_floats_have_six_significant_digits_and_missing_values_are_empty(tmp_path) -> None:
    results().write(tmp_path)
    lines = (tmp_path / "candidates" / "candidate_001.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    insert_row = dict(zip(header, lines[4].split("\t"), strict=True))
    assert insert_row["identity"] == ""
    assert insert_row["reference_gene_id"] == ""
    summary = (tmp_path / "summary.tsv").read_text().splitlines()
    score = dict(zip(summary[0].split("\t"), summary[1].split("\t"), strict=True))["score"]
    assert len(score.replace(".", "").lstrip("0")) <= 6


def test_candidate_tsv_has_family_names_and_blank_cells_without_one(tmp_path) -> None:
    results().write(tmp_path)
    lines = (tmp_path / "candidates" / "candidate_001.tsv").read_text().splitlines()
    rows = [dict(zip(lines[0].split("\t"), line.split("\t"), strict=True)) for line in lines[1:]]
    assert [row["family_name"] for row in rows] == [
        "sitA",
        "",
        "",
        "",
        "iron ABC transporter permease SitC",
    ]


def test_writing_twice_gives_identical_bytes(tmp_path) -> None:
    results().write(tmp_path / "a")
    results().write(tmp_path / "b")
    first = tree(tmp_path / "a")
    results().write(tmp_path / "a")
    assert tree(tmp_path / "a") == first == tree(tmp_path / "b")


def test_stale_candidate_files_are_removed_and_foreign_files_kept(tmp_path) -> None:
    (tmp_path / "candidates").mkdir()
    (tmp_path / "candidates" / "candidate_009.tsv").write_text("old")
    (tmp_path / "candidates" / "notes.txt").write_text("mine")
    (tmp_path / "plots").mkdir()
    (tmp_path / "plots" / "candidate_009.png").write_text("old")
    (tmp_path / "plots" / "notes.png").write_text("mine")
    (tmp_path / "README").write_text("mine too")
    results().write(tmp_path)
    assert not (tmp_path / "candidates" / "candidate_009.tsv").exists()
    assert (tmp_path / "candidates" / "notes.txt").read_text() == "mine"
    assert not (tmp_path / "plots" / "candidate_009.png").exists()
    assert (tmp_path / "plots" / "notes.png").read_text() == "mine"
    assert (tmp_path / "README").read_text() == "mine too"


def test_no_candidates_still_writes_headers(tmp_path) -> None:
    results(candidates=()).write(tmp_path)
    assert (tmp_path / "summary.tsv").read_text().startswith("rank\tcandidate_id\t")
    assert list((tmp_path / "candidates").iterdir()) == []
    assert [p.name for p in (tmp_path / "plots").iterdir()] == ["regions.png"]


def test_analysis_records_the_genome_display_name(tmp_path) -> None:
    shown = replace(GENOME, display_name="E. coli K-12 MG1655")
    results(genome=shown).write(tmp_path)
    data = json.loads((tmp_path / "analysis.json").read_text())
    assert data["target"]["display_name"] == "E. coli K-12 MG1655"


def analysis(tmp_path, **overrides) -> dict:
    results(**overrides).write(tmp_path)
    return json.loads((tmp_path / "analysis.json").read_text())


def test_analysis_records_versions_parameters_and_inputs(tmp_path) -> None:
    data = analysis(tmp_path)
    assert sorted(data) == [
        "families",
        "homolocus_version",
        "n_candidates",
        "parameters",
        "reference",
        "reference_genes_without_protein",
        "skipped_features",
        "target",
        "tool_versions",
        "untranslatable_genes",
        "warnings",
    ]
    assert data["tool_versions"]["diamond"] == "2.2.9"
    assert {"pyrodigal", "biopython"} <= set(data["tool_versions"])
    assert data["homolocus_version"]
    assert data["parameters"]["search"]["hsp_evalue"] == 1e-3
    assert data["parameters"]["detection"]["min_matches"] == 2
    assert data["parameters"]["family_grouping"]["family_min_identity"] == 90.0
    assert data["parameters"]["scoring"]["content"] == 0.4
    assert data["reference"] == {
        "contig": "ref",
        "start": 1,
        "end": 4500,
        "n_genes": 4,
        "source": None,
    }
    assert data["target"] == {
        "name": "target",
        "display_name": None,
        "source": "target.fna",
        "n_contigs": 1,
        "n_genes": 6,
    }
    names = ["sitA", None, "iron ABC transporter permease SitC", None]
    assert data["families"] == [
        {"id": f"F{n + 1:03d}", "name": names[n], "gene_ids": [g.id]}
        for n, g in enumerate(REF_GENES)
    ]
    assert data["n_candidates"] == 2
    assert data["warnings"] == []


def test_analysis_records_skipped_and_untranslatable_genes_one_based(tmp_path) -> None:
    no_protein = Gene("ref", 4000, 4300, 1, "A" * 300, None)
    reference = replace(
        REFERENCE,
        genes=(*REF_GENES, no_protein),
        skipped_features=(SkippedFeature("ref", 3700, 3790, "b9", SkipReason.PSEUDOGENE),),
    )
    lost_target = Gene("t1", 150_000, 150_300, -1, "A" * 300, None)
    genome = Genome(
        contigs={
            "t1": replace(GENOME.contigs["t1"], genes=(*GENOME.contigs["t1"].genes, lost_target))
        },
        skipped_features=(SkippedFeature("t1", 0, 30, None, SkipReason.COMPOUND_LOCATION),),
    )
    data = analysis(tmp_path, reference=reference, genome=genome, warnings=("careful",))
    assert data["reference_genes_without_protein"] == [no_protein.id]
    assert data["untranslatable_genes"] == {
        "reference": [no_protein.id],
        "target": [lost_target.id],
    }
    assert data["skipped_features"] == {
        "reference": [
            {"contig": "ref", "start": 3701, "end": 3790, "locus_tag": "b9", "reason": "pseudogene"}
        ],
        "target": [
            {
                "contig": "t1",
                "start": 1,
                "end": 30,
                "locus_tag": None,
                "reason": "compound_location",
            }
        ],
    }
    assert data["warnings"] == ["careful"]


def test_analysis_is_sorted_indented_and_has_no_timestamp(tmp_path) -> None:
    results().write(tmp_path)
    text = (tmp_path / "analysis.json").read_text()
    assert text.endswith("}\n")
    assert text.startswith('{\n  "families"')
    assert not re.search(r"\d{4}-\d{2}-\d{2}|\d{2}:\d{2}:\d{2}", text)  # no date or clock time
    assert not any(TIMESTAMP_KEY.fullmatch(key) for key in json_keys(json.loads(text)))


TIMESTAMP_KEY = re.compile(r"(timestamp|time|date|datetime|created|generated)(_at)?")


def json_keys(value) -> list[str]:
    """Every key at any depth of a parsed JSON value."""
    if isinstance(value, dict):
        return [key for key, inner in value.items() for key in (key, *json_keys(inner))]
    if isinstance(value, list):
        return [key for inner in value for key in json_keys(inner)]
    return []


def test_analysis_lists_parse_warnings_of_reference_and_target_before_run_warnings(
    tmp_path,
) -> None:
    reference = replace(REFERENCE, parse_warnings=("contig ref: pseudogene: 1",))
    genome = replace(
        GENOME, parse_warnings=("contig t1: 2 untranslatable", "contig t2: pseudogene: 1")
    )
    data = analysis(tmp_path, reference=reference, genome=genome, warnings=("careful",))
    assert data["warnings"] == [
        "reference: contig ref: pseudogene: 1",
        "target: contig t1: 2 untranslatable",
        "target: contig t2: pseudogene: 1",
        "careful",
    ]


from homolocus import PlotOptions  # noqa: E402


def many(n_candidates: int) -> RegionResults:
    """The sample results with their two candidates repeated to ``n_candidates``."""
    found = results()
    return replace(found, candidates=(found.candidates * n_candidates)[:n_candidates])


def titles(figure) -> list[str]:
    return [axis.get_title(loc="left") for axis in figure.axes if axis.get_title(loc="left")]


def test_plot_pages_split_a_long_overview_and_plot_returns_the_first() -> None:
    options = PlotOptions(candidates=0, lines_per_page=3)
    pages = many(5).plot_pages(options)
    assert [len(titles(page)) for page in pages] == [3, 3, 2]
    assert titles(many(5).plot(options)) == titles(pages[0])
    assert len(titles(many(5).plot())) == 4  # default: reference and three candidates


def test_write_adds_one_png_per_page_and_removes_stale_pages(tmp_path) -> None:
    many(5).write(tmp_path, plots=PlotOptions(candidates=0, lines_per_page=2))
    pages = sorted(p.name for p in (tmp_path / "plots").glob("regions*.png"))
    assert pages == ["regions.png", *(f"regions_00{n}.png" for n in range(2, 6))]
    many(5).write(tmp_path, plots=PlotOptions(candidates=0, lines_per_page=4))
    pages = sorted(p.name for p in (tmp_path / "plots").glob("regions*.png"))
    assert pages == ["regions.png", "regions_002.png"]


def test_writing_without_plots_removes_earlier_plots_only(tmp_path) -> None:
    many(3).write(tmp_path, plots=PlotOptions(candidates=0, lines_per_page=2))
    (tmp_path / "plots" / "notes.png").write_text("mine")
    many(3).write(tmp_path, plots=None)
    assert [p.name for p in (tmp_path / "plots").iterdir()] == ["notes.png"]
    assert (tmp_path / "candidates" / "candidate_003.tsv").is_file()
    results().write(tmp_path / "fresh", plots=None)
    assert not (tmp_path / "fresh" / "plots").exists()


def test_analysis_records_the_plot_options_or_none(tmp_path) -> None:
    results().write(tmp_path / "a", plots=PlotOptions(candidates=5))
    data = json.loads((tmp_path / "a" / "analysis.json").read_text())
    assert data["parameters"]["plots"] == {"candidates": 5, "lines_per_page": 30}
    results().write(tmp_path / "b", plots=None)
    assert json.loads((tmp_path / "b" / "analysis.json").read_text())["parameters"]["plots"] is None
