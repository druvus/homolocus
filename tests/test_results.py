"""RegionResults: candidates, summary and gene rows, data frames (spec 12, 14)."""

from dataclasses import replace

import pandas as pd
import pytest

from homolocus import (
    Alignment,
    Contig,
    DiamondSearch,
    FamilyGrouping,
    Gene,
    GeneFamily,
    GeneMatch,
    Genome,
    ReferenceFamilies,
    ReferenceRegion,
    RegionDetection,
    RegionResults,
    ScoringWeights,
)
from homolocus.detection.evaluation import DetectionContext, detect
from homolocus.models.cluster import FusionEvent, SplitEvent
from homolocus.search.protein_index import ProteinIndex

PROTEIN = 300
REF_GENES = tuple(
    Gene("ref", n * 1200, n * 1200 + 900, 1, "A" * 900, "M" * PROTEIN) for n in range(4)
)
REFERENCE = ReferenceRegion("ref", 0, 4500, "A" * 4500, REF_GENES)
FAMILIES = ReferenceFamilies(
    tuple(GeneFamily(f"F{n + 1:03d}", (g.id,)) for n, g in enumerate(REF_GENES))
)


def target(start: int, strand: int = 1, **metadata) -> Gene:
    return Gene("t1", start, start + 900, strand, "A" * 900, "M" * 300, **metadata)


def match(ref: Gene, tgt: Gene, *, query=(0, PROTEIN), identity=80.0, bits=300.0) -> GeneMatch:
    span = query[1] - query[0]
    hsp = Alignment(ref.id, tgt.id, identity, span, *query, 0, span, PROTEIN, PROTEIN, 1e-40, bits)
    return GeneMatch(ref.id, tgt.id, (hsp,))


A = target(10_000, locus_tag="tA", product="alpha")
B = target(11_200)
INSERT = target(12_400, strand=-1, product="insertion")
C = target(13_600)
FAR_A, FAR_B = target(100_000), target(101_200)
MATCHES = (
    match(REF_GENES[0], A),
    match(REF_GENES[1], B),
    match(REF_GENES[3], B, query=(0, 60), identity=50.0, bits=40.0),
    match(REF_GENES[2], C),
    match(REF_GENES[0], FAR_A),
    match(REF_GENES[1], FAR_B),
)
GENOME = Genome(
    contigs={"t1": Contig("t1", "A" * 200_000, (A, B, INSERT, C, FAR_A, FAR_B))},
    name="target",
    source="target.fna",
)
SETTINGS = (DiamondSearch(), FamilyGrouping(), RegionDetection(min_matches=2), ScoringWeights())


def results(**overrides) -> RegionResults:
    search, grouping, detection, scoring = SETTINGS
    context = DetectionContext(REFERENCE, FAMILIES, search, detection, scoring)
    fields = {
        "reference": REFERENCE,
        "genome": GENOME,
        "families": FAMILIES,
        "matches": MATCHES,
        "candidates": detect(context, GENOME, MATCHES),
        "search": search,
        "family_grouping": grouping,
        "detection": detection,
        "scoring": scoring,
        "diamond_version": "2.2.9",
        "raw_main": "r000001\tt000001\n",
        "raw_self": "",
        "protein_index": ProteinIndex.build(REF_GENES, GENOME.genes()),
    }
    fields.update(overrides)
    return RegionResults(**fields)


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
