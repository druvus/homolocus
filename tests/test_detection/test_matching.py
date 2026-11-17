"""Matching events: target genes that may seed candidates (spec 7.2)."""

from homolocus import Alignment, Contig, Gene, GeneMatch
from homolocus.detection.matching import matching_events

REF_A, REF_B = "ref|1-300|+", "ref|401-700|+"
THRESHOLDS = {"max_evalue": 1e-5, "min_identity": 30.0, "min_event_coverage": 0.5}


def target_gene(start: int, contig_id: str = "t1", strand: int = 1) -> Gene:
    return Gene(contig_id, start, start + 300, strand, "A" * 300, "M" * 100)


def contig(*genes: Gene, contig_id: str = "t1") -> Contig:
    return Contig(contig_id, "A" * 20_000, genes)


def match(
    reference_id: str,
    target: Gene,
    *,
    query: tuple[int, int] = (0, 100),
    subject: tuple[int, int] = (0, 100),
    identity: float = 60.0,
    evalue: float = 1e-30,
) -> GeneMatch:
    length = max(query[1] - query[0], subject[1] - subject[0])
    hsp = Alignment(
        reference_id, target.id, identity, length, *query, *subject, 100, 100, evalue, 200.0
    )
    return GeneMatch(reference_id, target.id, (hsp,))


def events(contig_: Contig, *matches: GeneMatch) -> list[str]:
    return [gene.id for gene in matching_events(contig_, matches, **THRESHOLDS)]


def test_full_length_match_is_an_event() -> None:
    gene = target_gene(0)
    assert events(contig(gene), match(REF_A, gene)) == [gene.id]


def test_split_fragment_qualifies_by_subject_coverage() -> None:
    gene = target_gene(0)
    fragment = match(REF_A, gene, query=(0, 40), subject=(0, 100))
    assert events(contig(gene), fragment) == [gene.id]


def test_fused_protein_qualifies_by_query_coverage() -> None:
    gene = target_gene(0)
    fused = match(REF_A, gene, query=(0, 100), subject=(0, 40))
    assert events(contig(gene), fused) == [gene.id]


def test_short_shared_domain_is_not_an_event() -> None:
    gene = target_gene(0)
    assert events(contig(gene), match(REF_A, gene, query=(0, 49), subject=(10, 59))) == []


def test_event_coverage_threshold_is_inclusive_on_each_side() -> None:
    gene = target_gene(0)
    query_side = match(REF_A, gene, query=(0, 50), subject=(0, 49))
    subject_side = match(REF_A, gene, query=(0, 49), subject=(0, 50))
    assert events(contig(gene), query_side) == [gene.id]
    assert events(contig(gene), subject_side) == [gene.id]


def test_pair_level_thresholds_must_pass() -> None:
    gene = target_gene(0)
    assert events(contig(gene), match(REF_A, gene, evalue=1e-4)) == []
    assert events(contig(gene), match(REF_A, gene, identity=29.0)) == []


def test_any_qualifying_match_is_enough_and_a_gene_counts_once() -> None:
    gene = target_gene(0)
    weak = match(REF_A, gene, query=(0, 10), subject=(0, 10))
    strong = match(REF_B, gene)
    assert events(contig(gene), weak, strong) == [gene.id]
    assert events(contig(gene), strong, match(REF_A, gene)) == [gene.id]


def test_duplicated_copies_are_separate_events_in_genomic_order() -> None:
    copy_2, unrelated, copy_1 = target_gene(5000, strand=-1), target_gene(2500), target_gene(1000)
    found = events(contig(copy_2, unrelated, copy_1), match(REF_A, copy_2), match(REF_A, copy_1))
    assert found == [copy_1.id, copy_2.id]


def test_matches_to_other_contigs_are_ignored() -> None:
    here, elsewhere = target_gene(0), target_gene(0, contig_id="t2")
    assert events(contig(here), match(REF_A, elsewhere)) == []


def test_contig_without_matches_has_no_events() -> None:
    assert events(contig(target_gene(0))) == []
