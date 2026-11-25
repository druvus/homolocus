"""Running DIAMOND searches (spec 6.1, 6.3, 6.5)."""

import stat
import time

import pytest

from fake_tools import fake_diamond
from homolocus import DiamondSearch, Gene
from homolocus.models.gene_match import aggregate
from homolocus.search.diamond import (
    DiamondError,
    DiamondNotFoundError,
    diamond_version,
    parse_tabular,
    search,
)
from homolocus.search.protein_index import ProteinIndex


def protein_gene(contig_id: str, number: int, protein: str | None) -> Gene:
    start = number * 10_000
    length = 3 * (len(protein) if protein else 10)
    return Gene(contig_id, start, start + length, 1, "A" * length, protein)


def genes(contig_id: str, *sequences: str | None) -> list[Gene]:
    return [protein_gene(contig_id, n, seq) for n, seq in enumerate(sequences)]


# Without DIAMOND


def test_missing_executable_gives_a_clear_error() -> None:
    with pytest.raises(DiamondNotFoundError, match="not found"):
        diamond_version("no-such-diamond-executable")


def test_version_is_read_from_the_executable(tmp_path) -> None:
    executable = fake_diamond(tmp_path, 'echo "diamond version 2.1.0"\n')
    assert diamond_version(executable) == "2.1.0"


def test_a_diamond_command_that_runs_over_its_timeout_is_stopped(tmp_path, proteins) -> None:
    executable = fake_diamond(
        tmp_path,
        'if [ "$1" = version ]; then echo "diamond version 2.1.0"; exit 0; fi\nexec sleep 30\n',
    )
    index = ProteinIndex.build(genes("r", proteins.random(50)), genes("t", proteins.random(50)))
    started = time.monotonic()
    with pytest.raises(DiamondError, match=r"DIAMOND makedb did not finish within 0\.5 s"):
        search(index, DiamondSearch(timeout=0.5), tmp_path / "raw.tsv", executable=executable)
    assert time.monotonic() - started < 10


def test_failing_diamond_raises_with_its_message(tmp_path, proteins) -> None:
    executable = fake_diamond(
        tmp_path,
        'if [ "$1" = version ]; then echo "diamond version 2.1.0"; exit 0; fi\n'
        'echo "database error" >&2; exit 3\n',
    )
    index = ProteinIndex.build(genes("r", proteins.random(50)), genes("t", proteins.random(50)))
    with pytest.raises(DiamondError, match=r"DIAMOND makedb failed \(exit 3\): database error"):
        search(index, DiamondSearch(), tmp_path / "raw.tsv", executable=executable)


def test_executable_that_cannot_start_raises_diamond_error(tmp_path) -> None:
    path = tmp_path / "broken-diamond"
    path.write_text("#!/no/such/interpreter\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    with pytest.raises(DiamondError, match="could not run DIAMOND version"):
        diamond_version(str(path))


@pytest.mark.parametrize(
    ("reference", "target"),
    [((None,), ("MKV",)), (("MKV",), (None,)), ((), ())],
)
def test_nothing_to_search_returns_no_alignments_and_an_empty_raw_file(
    tmp_path, reference, target
) -> None:
    index = ProteinIndex.build(genes("r", *reference), genes("t", *target))
    raw = tmp_path / "raw.tsv"
    assert search(index, DiamondSearch(), raw, executable="no-such-diamond") == ()
    assert raw.read_text(encoding="utf-8") == ""


def test_self_search_with_no_reference_proteins_needs_no_diamond(tmp_path) -> None:
    index = ProteinIndex.build(genes("r", None), genes("t", "MKV"))
    raw = tmp_path / "self.tsv"
    assert search(index, DiamondSearch(), raw, self_search=True, executable="none") == ()


# With DIAMOND


@pytest.mark.diamond
def test_homolog_is_found_with_canonical_ids(tmp_path, proteins) -> None:
    query = proteins.random(300)
    reference = genes("r", query)
    target = genes("t", proteins.random(250), proteins.mutate(query, 0.7), proteins.random(280))
    index = ProteinIndex.build(reference, target)
    alignments = search(index, DiamondSearch(), tmp_path / "raw.tsv")
    assert {(a.query_id, a.subject_id) for a in alignments} == {(reference[0].id, target[1].id)}
    assert 50 < alignments[0].identity < 90


@pytest.mark.diamond
def test_weak_secondary_hsp_between_the_cutoffs_is_kept(tmp_path, proteins) -> None:
    strong, weak = proteins.random(200), proteins.random(60)
    reference = genes("r", strong + proteins.random(40) + weak)
    target = genes(
        "t", proteins.mutate(weak, 0.4) + proteins.random(60) + proteins.mutate(strong, 0.8)
    )
    index = ProteinIndex.build(reference, target)
    loose = DiamondSearch(evalue=1e-30, hsp_evalue=1e-3)
    (pair,) = aggregate(search(index, loose, tmp_path / "loose.tsv"))
    assert len(pair.alignments) == 2
    assert max(hsp.evalue for hsp in pair.alignments) > loose.evalue
    assert pair.passes(max_evalue=loose.evalue, min_identity=30)
    strict = DiamondSearch(evalue=1e-30, hsp_evalue=1e-30)
    (strict_pair,) = aggregate(search(index, strict, tmp_path / "strict.tsv"))
    assert len(strict_pair.alignments) == 1
    assert pair.query_coverage > strict_pair.query_coverage


@pytest.mark.diamond
@pytest.mark.parametrize(("limit", "expected"), [(2, 2), (0, 5)])
def test_max_target_seqs_limits_subjects_per_query(tmp_path, proteins, limit, expected) -> None:
    query = proteins.random(300)
    target = genes("t", *(proteins.mutate(query, 0.8) for _ in range(5)))
    index = ProteinIndex.build(genes("r", query), target)
    alignments = search(index, DiamondSearch(max_target_seqs=limit), tmp_path / "raw.tsv")
    assert len({a.subject_id for a in alignments}) == expected


@pytest.mark.diamond
def test_self_search_excludes_self_hits_but_keeps_paralogs(tmp_path, proteins) -> None:
    family = proteins.random(300)
    reference = genes("r", family, proteins.mutate(family, 0.95), proteins.random(300))
    index = ProteinIndex.build(reference, [])
    alignments = search(index, DiamondSearch(), tmp_path / "self.tsv", self_search=True)
    assert {(a.query_id, a.subject_id) for a in alignments} == {
        (reference[0].id, reference[1].id),
        (reference[1].id, reference[0].id),
    }


@pytest.mark.diamond
def test_threads_setting_runs_and_gives_the_same_alignments(tmp_path, proteins) -> None:
    query = proteins.random(300)
    target = genes("t", proteins.mutate(query, 0.7), proteins.mutate(query, 0.6))
    index = ProteinIndex.build(genes("r", query), target)
    single = search(index, DiamondSearch(threads=1), tmp_path / "one.tsv")
    default = search(index, DiamondSearch(), tmp_path / "default.tsv")
    assert single == default
    assert single


@pytest.mark.diamond
def test_raw_output_is_diamonds_unchanged_token_table(tmp_path, proteins) -> None:
    query = proteins.random(300)
    index = ProteinIndex.build(genes("r", query), genes("t", proteins.mutate(query, 0.7)))
    raw = tmp_path / "raw.tsv"
    alignments = search(index, DiamondSearch(), raw)
    lines = raw.read_text(encoding="utf-8").splitlines()
    assert lines
    assert all(line.split("\t")[:2] == ["r000001", "t000001"] for line in lines)
    assert parse_tabular(lines, index) == alignments


@pytest.mark.diamond
def test_installed_diamond_is_supported() -> None:
    major = int(diamond_version().split(".")[0])
    assert major >= 2
