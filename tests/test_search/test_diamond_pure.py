"""DIAMOND command lines, tabular parsing and version checks, without DIAMOND (spec 6)."""

from pathlib import Path

import pytest

from homolocus import DiamondSearch, Gene
from homolocus.search.diamond import (
    OUTPUT_FIELDS,
    DiamondError,
    blastp_command,
    makedb_command,
    parse_tabular,
    parse_version,
)
from homolocus.search.protein_index import ProteinIndex


def gene(start: int, contig_id: str) -> Gene:
    return Gene(contig_id, start, start + 9, 1, "ATGAAAGTA", "MKV")


INDEX = ProteinIndex.build([gene(0, "r"), gene(20, "r")], [gene(0, "t"), gene(20, "t")])
R1, R2, T1, T2 = "r|1-9|+", "r|21-29|+", "t|1-9|+", "t|21-29|+"


def row(tokens: tuple[str, str], query: tuple[int, int], subject: tuple[int, int], **kw) -> str:
    qstart, qend = query
    sstart, send = subject
    values = {
        "pident": "55.5",
        "length": str(max(qend - qstart, send - sstart) + 1),
        "evalue": "1.2e-30",
        "bitscore": "110.5",
        "qlen": "300",
        "slen": "320",
    }
    values.update(kw)
    return "\t".join(
        [*tokens, values["pident"], values["length"], "3", "1", str(qstart), str(qend)]
        + [str(sstart), str(send), values["evalue"], values["bitscore"], values["qlen"]]
        + [values["slen"]]
    )


# Commands


def test_output_fields_follow_the_specification() -> None:
    assert OUTPUT_FIELDS == (
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
    )


def test_blastp_command_passes_every_setting_explicitly() -> None:
    command = blastp_command("diamond", Path("q.faa"), Path("db"), Path("out.tsv"), DiamondSearch())
    assert command[:8] == [
        "diamond",
        "blastp",
        "--query",
        "q.faa",
        "--db",
        "db",
        "--out",
        "out.tsv",
    ]
    assert command[8:10] == ["--outfmt", "6"]
    assert command[10:24] == list(OUTPUT_FIELDS)
    assert command[24:] == [
        "--evalue",
        "0.001",
        "--max-target-seqs",
        "100",
        "--max-hsps",
        "0",
        "--sensitive",
        "--quiet",
    ]


def test_blastp_never_filters_hsps_by_identity_or_coverage() -> None:
    command = blastp_command("diamond", Path("q"), Path("d"), Path("o"), DiamondSearch())
    assert not {"--id", "--query-cover", "--subject-cover", "--min-score"} & set(command)


def test_hsp_evalue_not_pair_evalue_goes_to_diamond() -> None:
    config = DiamondSearch(evalue=1e-10, hsp_evalue=1e-2)
    command = blastp_command("diamond", Path("q"), Path("d"), Path("o"), config)
    assert command[command.index("--evalue") + 1] == "0.01"


def test_default_sensitivity_adds_no_flag() -> None:
    command = blastp_command(
        "d", Path("q"), Path("d"), Path("o"), DiamondSearch(sensitivity="default")
    )
    assert not [arg for arg in command if arg.endswith("sensitive") or arg == "--fast"]


@pytest.mark.parametrize("sensitivity", ["fast", "mid-sensitive", "ultra-sensitive"])
def test_sensitivity_flag(sensitivity) -> None:
    command = blastp_command(
        "d", Path("q"), Path("d"), Path("o"), DiamondSearch(sensitivity=sensitivity)
    )
    assert f"--{sensitivity}" in command


def test_threads_go_to_both_makedb_and_blastp() -> None:
    config = DiamondSearch(threads=4)
    makedb = makedb_command("diamond", Path("t.faa"), Path("db"), config)
    blastp = blastp_command("diamond", Path("q"), Path("db"), Path("o"), config)
    assert makedb == [
        "diamond",
        "makedb",
        "--in",
        "t.faa",
        "--db",
        "db",
        "--quiet",
        "--threads",
        "4",
    ]
    assert blastp[-2:] == ["--threads", "4"]


def test_no_threads_flag_by_default() -> None:
    makedb = makedb_command("diamond", Path("t.faa"), Path("db"), DiamondSearch())
    blastp = blastp_command("diamond", Path("q"), Path("db"), Path("o"), DiamondSearch())
    assert "--threads" not in makedb + blastp


def test_counts_are_passed_as_given() -> None:
    config = DiamondSearch(max_target_seqs=0, max_hsps=3)
    command = blastp_command("d", Path("q"), Path("d"), Path("o"), config)
    assert command[command.index("--max-target-seqs") + 1] == "0"
    assert command[command.index("--max-hsps") + 1] == "3"


# Parsing


def test_row_becomes_an_alignment_with_internal_coordinates() -> None:
    (hsp,) = parse_tabular([row(("r000001", "t000002"), (11, 60), (21, 75))], INDEX)
    assert (hsp.query_id, hsp.subject_id) == (R1, T2)
    assert (hsp.query_start, hsp.query_end, hsp.subject_start, hsp.subject_end) == (10, 60, 20, 75)
    assert (hsp.identity, hsp.alignment_length) == (55.5, 55)
    assert (hsp.evalue, hsp.bit_score) == (1.2e-30, 110.5)
    assert (hsp.query_length, hsp.subject_length) == (300, 320)


def test_several_hsps_and_hits_are_all_kept_in_deterministic_order() -> None:
    lines = [
        row(("r000002", "t000001"), (1, 50), (1, 50)),
        row(("r000001", "t000002"), (101, 150), (101, 150)),
        row(("r000001", "t000002"), (1, 50), (1, 50)),
        row(("r000001", "t000001"), (1, 50), (1, 50)),
        "",
    ]
    parsed = parse_tabular(lines, INDEX)
    assert [(h.query_id, h.subject_id, h.query_start) for h in parsed] == [
        (R1, T1, 0),
        (R1, T2, 0),
        (R1, T2, 100),
        (R2, T1, 0),
    ]
    assert parse_tabular(reversed(lines), INDEX) == parsed


def test_query_and_subject_may_both_be_reference_tokens() -> None:
    (hsp,) = parse_tabular([row(("r000001", "r000002"), (1, 50), (1, 50))], INDEX)
    assert (hsp.query_id, hsp.subject_id) == (R1, R2)


@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("r000001\tt000001\t55.5", "line 1: expected 14 fields, found 3"),
        (row(("r000009", "t000001"), (1, 50), (1, 50)), "line 1: unknown token"),
        (row(("r000001", "t000001"), (1, 50), (1, 50), pident="high"), "line 1"),
        (row(("r000001", "t000001"), (1, 50), (1, 50), qlen="40"), "line 1"),
        (row(("r000001", "t000001"), (0, 50), (1, 50)), "line 1"),
    ],
)
def test_malformed_rows_name_the_line(line: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_tabular([line], INDEX)


def test_line_numbers_count_from_one_including_blank_lines() -> None:
    with pytest.raises(ValueError, match="line 3"):
        parse_tabular([row(("r000001", "t000001"), (1, 50), (1, 50)), "", "bad"], INDEX)


# Versions


@pytest.mark.parametrize(
    ("text", "version"),
    [("diamond version 2.2.9\n", "2.2.9"), ("diamond version 2.0.0", "2.0.0")],
)
def test_supported_versions(text: str, version: str) -> None:
    assert parse_version(text) == version


@pytest.mark.parametrize("text", ["diamond version 0.9.36", "diamond version 1.99.0"])
def test_old_versions_are_rejected(text: str) -> None:
    with pytest.raises(DiamondError, match="2.0 or newer"):
        parse_version(text)


def test_unreadable_version_is_rejected() -> None:
    with pytest.raises(DiamondError, match="version"):
        parse_version("something else")
