"""CDS translation rules (spec 3.2)."""

import pytest

from homolocus.annotation.translation import (
    InternalStop,
    InvalidSupplied,
    Protein,
    Untranslatable,
    check_supplied_translation,
    translate_cds,
)


def test_complete_cds_drops_trailing_stop() -> None:
    assert translate_cds("ATGAAACCCTAA") == Protein("MKP")


def test_translation_table_is_respected() -> None:
    assert translate_cds("ATGTGATAA", table=11) == InternalStop()
    assert translate_cds("ATGTGATAA", table=4) == Protein("MW")


@pytest.mark.parametrize(("codon_start", "prefix"), [(2, "G"), (3, "GG")])
def test_translation_starts_at_codon_start(codon_start: int, prefix: str) -> None:
    result = translate_cds(prefix + "AAACCCTAA", codon_start=codon_start)
    assert result == Protein("KP")


def test_incomplete_last_codon_is_trimmed_and_flagged() -> None:
    assert translate_cds("ATGAAACC") == Protein("MK", trimmed=True)


def test_alternative_start_codon_becomes_methionine_when_five_prime_complete() -> None:
    assert translate_cds("GTGAAATAA") == Protein("MK")
    assert translate_cds("TTGAAATAA") == Protein("MK")


def test_alternative_start_is_literal_when_five_prime_partial() -> None:
    assert translate_cds("GTGAAATAA", complete_5prime=False) == Protein("VK")


def test_alternative_start_is_literal_when_codon_start_is_not_one() -> None:
    assert translate_cds("AGTGAAATAA", codon_start=2) == Protein("VK")


def test_non_start_first_codon_is_literal() -> None:
    assert translate_cds("AAAGTGTAA") == Protein("KV")


def test_lower_case_sequence_is_translated() -> None:
    assert translate_cds("gtgaaataa") == Protein("MK")


def test_ambiguous_base_translates_to_x() -> None:
    assert translate_cds("ATGNNNTAA") == Protein("MX")


def test_internal_stop_is_reported() -> None:
    assert translate_cds("ATGTAAAAATAA") == InternalStop()


@pytest.mark.parametrize(
    ("nucleotides", "kwargs"),
    [
        ("ATGAZATAA", {}),
        ("TAA", {}),
        ("AT", {}),
        ("ATGAAATAA", {"table": 99}),
        ("ATGAAATAA", {"codon_start": 4}),
        ("ATGAAATAA", {"codon_start": 0}),
    ],
)
def test_untranslatable_input(nucleotides: str, kwargs: dict[str, int]) -> None:
    assert isinstance(translate_cds(nucleotides, **kwargs), Untranslatable)


def test_supplied_translation_is_trusted() -> None:
    assert check_supplied_translation("MKPV") == Protein("MKPV")


def test_supplied_translation_with_stop_is_internal_stop() -> None:
    assert check_supplied_translation("MK*P") == InternalStop()


@pytest.mark.parametrize("text", ["", "MK1P", "MK P"])
def test_invalid_supplied_translation(text: str) -> None:
    assert check_supplied_translation(text) == InvalidSupplied()
