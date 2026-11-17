"""CDS translation rules (spec 3.2). Pure functions, no I/O."""

import re
from dataclasses import dataclass
from typing import Literal, NewType

from Bio.Data import CodonTable
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

_VALID_SUPPLIED = re.compile(r"[A-Za-z*]+")

CodonStart = Literal[1, 2, 3]
"""Position of the first base of the first codon (``/codon_start``)."""

TranslationTable = NewType("TranslationTable", int)
"""An NCBI translation table ID known to Biopython. Build it with parse_translation_table."""

TABLE_11 = TranslationTable(11)

_CODON_STARTS: dict[str, CodonStart] = {"1": 1, "2": 2, "3": 3}


def parse_codon_start(text: str | None) -> CodonStart | None:
    """Parse ``/codon_start``: 1 when absent, None when not 1, 2 or 3."""
    if text is None:
        return 1
    return _CODON_STARTS.get(text)


def parse_translation_table(text: str | None) -> TranslationTable | None:
    """Parse ``/transl_table``: table 11 when absent, None when not a known table."""
    if text is None:
        return TABLE_11
    if not text.isdigit() or int(text) not in CodonTable.unambiguous_dna_by_id:
        return None
    return TranslationTable(int(text))


@dataclass(frozen=True)
class Protein:
    """A usable protein. ``trimmed`` is set when an incomplete last codon was dropped."""

    sequence: str
    trimmed: bool = False


@dataclass(frozen=True)
class InternalStop:
    """The translation contains a stop before its end; the CDS is a pseudogene."""


@dataclass(frozen=True)
class Untranslatable:
    """No protein could be derived."""

    reason: str


@dataclass(frozen=True)
class InvalidSupplied:
    """A supplied ``/translation`` is empty or malformed and must be ignored."""


Translation = Protein | InternalStop | Untranslatable
SuppliedTranslation = Protein | InternalStop | InvalidSupplied


def translate_cds(
    nucleotides: str,
    *,
    table: TranslationTable = TABLE_11,
    codon_start: CodonStart = 1,
    complete_5prime: bool = True,
) -> Translation:
    """Translate a CDS given in coding orientation (spec 3.2).

    ``table`` and ``codon_start`` are already valid; parse outside values with
    parse_translation_table and parse_codon_start.
    """
    coding = nucleotides[codon_start - 1 :].upper()
    whole_codons = coding[: len(coding) - len(coding) % 3]
    try:
        protein = str(Seq(whole_codons).translate(table=table))
    except TranslationError as error:
        return Untranslatable(str(error))
    if complete_5prime and codon_start == 1 and _is_start_codon(whole_codons[:3], table):
        protein = "M" + protein[1:]
    return _classify(protein.removesuffix("*"), trimmed=whole_codons != coding)


def check_supplied_translation(text: str) -> SuppliedTranslation:
    """Classify a ``/translation`` qualifier value (spec 3.2)."""
    if not _VALID_SUPPLIED.fullmatch(text):
        return InvalidSupplied()
    if "*" in text:
        return InternalStop()
    return Protein(text)


def _is_start_codon(codon: str, table: TranslationTable) -> bool:
    start_codons: list[str] = CodonTable.unambiguous_dna_by_id[table].start_codons
    return codon in start_codons


def _classify(protein: str, *, trimmed: bool) -> Translation:
    if not protein:
        return Untranslatable("no complete codon before the stop")
    if "*" in protein:
        return InternalStop()
    return Protein(protein, trimmed=trimmed)
