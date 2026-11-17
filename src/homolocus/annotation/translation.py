"""CDS translation rules (spec 3.2). Pure functions, no I/O."""

import re
from dataclasses import dataclass

from Bio.Data import CodonTable
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

_VALID_SUPPLIED = re.compile(r"[A-Za-z*]+")


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
    table: int = 11,
    codon_start: int = 1,
    complete_5prime: bool = True,
) -> Translation:
    """Translate a CDS given in coding orientation (spec 3.2)."""
    if codon_start not in (1, 2, 3):
        return Untranslatable(f"invalid codon_start {codon_start}")
    if table not in CodonTable.unambiguous_dna_by_id:
        return Untranslatable(f"unknown translation table {table}")
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


def _is_start_codon(codon: str, table: int) -> bool:
    start_codons: list[str] = CodonTable.unambiguous_dna_by_id[table].start_codons
    return codon in start_codons


def _classify(protein: str, *, trimmed: bool) -> Translation:
    if not protein:
        return Untranslatable("no complete codon before the stop")
    if "*" in protein:
        return InternalStop()
    return Protein(protein, trimmed=trimmed)
