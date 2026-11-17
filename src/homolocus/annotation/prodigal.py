"""Gene prediction with Pyrodigal (spec 2.2). Computation only, no file I/O."""

from collections.abc import Mapping

import pyrodigal

from homolocus.models.gene import Gene, as_strand

TRAINING_MINIMUM = pyrodigal.MIN_SINGLE_GENOME
"""Contigs at least this long (20 kb) use single mode; shorter ones use meta mode."""

_TABLE_11_BINS = pyrodigal.MetagenomicBins(
    [b for b in pyrodigal.METAGENOMIC_BINS if b.training_info.translation_table == 11]
)


def predict_genes(sequences: Mapping[str, str]) -> dict[str, tuple[Gene, ...]]:
    """Predict protein-coding genes on each contig, keeping input order.

    Contigs of at least 20 kb use single mode, trained once on all contigs.
    Shorter contigs use meta mode with translation-table-11 models only.
    """
    single = _single_finder(sequences)
    meta = pyrodigal.GeneFinder(meta=True, metagenomic_bins=_TABLE_11_BINS)
    predicted: dict[str, tuple[Gene, ...]] = {}
    for contig_id, sequence in sequences.items():
        finder = single if single and len(sequence) >= TRAINING_MINIMUM else meta
        genes = finder.find_genes(sequence)
        predicted[contig_id] = tuple(_gene(contig_id, gene) for gene in genes)
    return predicted


def _single_finder(sequences: Mapping[str, str]) -> pyrodigal.GeneFinder | None:
    if all(len(sequence) < TRAINING_MINIMUM for sequence in sequences.values()):
        return None
    finder = pyrodigal.GeneFinder()
    finder.train(*sequences.values())
    return finder


def _gene(contig_id: str, predicted: pyrodigal.Gene) -> Gene:
    return Gene(
        contig_id=contig_id,
        start=predicted.begin - 1,
        end=predicted.end,
        strand=as_strand(predicted.strand),
        nucleotide_sequence=predicted.sequence(),
        protein_sequence=predicted.translate(include_stop=False),
        partial=predicted.partial_begin or predicted.partial_end,
    )
