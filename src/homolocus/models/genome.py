"""Genome model (spec 5.3)."""

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING

from homolocus.models.contig import Contig
from homolocus.models.gene import Gene
from homolocus.models.skipped_feature import SkippedFeature

if TYPE_CHECKING:
    from homolocus.annotation.genome_file import Annotation


@dataclass(frozen=True)
class Genome:
    """One assembly: every contig, in input order, keyed by contig ID.

    ``skipped_features`` lists CDS features excluded while parsing (spec 5.10).
    ``display_name`` is an optional presentation name; ``label`` shows it, or
    ``name`` without one (spec 5.3). ``parse_warnings`` keeps the text of
    each per-contig GenBank parse warning (spec 3.2, 14).
    """

    contigs: Mapping[str, Contig]
    name: str | None = None
    source: str | None = None
    skipped_features: tuple[SkippedFeature, ...] = ()
    display_name: str | None = None
    parse_warnings: tuple[str, ...] = ()
    _genes_by_id: Mapping[str, Gene] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not self.contigs:
            raise ValueError("a genome needs at least one contig")
        if self.display_name == "":
            raise ValueError("a genome's display name must not be empty; use None")
        for key, contig in self.contigs.items():
            if key != contig.id:
                raise ValueError(f"contig key {key!r} does not match contig ID {contig.id!r}")
        object.__setattr__(self, "contigs", MappingProxyType(dict(self.contigs)))
        object.__setattr__(self, "_genes_by_id", {gene.id: gene for gene in self.genes()})

    @property
    def label(self) -> str | None:
        """The name to show: the display name when set, else the name."""
        return self.display_name or self.name

    def __repr__(self) -> str:
        return (
            f"Genome(name={self.name!r}, label={self.label!r}, "
            f"contigs={len(self.contigs)}, genes={len(self._genes_by_id)})"
        )

    def genes(self) -> Iterator[Gene]:
        """Iterate over all genes, contig by contig, each in start order."""
        for contig in self.contigs.values():
            yield from contig.genes

    def gene(self, gene_id: str) -> Gene:
        """Return the gene with this canonical ID. Raises ``KeyError`` if absent."""
        return self._genes_by_id[gene_id]

    @classmethod
    def from_file(cls, path: str | Path, annotation: "Annotation" = "auto") -> "Genome":
        """Read a target genome from a GenBank or FASTA file (spec 2.2).

        ``"auto"`` uses supplied CDS annotation when a GenBank file has any,
        and Pyrodigal otherwise. ``"pyrodigal"`` always predicts genes.
        """
        from homolocus.annotation.genome_file import read_genome  # noqa: PLC0415

        return read_genome(path, annotation)
