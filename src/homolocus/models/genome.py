"""Genome model (spec 5.3)."""

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from homolocus.models.contig import Contig
from homolocus.models.gene import Gene


@dataclass(frozen=True)
class Genome:
    """One assembly: every contig, in input order, keyed by contig ID."""

    contigs: Mapping[str, Contig]
    name: str | None = None
    source: str | None = None
    _genes_by_id: Mapping[str, Gene] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not self.contigs:
            raise ValueError("a genome needs at least one contig")
        for key, contig in self.contigs.items():
            if key != contig.id:
                raise ValueError(f"contig key {key!r} does not match contig ID {contig.id!r}")
        object.__setattr__(self, "contigs", MappingProxyType(dict(self.contigs)))
        object.__setattr__(self, "_genes_by_id", {gene.id: gene for gene in self.genes()})

    def genes(self) -> Iterator[Gene]:
        """Iterate over all genes, contig by contig, each in start order."""
        for contig in self.contigs.values():
            yield from contig.genes

    def gene(self, gene_id: str) -> Gene:
        """Return the gene with this canonical ID. Raises ``KeyError`` if absent."""
        return self._genes_by_id[gene_id]
