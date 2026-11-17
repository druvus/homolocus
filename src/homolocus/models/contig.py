"""Contig model (spec 5.2)."""

from dataclasses import dataclass

from homolocus.models.gene import Gene


def _gene_order(gene: Gene) -> tuple[int, int, int]:
    return gene.start, gene.end, gene.strand


@dataclass(frozen=True)
class Contig:
    """One linear sequence of an assembly and the genes annotated on it.

    Genes are sorted by start coordinate, irrespective of strand, when the
    contig is built. A contig may have no genes.
    """

    id: str
    sequence: str
    genes: tuple[Gene, ...] = ()
    description: str | None = None

    def __post_init__(self) -> None:
        if not self.id:
            raise ValueError("contig ID must not be empty")
        if not self.sequence:
            raise ValueError(f"contig {self.id} has an empty sequence")
        for gene in self.genes:
            self._check_gene(gene)
        self._check_unique_ids()
        object.__setattr__(self, "genes", tuple(sorted(self.genes, key=_gene_order)))

    def _check_gene(self, gene: Gene) -> None:
        if gene.contig_id != self.id:
            raise ValueError(f"gene {gene.id} belongs to contig {gene.contig_id}, not {self.id}")
        if gene.end > self.length:
            raise ValueError(f"gene {gene.id} extends past the end of contig {self.id}")

    def _check_unique_ids(self) -> None:
        ids = [gene.id for gene in self.genes]
        if len(set(ids)) != len(ids):
            raise ValueError(f"contig {self.id} has duplicate gene IDs")

    @property
    def length(self) -> int:
        """Length of the contig in base pairs."""
        return len(self.sequence)

    @property
    def n_genes(self) -> int:
        """Number of genes on the contig."""
        return len(self.genes)
