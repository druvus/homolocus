"""Contig model (spec 5.2)."""

from dataclasses import dataclass

from homolocus.models.gene import Gene, sorted_genes


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
        genes = sorted_genes(
            self.genes, contig_id=self.id, start=0, end=self.length, label=f"contig {self.id}"
        )
        object.__setattr__(self, "genes", genes)

    @property
    def length(self) -> int:
        """Length of the contig in base pairs."""
        return len(self.sequence)

    @property
    def n_genes(self) -> int:
        """Number of genes on the contig."""
        return len(self.genes)
