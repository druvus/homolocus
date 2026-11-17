"""Paralogous reference gene families (spec 5.5, 7.6)."""

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field


@dataclass(frozen=True)
class GeneFamily:
    """Reference genes grouped as paralogs. Members are in genomic order."""

    id: str
    gene_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.id:
            raise ValueError("family ID must not be empty")
        if not self.gene_ids:
            raise ValueError(f"family {self.id} needs at least one gene")
        if len(set(self.gene_ids)) != len(self.gene_ids):
            raise ValueError(f"family {self.id} has duplicate genes")

    @property
    def representative_id(self) -> str:
        """The earliest member in the reference, used for ordering."""
        return self.gene_ids[0]


@dataclass(frozen=True)
class ReferenceFamilies:
    """All families of a reference region; each searchable gene is in exactly one.

    ``len()`` is the number of families, the reference-coverage denominator.
    """

    families: tuple[GeneFamily, ...]
    _by_gene: Mapping[str, GeneFamily] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if len({family.id for family in self.families}) != len(self.families):
            raise ValueError("duplicate family ID")
        by_gene: dict[str, GeneFamily] = {}
        for family in self.families:
            for gene_id in family.gene_ids:
                if gene_id in by_gene:
                    raise ValueError(f"gene {gene_id} is in more than one family")
                by_gene[gene_id] = family
        object.__setattr__(self, "_by_gene", by_gene)

    def family_of(self, gene_id: str) -> GeneFamily:
        """The family of a searchable reference gene. Raises KeyError otherwise."""
        return self._by_gene[gene_id]

    def __len__(self) -> int:
        return len(self.families)

    def __iter__(self) -> Iterator[GeneFamily]:
        return iter(self.families)
