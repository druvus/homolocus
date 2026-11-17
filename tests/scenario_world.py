"""Synthetic reference regions and targets for real-DIAMOND scenarios (spec 1, 17).

Genes are laid out from proteins: each gene is ``3 * len(protein) + 3`` bp
(stop codon included) with ``GAP`` bp between neighbours. Nucleotide
sequences are placeholders; DIAMOND only sees the proteins.
"""

from collections.abc import Sequence
from dataclasses import dataclass

from conftest import SyntheticProteins
from homolocus import (
    Contig,
    DiamondSearch,
    FamilyGrouping,
    Gene,
    Genome,
    ReferenceRegion,
    RegionDetection,
    RegionFinder,
    RegionResults,
)

LENGTH = 300
GAP = 300
OFFSET = 10_000


@dataclass(frozen=True)
class Slot:
    """One target gene: its protein and strand."""

    protein: str
    strand: int = 1


def lay_out(contig_id: str, slots: Sequence[Slot], *, offset: int) -> tuple[Gene, ...]:
    """Genes for ``slots`` placed left to right from ``offset``."""
    genes, position = [], offset
    for slot in slots:
        length = 3 * len(slot.protein) + 3
        genes.append(
            Gene(contig_id, position, position + length, slot.strand, "A" * length, slot.protein)
        )
        position += length + GAP
    return tuple(genes)


class World:
    """A reference region of random proteins and ways to derive target genes."""

    def __init__(self, n_genes: int, *, seed: int = 4) -> None:
        self._proteins = SyntheticProteins(seed)
        self.reference_proteins = [self._proteins.random(LENGTH) for _ in range(n_genes)]

    def reference(self) -> ReferenceRegion:
        """The reference region on contig ``ref``, genes in protein order."""
        genes = lay_out("ref", [Slot(p) for p in self.reference_proteins], offset=0)
        end = genes[-1].end
        return ReferenceRegion("ref", 0, end, "A" * end, genes)

    def copy(self, index: int, identity: float = 0.8) -> str:
        """A homologue of reference protein ``index`` at about ``identity``."""
        return self._proteins.mutate(self.reference_proteins[index], identity)

    def part(self, index: int, start: int, end: int, identity: float = 0.8) -> str:
        """A homologue of residues ``start:end`` of reference protein ``index``."""
        return self._proteins.mutate(self.reference_proteins[index][start:end], identity)

    def unrelated(self, length: int = LENGTH) -> str:
        """A random protein with no homologue in the reference."""
        return self._proteins.random(length)

    def contig(
        self, contig_id: str, slots: Sequence[Slot], *, offset: int = OFFSET, tail: int = OFFSET
    ) -> Contig:
        """A contig holding ``slots`` after ``offset`` bp, with ``tail`` bp after them."""
        genes = lay_out(contig_id, slots, offset=offset)
        length = (genes[-1].end if genes else offset) + tail
        return Contig(contig_id, "A" * length, genes)


def genome(*contigs: Contig) -> Genome:
    """A target genome of the given contigs."""
    return Genome(contigs={c.id: c for c in contigs}, name="target", source="target.fna")


def run(
    reference: ReferenceRegion,
    target: Genome,
    *,
    detection: RegionDetection | None = None,
    family_grouping: FamilyGrouping | None = None,
) -> RegionResults:
    """RegionFinder with real DIAMOND on one thread and default thresholds."""
    finder = RegionFinder(
        reference,
        search=DiamondSearch(threads=1),
        family_grouping=family_grouping or FamilyGrouping(),
        detection=detection or RegionDetection(),
    )
    return finder.find(target)
