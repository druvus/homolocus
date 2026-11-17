"""Reference region model (spec 2.1, 5.4)."""

from dataclasses import dataclass
from pathlib import Path

from homolocus.coordinates import to_one_based
from homolocus.models.gene import Gene, sorted_genes
from homolocus.models.skipped_feature import SkippedFeature


@dataclass(frozen=True)
class ReferenceRegion:
    """The reference interval whose proteins define what is searched for.

    ``start`` and ``end`` are 0-based half-open on the original contig. Genes
    keep their contig coordinates (not rebased) and lie fully inside the
    region. ``skipped_features`` lists excluded CDS inside the region.
    """

    contig_id: str
    start: int
    end: int
    sequence: str
    genes: tuple[Gene, ...] = ()
    name: str | None = None
    source_genome: str | None = None
    skipped_features: tuple[SkippedFeature, ...] = ()

    def __post_init__(self) -> None:
        if not self.contig_id:
            raise ValueError("reference contig ID must not be empty")
        to_one_based(self.start, self.end)
        if len(self.sequence) != self.end - self.start:
            raise ValueError(
                f"reference sequence length {len(self.sequence)} does not match "
                f"region length {self.end - self.start}"
            )
        self._check_skipped_features()
        genes = sorted_genes(
            self.genes,
            contig_id=self.contig_id,
            start=self.start,
            end=self.end,
            label="the reference region",
        )
        object.__setattr__(self, "genes", genes)

    def _check_skipped_features(self) -> None:
        for skipped in self.skipped_features:
            inside = self.start <= skipped.start and skipped.end <= self.end
            if skipped.contig_id != self.contig_id or not inside:
                raise ValueError("skipped feature lies outside the reference region")

    @classmethod
    def from_genbank(
        cls,
        path: str | Path,
        contig: str | None = None,
        start: int | None = None,
        end: int | None = None,
    ) -> "ReferenceRegion":
        """Read a reference region from a GenBank file (spec 2.1).

        ``start`` and ``end`` are 1-based inclusive and must be given together.
        Without them the whole record is the region. ``contig`` selects the
        record and may be omitted when the file holds a single record.
        """
        from homolocus.annotation.genbank import read_reference_region  # noqa: PLC0415

        return read_reference_region(path, contig=contig, start=start, end=end)
