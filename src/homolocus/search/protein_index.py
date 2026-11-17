"""Short tokens for DIAMOND input and their mapping to gene IDs (spec 4.1, 6.5).

Reference proteins get ``r000001``, ... and target proteins ``t000001``, ...
Tokens are looked up per side, so the same gene may appear on both sides.
Genes without a protein are not searchable and get no token.
"""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

from homolocus.models.gene import Gene


@dataclass(frozen=True)
class ProteinIndex:
    """Token <-> gene ID mapping and the protein sequences sent to DIAMOND."""

    reference_proteins: Mapping[str, str]
    target_proteins: Mapping[str, str]
    gene_ids: Mapping[str, str]

    @classmethod
    def build(cls, reference_genes: Iterable[Gene], target_genes: Iterable[Gene]) -> "ProteinIndex":
        """Assign tokens, in input order, to every gene that has a protein."""
        gene_ids: dict[str, str] = {}
        reference = _tokenise("r", reference_genes, gene_ids)
        target = _tokenise("t", target_genes, gene_ids)
        return cls(reference, target, gene_ids)

    @property
    def n_reference(self) -> int:
        """Number of searchable reference proteins."""
        return len(self.reference_proteins)

    @property
    def n_target(self) -> int:
        """Number of searchable target proteins."""
        return len(self.target_proteins)

    def gene_id_for(self, token: str) -> str:
        """Canonical gene ID of a token. Raises KeyError for an unknown token."""
        return self.gene_ids[token]

    def reference_fasta(self) -> str:
        """Reference proteins as FASTA text keyed by token."""
        return _fasta(self.reference_proteins)

    def target_fasta(self) -> str:
        """Target proteins as FASTA text keyed by token."""
        return _fasta(self.target_proteins)

    def mapping_rows(self) -> tuple[tuple[str, str], ...]:
        """(token, gene ID) pairs: reference tokens first, then target tokens."""
        return tuple(self.gene_ids.items())

    def write_mapping(self, path: str | Path) -> None:
        """Write the token -> gene ID table (``alignments/protein_ids.tsv``)."""
        rows = [
            "token\tgene_id",
            *(f"{token}\t{gene_id}" for token, gene_id in self.mapping_rows()),
        ]
        Path(path).write_text("\n".join(rows) + "\n", encoding="utf-8")


def _tokenise(prefix: str, genes: Iterable[Gene], gene_ids: dict[str, str]) -> dict[str, str]:
    proteins: dict[str, str] = {}
    for gene in genes:
        if gene.protein_sequence is None:
            continue
        token = f"{prefix}{len(proteins) + 1:06d}"
        proteins[token] = gene.protein_sequence
        gene_ids[token] = gene.id
    return proteins


def _fasta(proteins: Mapping[str, str]) -> str:
    return "".join(f">{token}\n{sequence}\n" for token, sequence in proteins.items())
