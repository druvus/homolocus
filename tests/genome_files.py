"""GenBank and FASTA input files with planted genes, for CLI and README tests."""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def planted(dna, orfs):
    """Sequence of ``orfs`` separated by 60 random bases, and their CDS features."""
    parts, features, position = [dna.random(60)], [], 60
    for orf in orfs:
        features.append(SeqFeature(FeatureLocation(position, position + len(orf), 1), type="CDS"))
        parts += [orf, dna.random(60)]
        position += len(orf) + 60
    return "".join(parts), features


def genbank(path, record_id, dna, make_record, orfs=None):
    """A GenBank record of planted ORFs (three random ones by default)."""
    sequence, features = planted(dna, orfs or [dna.orf(150) for _ in range(3)])
    SeqIO.write([make_record(record_id, sequence, features)], path, "genbank")
    return path


def fasta(path, contig_id, dna, orfs, flank=3000):
    """A one-contig FASTA with ``orfs`` planted between ``flank`` random bases."""
    sequence, _ = planted(dna, orfs)
    contig = dna.random(flank) + sequence + dna.random(flank)
    SeqIO.write([SeqRecord(Seq(contig), id=contig_id, description="")], path, "fasta")
    return path
