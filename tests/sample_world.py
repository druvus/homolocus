"""A small synthetic run shared by the results, output and visualization tests.

Reference genes R1-R4 (one family each) and a target contig t1 with two
candidates: A, B, an inserted gene, C (families F001-F003, F004 missing)
and a weaker copy FAR_A, FAR_B far away. B also has a weak match to R4.
"""

from homolocus import (
    Alignment,
    Contig,
    DiamondSearch,
    FamilyGrouping,
    Gene,
    GeneFamily,
    GeneMatch,
    Genome,
    ReferenceFamilies,
    ReferenceRegion,
    RegionDetection,
    RegionResults,
    ScoringWeights,
)
from homolocus.detection.evaluation import DetectionContext, detect
from homolocus.search.protein_index import ProteinIndex

PROTEIN = 300
REF_GENES = tuple(
    Gene("ref", n * 1200, n * 1200 + 900, 1, "A" * 900, "M" * PROTEIN) for n in range(4)
)
REFERENCE = ReferenceRegion("ref", 0, 4500, "A" * 4500, REF_GENES)
FAMILIES = ReferenceFamilies(
    tuple(GeneFamily(f"F{n + 1:03d}", (g.id,)) for n, g in enumerate(REF_GENES))
)


def target(start: int, strand: int = 1, **metadata) -> Gene:
    return Gene("t1", start, start + 900, strand, "A" * 900, "M" * 300, **metadata)


def match(ref: Gene, tgt: Gene, *, query=(0, PROTEIN), identity=80.0, bits=300.0) -> GeneMatch:
    span = query[1] - query[0]
    hsp = Alignment(ref.id, tgt.id, identity, span, *query, 0, span, PROTEIN, PROTEIN, 1e-40, bits)
    return GeneMatch(ref.id, tgt.id, (hsp,))


A = target(10_000, locus_tag="tA", product="alpha")
B = target(11_200)
INSERT = target(12_400, strand=-1, product="insertion")
C = target(13_600)
FAR_A, FAR_B = target(100_000), target(101_200)
MATCHES = (
    match(REF_GENES[0], A),
    match(REF_GENES[1], B),
    match(REF_GENES[3], B, query=(0, 60), identity=50.0, bits=40.0),
    match(REF_GENES[2], C),
    match(REF_GENES[0], FAR_A),
    match(REF_GENES[1], FAR_B),
)
GENOME = Genome(
    contigs={"t1": Contig("t1", "A" * 200_000, (A, B, INSERT, C, FAR_A, FAR_B))},
    name="target",
    source="target.fna",
)
SETTINGS = (DiamondSearch(), FamilyGrouping(), RegionDetection(min_matches=2), ScoringWeights())


def results(**overrides) -> RegionResults:
    """The sample results; any field can be overridden."""
    search, grouping, detection, scoring = SETTINGS
    context = DetectionContext(REFERENCE, FAMILIES, search, detection, scoring)
    matches = overrides.get("matches", MATCHES)
    genome = overrides.get("genome", GENOME)
    fields = {
        "reference": REFERENCE,
        "genome": genome,
        "families": FAMILIES,
        "matches": matches,
        "candidates": detect(context, genome, matches),
        "search": search,
        "family_grouping": grouping,
        "detection": detection,
        "scoring": scoring,
        "diamond_version": "2.2.9",
        "raw_main": "r000001\tt000001\n",
        "raw_self": "",
        "protein_index": ProteinIndex.build(REF_GENES, GENOME.genes()),
    }
    fields.update(overrides)
    return RegionResults(**fields)
