# homolocus — v0.1 Specification

- Status: Working specification (revision 2)
- Package name: `homolocus`
- Python: 3.10+
- Build backend: Hatchling
- Primary use case: Identify and visualize a homologous genomic region in a bacterial genome using the proteins encoded by a reference region.

## 0. Decision status and conventions

Every behaviour-affecting statement carries one of three tags where it matters:

- **[DECIDED]** Confirmed by the user.
- **[PROPOSED]** Default suggested during design, not yet explicitly confirmed. Implement as written, but keep it easy to change, and list it in the first review.
- **[OPEN]** Not settled. Do not resolve silently. Record a proposed default, and ask when the choice materially affects behaviour.

Untagged text carries over from revision 1, which the user accepted.

### 0.1 Decision log (revision 2)

| # | Topic | Decision |
|---|-------|----------|
| 1 | Duplicated target homologues in candidate generation | **[DECIDED]** Each copy is a separate matching event for distance calculations. All copies are retained in reconstruction and scoring. |
| 2 | Where coverage thresholds apply | **[DECIDED]** Two levels. Pair level filters on e-value and identity only. Coverage is evaluated at reference-gene level as a union over target proteins. |
| 3 | User-facing coordinates | **[DECIDED]** 1-based, inclusive in CLI, API inputs, TSV and plots. Internally 0-based half-open. Conversion only in the I/O layer. |
| 4 | Candidate grouping | **[DECIDED]** Single-linkage chains on the distance limits. Chains longer than `max_region_bp` are split into overlapping windows, then deduplicated. |
| 5 | Gene ID scheme | **[DECIDED]** `{contig_id}\|{start_1based}-{end_1based}\|{strand_char}`, e.g. `NC_000913.3\|1900001-1901200\|+`. `locus_tag` is metadata only. |
| 6 | Circular genomes and compound CDS | **[DECIDED]** Not supported in v0.1. Contigs are treated as linear. CDS with compound (`join`) locations are skipped with a warning. |
| 7 | `auto` annotation | **[DECIDED]** GenBank with at least one CDS is used as is. GenBank without CDS and FASTA use Pyrodigal. Contigs under 20 kb use meta mode. Partial genes at contig edges are flagged. |
| 8 | Translation errors | **[DECIDED]** See 3.2. |
| 9 | Paralogous reference genes | **[DECIDED]** Reference genes in one family (self-search, identity at or above 90 %) count as one unit in reference coverage. A fusion requires distinct families. |
| 10 | Scoring | **[DECIDED]** Simple weighted sum with visible components. Starting weights 0.40 / 0.30 / 0.15 / 0.10 / 0.05, configurable and not validated. |

### 0.2 Units

- **[PROPOSED]** Identity is a percentage (0–100), as in DIAMOND output. All coverage values are fractions (0–1).
- Earlier examples used `min_query_coverage=50`. Those parameters are replaced as described in 6.2.

────────

## 1. Purpose and scope

The package takes a reference genomic region and a target genome assembly. It identifies proteins in the target genome homologous to proteins encoded in the reference region, then detects candidate genomic regions containing a related collection of proteins.

The approach must not depend on predefined flanking genes. The reference region is characterized primarily by its protein content. Gene order and orientation provide supporting evidence, but neither is an absolute requirement.

The package should tolerate:

- Missing reference genes.
- Insertions of unrelated genes.
- Gene duplications and multiple homologues.
- Rearrangements in gene order.
- Reverse-oriented versions of the region in the target genome.
- Potential split proteins and fused proteins (detected and flagged in v0.1).
- Repeated reference gene families such as IS-element transposases (see 7.6).

The principal example is an E. coli pathogenicity island (PAI), but the package is intended for bacterial genomes more generally.

### v0.1 goals

- Parse reference and target genome input.
- Annotate target FASTA input with Pyrodigal when needed.
- Group paralogous reference proteins into families.
- Search reference proteins against target proteins using DIAMOND.
- Generate candidate regions independently on each contig.
- Rank all candidates passing configurable detection thresholds.
- Produce Matplotlib plots, including multiple candidates as separate lines.
- Provide both a Python API and a CLI.

### Explicitly out of scope for v0.1

- GFF input.
- Non-CDS feature modelling (RNA genes, etc.).
- Pseudogene retention.
- CDS with compound (`join`) locations, and features wrapping the origin of a circular contig.
- Combining matches from different contigs into one candidate.
- Interactive plotting.
- Saving and reloading a serialized analysis/results object.
- Full reconstruction of split/fused protein architecture.
- General assembly QC statistics such as N50 and GC content.

────────

## 2. Input formats

### 2.1 Reference input

Primary mode: a complete GenBank genome record plus a contig identifier and region coordinates.

```python
reference = ReferenceRegion.from_genbank(
    "ecoli.gb",
    contig="NC_000913.3",
    start=1_900_000,   # 1-based, inclusive  [DECIDED]
    end=1_920_000,     # 1-based, inclusive
)
```

Alternative mode: a GenBank file containing only the reference region. The record sequence boundaries define the region.

```python
reference = ReferenceRegion.from_genbank("ecoli_pai.gbk")
```

Rules:

- Reference region boundaries are independent of gene boundaries.
- Only CDS features completely contained within the requested region are included. A CDS overlapping either region boundary is excluded.
- CDS with compound locations are skipped with a warning and listed in the run metadata.
- The reference sequence and gene coordinates retain the original contig orientation. The package does not reverse-complement or reorient the reference region.
- Reference genes with `protein_sequence=None` stay in the model but cannot be searched. They are excluded from the reference-coverage denominator and reported separately (8.1).
- **[PROPOSED]** The reference file may be gzip-compressed, as for target genomes (2.2).

### 2.2 Target genome input

The target is one assembly per analysis. The assembly may contain one or more contigs.

- **[PROPOSED]** A run may take several target assemblies. Each is searched and detected on its own; the reference families and the self-search are computed once and shared. A genome's results are identical to a run with that genome alone.

Supported inputs for v0.1:

- GenBank (`.gb`, `.gbk`): use the supplied CDS annotations.
- FASTA (`.fa`, `.fasta`, `.fna`): run Pyrodigal to predict protein-coding genes.

```python
genome = Genome.from_file("target.gbk", annotation="auto")
genome = Genome.from_file("target.fna", annotation="auto")
```

`auto` behaviour **[DECIDED]**:

- GenBank input with at least one CDS feature: use the supplied annotation. Do not run Pyrodigal.
- GenBank input with no CDS features, and all FASTA input: run Pyrodigal.
- Pyrodigal mode: contigs shorter than 20 kb use meta mode. Training details for longer contigs (single genome training on the whole assembly vs. per contig) are **[PROPOSED]**: when at least one contig is 20 kb or longer, train once on all contigs of the assembly. Translation table 11.
- Partial genes at contig edges (Pyrodigal partial flags, or `<`/`>` in GenBank locations) are kept and flagged (`Gene.partial`).

Input details:

- **[DECIDED]** The file format is detected from the content, not the extension: a first non-empty line starting with `LOCUS` is GenBank, one starting with `>` is FASTA, anything else is an error.
- **[PROPOSED]** Any input file may be gzip-compressed. Compression is detected from the content (the gzip magic bytes `1f 8b`), not the extension. The file is decompressed while it is read, and format detection then applies to the decompressed text. Other compression formats are not supported.
- **[DECIDED]** `annotation` accepts `"auto"` (above) and `"pyrodigal"`, which always predicts genes and ignores any supplied CDS annotation.
- **[PROPOSED]** Meta mode uses only Pyrodigal's translation-table-11 models, so table 11 holds in meta mode too.
- **[PROPOSED]** Training happens only when at least one contig is 20 kb or longer; it then uses all contigs of the assembly. Contigs under 20 kb use meta mode either way.

GFF is not supported in v0.1.

### 2.3 Multiple contigs

`Genome` represents an assembly containing one or more contigs.

All contigs must be preserved, including contigs with no CDS. Candidate generation and region detection operate independently on each contig. Matches on different contigs are never combined into one candidate in v0.1.

All contigs are treated as linear **[DECIDED]**.

────────

## 3. Annotation rules

### 3.1 CDS-only model

For v0.1, `Gene` represents protein-coding CDS features only.

- RNA genes and other feature types are not represented.
- Explicitly annotated pseudogenes are excluded entirely (`/pseudo`, `/pseudogene`).
- Excluded pseudogenes are not in protein FASTA files, DIAMOND searches, gene counts or scoring.
- CDS with compound locations are skipped with a warning **[DECIDED]**.

### 3.2 Protein translation

For GenBank CDS features:

- Use `/translation` when available.
- If `/translation` is absent, derive the protein from the CDS nucleotide sequence.
- Respect `/transl_table` (default: bacterial table 11) and `/codon_start` (default 1).
- Handle reverse-strand CDS correctly.
- Use the CDS nucleotide sequence in coding orientation (5′ to 3′).

Error policy **[DECIDED]**:

- Missing or invalid `/translation` followed by a failed derived translation: `protein_sequence=None`, a warning, and the gene is recorded as untranslatable. Such a gene cannot be searched.
- Internal stop codon in a translation: warning, and the CDS is treated as a pseudogene and excluded.
- **[PROPOSED]** The internal-stop check applies to a translation derived by the package. A supplied `/translation` is trusted, except that one containing `*` is handled like an internal stop. `/transl_except` is not interpreted in v0.1.
- **[PROPOSED]** A supplied `/translation` is invalid when it is empty or contains characters other than letters and `*`. An invalid one is ignored and the protein is derived instead.

Derived translations **[DECIDED]**:

- Translation starts at `/codon_start`. If the remaining length is not a multiple of 3, the incomplete last codon is dropped, the rest is translated, and a warning is issued. This applies to complete and partial CDS alike.
- An alternative start codon of the translation table (e.g. GTG, TTG) is translated as `M` when the CDS is complete at its 5′ end and `/codon_start` is 1. Otherwise the first codon is translated literally.
- One trailing stop is removed. Any other `*` is an internal stop.

Reporting:

- **[DECIDED]** Features excluded as pseudogenes (explicit or by internal stop) and compound CDS are recorded as `SkippedFeature` entries on the resulting `Genome` or `ReferenceRegion` (5.10), in addition to a warning. Untranslatable genes stay in the model with `protein_sequence=None`.
- **[PROPOSED]** Warnings are summarised as one per contig, with counts by reason (skipped, untranslatable, trimmed). The details are in `skipped_features`.

For FASTA input, Pyrodigal supplies the predicted CDS coordinates and protein translations.

Tests must cover internal stops, partial CDS features and non-standard `/codon_start`.

────────

## 4. Coordinate and orientation conventions

Internal models **[DECIDED]**:

- Coordinates are 0-based, half-open: `[start, end)`.
- `strand` is +1 or -1.
- Gene coordinates refer to the original contig.
- Gene nucleotide sequences are stored in coding orientation.
- Contig sequences remain in their original genomic orientation.
- Genes within a contig are sorted by genomic start coordinate, irrespective of strand.

User-facing coordinates **[DECIDED]**:

- Everything the user types or reads is 1-based and inclusive: `ReferenceRegion.from_genbank(start=, end=)`, CLI `--start/--end`, TSV columns, plot axes and labels.
- Conversion happens only in the I/O layer. Helper functions live in one module (`coordinates.py`) with tests.
- DIAMOND coordinates (1-based inclusive) are converted to the internal convention when parsed. The raw DIAMOND file is preserved unchanged.

The region detection algorithm must find the reference region even if its orientation in the target genome is reversed. Do not require matching absolute strand orientation.

### 4.1 Gene IDs [DECIDED]

- Canonical ID: `{contig_id}|{start_1based}-{end_1based}|{strand_char}`, with `+` or `-`.
- `locus_tag`, `gene_name` and `protein_id` are metadata and never used as keys.
- **[PROPOSED]** DIAMOND FASTA headers use short sequential tokens, with a stored mapping back to canonical IDs. This avoids problems with special characters in contig IDs. Reference proteins get `r000001`, `r000002`, …, and target proteins `t000001`, …; one index serves both searches.
- **[DECIDED]** The mapping is written to `alignments/protein_ids.tsv` (columns `token`, `gene_id`), so the unchanged raw DIAMOND files can be read (6.4, 14).

────────

## 5. Data models

Exact field names and validation can be refined during implementation.

- **[PROPOSED]** Printing a model never prints sequences or raw search output. `Gene` and `ReferenceRegion` leave `nucleotide_sequence`, `protein_sequence`, `qualifiers` and `sequence` out of their repr. Collections print a one-line summary: `Contig` its ID, length and number of genes; `Genome` its name, label, contigs and genes; `ReferenceRegion` its contig, 1-based span and genes; `RegionResults` the genome label, families, candidates and best score; `GenomeSetResults` the genomes and total candidates. Equality is unchanged.

### 5.1 Gene

```python
@dataclass
class Gene:
    id: str                       # canonical ID, see 4.1
    contig_id: str

    start: int                    # 0-based, inclusive
    end: int                      # exclusive
    strand: int

    nucleotide_sequence: str      # coding orientation
    protein_sequence: str | None

    gene_name: str | None = None
    product: str | None = None
    locus_tag: str | None = None
    protein_id: str | None = None
    partial: bool = False         # partial at a contig edge or marked partial
    qualifiers: dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.end - self.start
```

### 5.2 Contig

```python
@dataclass
class Contig:
    id: str
    sequence: str
    genes: list[Gene]
    description: str | None = None
```

Preserve contigs without genes. Sort genes by start. Provide `length` and `n_genes`.

### 5.3 Genome

```python
@dataclass
class Genome:
    contigs: dict[str, Contig]
    name: str | None = None
    source: str | None = None
    skipped_features: tuple[SkippedFeature, ...] = ()
    parse_warnings: tuple[str, ...] = ()   # GenBank parse warnings, one per contig (3.2)
```

Preserve all contigs. Provide convenient access to genes across contigs. No assembly QC statistics.

- **[PROPOSED]** `Genome.from_file` sets `name` to the file name without its extension. For a compressed file, a trailing `.gz` (any case) is removed first, so `strain.fna.gz` gives `strain`. Within one run, genome names must be unique and usable as a directory name (not empty, no `/`, not starting with `.`).
- **[PROPOSED]** `Genome.display_name` is an optional presentation name (never empty). `Genome.label` is the display name when set, else `name`. Labels are what plots, the combined summary, stdout and `analysis.json` show; `name` stays the key for output folders and lookups. Display names, when given, must be unique within a run.
- **[PROPOSED]** `parse_warnings` keeps the text of each per-contig GenBank parse warning (3.2) on the `Genome` or `ReferenceRegion` it was raised for, so it reaches `analysis.json` (14). Genomes annotated by Pyrodigal have none.

### 5.4 ReferenceRegion

```python
@dataclass
class ReferenceRegion:
    contig_id: str
    start: int                    # internal 0-based half-open
    end: int
    sequence: str
    genes: list[Gene]
    name: str | None = None
    source_genome: str | None = None
    skipped_features: tuple[SkippedFeature, ...] = ()  # those inside the region
    parse_warnings: tuple[str, ...] = ()   # GenBank parse warnings, one per contig (3.2)
```

Preserve original contig coordinates and orientation. Gene coordinates are not rebased. Include only fully contained CDS. Support both input modes.

### 5.5 GeneFamily (new)

```python
@dataclass
class GeneFamily:
    id: str                       # e.g. "F001"
    gene_ids: list[str]           # reference gene IDs
    representative_id: str        # earliest gene in the reference, for ordering
    name: str | None              # from the reference annotation (7.6), never empty
```

Built from a reference self-search (7.6). Every searchable reference gene belongs to exactly one family. A gene with no paralogs is a family of one. Genes without a protein belong to no family. `gene_ids` are in genomic order and the representative is the first of them.

```python
@dataclass
class ReferenceFamilies:
    families: tuple[GeneFamily, ...]   # every searchable gene in exactly one family

    def family_of(self, gene_id: str) -> GeneFamily: ...
    def __len__(self) -> int: ...       # number of families, the reference-coverage denominator
```

### 5.6 Alignment

One DIAMOND HSP.

```python
@dataclass
class Alignment:
    query_id: str
    subject_id: str

    identity: float
    alignment_length: int

    query_start: int              # 0-based half-open, internal
    query_end: int
    subject_start: int
    subject_end: int

    query_length: int
    subject_length: int

    evalue: float
    bit_score: float
```

Keep HSPs separate. Preserve coordinates and lengths. Retain enough evidence to calculate query and subject coverage and identify discontinuous matches.

### 5.7 GeneMatch

Aggregated evidence between one reference protein and one target protein.

```python
@dataclass
class GeneMatch:
    reference_gene_id: str
    target_gene_id: str
    alignments: list[Alignment]

    query_coverage: float         # fraction of reference protein, union of HSPs
    subject_coverage: float       # fraction of target protein, union of HSPs
    identity: float               # alignment-length-weighted, percent
    bit_score: float              # sum over HSPs
```

Aggregation:

1. Collect all HSPs for a query-subject pair.
2. Union of aligned query positions. Union of aligned subject positions.
3. Query and subject coverage from those unions, without double-counting overlaps.
4. Alignment-length-weighted identity.
5. Sum of bit scores.
6. Pair-level thresholds (e-value, identity) are evaluated on the aggregated evidence. No coverage threshold is applied at pair level **[DECIDED]**.

- **[DECIDED]** The e-value of a pair is the smallest e-value of its HSPs. A pair passes when its best HSP is significant; weaker HSPs of the same pair are kept and count toward coverage.
- **[PROPOSED]** Threshold comparisons are inclusive: `evalue <= evalue threshold` and `identity >= min_identity`.
- Coverage, identity, bit score and e-value are derived from `alignments` rather than stored separately, so they cannot disagree with the HSPs.

### 5.8 Cluster

One candidate genomic region on one target contig.

Retain:

- Target contig ID and candidate start/end (internal 0-based half-open).
- Target genes in the interval, each tagged as matched or additional.
- All relevant `GeneMatch` objects.
- Families represented and missing.
- Score and `ScoreComponents` (section 8).
- Reference-gene-level coverage per family (8.1).
- Flags: potential split, potential fusion, contains partial genes, touches contig edge. The split and fusion events behind the first two flags are kept as `splits` and `fusions` (9.4).
- Number of duplicated homologue copies per family.
- **[PROPOSED]** `reversed`: whether the candidate runs opposite to the reference (direction, section 8).

The cluster is a result of detection, not a place to run DIAMOND.

### 5.9 Relationships

The model supports many-to-many relationships: one reference protein to many target proteins (duplication or split), one target protein to many reference proteins (fusion or shared family), and many HSPs per pair. No forced one-to-one mapping and no mutually exclusive `match_type`. Higher-level events are derived from alignment evidence.

### 5.10 SkippedFeature (new)

A CDS feature excluded during annotation parsing (3.1, 3.2).

```python
class SkipReason(str, Enum):
    PSEUDOGENE = "pseudogene"                # /pseudo or /pseudogene
    INTERNAL_STOP = "internal_stop"          # treated as a pseudogene
    COMPOUND_LOCATION = "compound_location"  # join(...), unsupported in v0.1

@dataclass
class SkippedFeature:
    contig_id: str
    start: int                    # internal 0-based half-open, whole location span
    end: int
    locus_tag: str | None
    reason: SkipReason
```

Reported in `analysis.json` (14).

────────

## 6. Protein homology search

### 6.1 Backend

DIAMOND is the default and only backend in v0.1. It is an external executable: check availability at runtime and give a clear error. Pyrodigal is a mandatory Python dependency. The search code lives in its own module so detection does not depend on DIAMOND-specific structures.

### 6.2 Configuration

Pair-level and reference-gene-level thresholds are separate **[DECIDED]**.

```python
search = DiamondSearch(
    evalue=1e-5,
    min_identity=30,            # percent, on aggregated GeneMatch identity
    sensitivity="sensitive",
    max_target_seqs=100,
    max_hsps=0,                 # explicit; 0 = report all HSPs per pair
    hsp_evalue=1e-3,            # per-HSP cutoff passed to DIAMOND
    threads=None,               # DIAMOND threads; None = DIAMOND's default
)
```

- **[DECIDED]** `evalue` is the pair-level threshold applied to a pair's best HSP (5.7). DIAMOND's own `--evalue` receives the looser per-HSP cutoff `hsp_evalue` (default 1e-3, never stricter than `evalue`), so weak secondary HSPs of real matches survive and count toward coverage.
- **[DECIDED]** `threads` sets the number of threads for both `diamond makedb` and `diamond blastp`. `None` passes no flag, leaving DIAMOND's default.
- **[PROPOSED]** DIAMOND 2.0 or newer is required. The installed version is checked at runtime and recorded in the run metadata.

- `min_query_coverage` and `min_subject_coverage` are removed from the search configuration.
- **Do not use DIAMOND's own `--query-cover`, `--subject-cover` or `--id` filters.** They act per HSP and would contradict the aggregation rule. Filter in Python after aggregation.
- Always pass `--max-hsps` and `-k` explicitly and verify behaviour against the installed DIAMOND version. Do not rely on defaults.
- Protein-level evidence stays distinct from region-level confidence.
- All values are provisional and not empirically validated.

### 6.3 Retaining hits and HSPs

- Retain multiple target hits per reference protein. Do not keep only the best.
- Retain multiple HSPs per query-subject pair.
- `max_target_seqs` limits target sequences reported per query. It is not the number of HSPs per pair. Record the configured limit in the run metadata. Add a test that a query with more hits than the limit behaves as documented.

### 6.4 DIAMOND tabular fields [PROPOSED]

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
```

- Raw output is written unchanged to `alignments/diamond.tsv`. It contains tokens; `alignments/protein_ids.tsv` maps them to gene IDs (4.1).
- **[PROPOSED]** Parsed alignments are sorted by query token, subject token, query start and subject start, because DIAMOND's output order can vary with threads (14, idempotence).
- Coordinates are converted from 1-based inclusive when parsed.

### 6.5 Searches performed

1. Reference proteins against target proteins (main search).
2. Reference proteins against themselves (family grouping, 7.6). Self-hits are excluded.

Genes with `protein_sequence=None` are excluded from both.

────────

## 7. Candidate generation and region detection

### 7.1 Stages

1. Family grouping of reference genes (7.6).
2. Matching events: target genes eligible to seed candidates (7.2).
3. Candidate generation: chains and windows on each contig (7.3).
4. Candidate evaluation: coverage, scoring and flags (7.7, section 8, section 9).

No predefined flanking genes. Each contig is processed independently.

### 7.2 Matching events [DECIDED rule, PROPOSED default]

A target gene is a matching event if it has at least one `GeneMatch` passing pair-level thresholds **and** that match satisfies:

```text
query_coverage >= min_event_coverage   OR   subject_coverage >= min_event_coverage
```

- The OR admits split fragments (high subject coverage, low query coverage) and fused proteins (high query coverage, low subject coverage), while excluding short shared-domain hits.
- `min_event_coverage` default **[PROPOSED]**: 0.5.
- Weak matches below the event threshold still remain in the evidence stored in the cluster, if the gene lies inside a candidate.

### 7.3 Chains and windows [DECIDED]

On each contig, take matching events sorted by genomic position.

**Chains.** Consecutive events are linked when both hold:

- Intervening genes (genes between them in the contig's sorted gene list) at most `max_intervening_genes`.
- Base-pair gap (start of the next event minus end of the previous, minimum 0 for overlapping genes) at most `max_inter_match_bp_gap`.

Linked events form a chain (single linkage). Each duplicated copy is its own event **[DECIDED]**.

**Windows.** If a chain's span (start of first event to end of last event) is at most `max_region_bp`, the chain is one candidate. Otherwise, split it into windows:

- For each event as the left anchor, extend right while the span stays at most `max_region_bp`.
- Discard windows whose set of events is a subset of another window's.

Measurement and ordering:

- **[DECIDED]** Events are sorted by start, so with nested or overlapping genes the previous event may not reach furthest. The bp gap is the next event's start minus the largest end reached so far in the chain (minimum 0), and a span runs from the first event's start to the largest end. Without nesting this equals the wording above, and a candidate always contains every event gene.
- **[PROPOSED]** Intervening genes are counted from the contig's sorted gene list, between the positions of the two consecutive events.
- **[PROPOSED]** A single event longer than `max_region_bp` still forms a window of one; `min_matches` removes it later.
- **[PROPOSED]** Windows are reported per contig, ordered by start and then end. Windows with the same event set are kept once.

**Evaluation and deduplication.**

- Each chain or window is evaluated independently (7.7).
- Candidates that fail `min_matches` or `min_reference_coverage` are discarded.
- **[PROPOSED]** After ranking, a candidate whose span overlaps a higher-ranked candidate on the same contig by more than `max_candidate_overlap` (default 0.5, fraction of the shorter span) is discarded as an alternative window of the same locus. Set to 1.0 to keep all.

### 7.4 Duplicated homologues [DECIDED]

Each copy is a separate matching event for distance calculations. Reconstruction and scoring retain all copies. For order and orientation scoring, a family uses its best copy in the candidate (highest bit score). Which copy that is must be recorded.

```text
Reference: A - B - C - D
Target:    A - B1 - X - B2 - C - D   (B1 and B2 both retained; X counted as intervening)
```

### 7.5 Detection configuration

```python
detection = RegionDetection(
    min_matches=3,                    # distinct represented families (>= 2)
    min_reference_coverage=0.3,       # fraction of families
    min_reference_gene_coverage=0.5,  # per-family combined coverage to count as represented
    min_event_coverage=0.5,
    max_intervening_genes=10,
    max_inter_match_bp_gap=20_000,
    max_region_bp=50_000,
    max_candidate_overlap=0.5,
    padding_genes=0,
    padding_bp=0,
    edge_margin_bp=500,               # contig-edge flag distance (10)
    min_split_contribution=0.2,       # split and fusion heuristics (9.4)
    max_fragment_overlap=0.2,
    min_fusion_part_coverage=0.2,
)
```

- All values are provisional and configurable.
- `min_matches` counts distinct represented families, never copies. It must be at least 2, because order and orientation need at least two points.
- **[PROPOSED]** Both `min_matches` and `min_reference_coverage` must hold (AND). Validate that `min_matches` does not exceed the number of families, with a clear error.
- Distance concepts remain distinct: intervening genes, base-pair gap between consecutive events, total span.

### 7.6 Paralogous reference families [DECIDED]

Reference genes are grouped before the main search:

- Search reference proteins against themselves with DIAMOND, excluding self-hits.
- Two genes are linked when identity is at least `family_min_identity` (default 90) **and** both coverage values are at least `family_min_coverage` (**[PROPOSED]** default 0.8, to avoid grouping genes that share only a domain).
- Families are connected components of that graph (**[PROPOSED]**).
- **[DECIDED]** DIAMOND reports A→B and B→A separately. A and B are linked when either direction passes.
- **[PROPOSED]** Identity and coverage come from the aggregated self-search `GeneMatch` (5.7). Thresholds are inclusive. No e-value condition is added beyond DIAMOND's `hsp_evalue`.
- **[PROPOSED]** Family IDs are `F001`, `F002`, … numbered in genomic order of their representatives (3 digits, more beyond 999).
- **[PROPOSED]** A family's `name` comes from the reference annotation of its members: their `/gene` qualifiers, unique and in genomic order, joined with `/` (`tnpA/tnpB`); when no member has one, their `/product` qualifiers the same way; otherwise no name. Family IDs stay the keys in every model and file; the name is a label.
- **[PROPOSED]** A `/product` made only of the words `hypothetical`, `uncharacterized`, `uncharacterised`, `conserved`, `putative` and `protein` (any case) does not name a family: "hypothetical protein", "conserved hypothetical protein" and "putative uncharacterized protein" are ignored, while "hypothetical protein YfaZ" or "type VI secretion system protein" still name one.
- Configuration: `FamilyGrouping(family_min_identity=90, family_min_coverage=0.8)`. Setting `family_min_identity=None` disables grouping, so every gene is its own family.

Effects:

- `reference_coverage`, `min_matches` and the order/orientation score count a family once.
- A target gene matching several members of one family is not a fusion.
- Family members all remain in the output. The member with the best bit score to a given target gene is used for display. All `GeneMatch` objects are kept.
- Plots colour by family (section 11).

### 7.7 Reference-gene-level coverage and "represented" [DECIDED]

For a candidate and a reference gene R:

- Take all target genes in the candidate that have a `GeneMatch` to R passing pair-level thresholds.
- Coverage of R is the union of aligned reference positions across those target proteins, divided by R's length. No double counting.
- A family is represented when the best coverage among its member genes is at least `min_reference_gene_coverage`.
- This lets several split fragments jointly represent one reference protein.
- **[DECIDED]** A candidate's evidence genes are the target genes lying fully inside its unpadded span. Padding never adds evidence (10).
- **[DECIDED]** A gene's role is `matched` when it is a matching event and `additional` otherwise. Genes with only weak, sub-event matches are `additional`; their GeneMatches are still kept and count toward 7.7 coverage.
- **[PROPOSED]** Per family, the best member is the member with the highest 7.7 coverage (ties: genomic order). The best copy is the evidence gene with the highest bit score to that best member (ties: genomic order). Similarity, order, orientation and spacing all use this pairing.

### 7.8 Gene order

Gene order is supporting evidence, not a requirement. Rearrangements do not disqualify a candidate. Do not implement rigid global synteny alignment as the only method. Order and orientation scores are reversal-invariant (section 8).

### 7.9 Multiple candidates

Return all candidates passing thresholds, ranked by composite score (ties: higher reference coverage, then lower contig start). Do not return only the top one. Candidates on separate contigs remain separate.

- **[PROPOSED]** Full ranking key: score descending, reference coverage descending, then contig order in the genome, start, end and event IDs, so the order is always deterministic.
- **[PROPOSED]** Scores are compared at 12-decimal precision when ranking, so floating-point noise never decides a tie that belongs to reference coverage. The reported score is not rounded.
- **[PROPOSED]** Overlap deduplication (7.3) is greedy in rank order: a candidate is dropped when its unpadded span overlaps an already kept candidate on the same contig by more than `max_candidate_overlap` times the shorter span.
- **[PROPOSED]** The checks run in `detection/evaluation.py`: events, windows, evaluation, filtering, ranking, deduplication, per contig. `min_matches` greater than the number of families raises an error before any contig is processed.

────────

## 8. Candidate scoring [DECIDED: simple weighted sum]

Composite score in `[0, 1]`:

```text
score = w_content*content + w_similarity*similarity + w_order*order
      + w_orientation*orientation + w_spacing*spacing
```

Default weights **[DECIDED start values, not validated]**:

| Component | Weight |
|-----------|--------|
| content | 0.40 |
| similarity | 0.30 |
| order | 0.15 |
| orientation | 0.10 |
| spacing | 0.05 |

```python
ScoringWeights(content=0.40, similarity=0.30, order=0.15, orientation=0.10, spacing=0.05)
```

Weights must be non-negative and are normalized to sum to 1. Gene content and protein similarity carry the most weight. All components are exposed on the cluster as `ScoreComponents`, never just the single score.

Component definitions **[PROPOSED]**, each in `[0, 1]`:

- **content** = `reference_coverage`: represented families / total families (families with at least one searchable gene).
- **similarity** = mean over represented families of `identity/100 × coverage`, using the best member and best copy (identity from the aggregated `GeneMatch`, coverage from 7.7).
- **order**: take represented families in reference order and the best copy of each in target order. order = longest monotonic subsequence (increasing or decreasing, whichever is longer) divided by number of represented families. A fully reversed region scores 1.
  - **[PROPOSED]** Reference order is the genomic order of the best members. The subsequence is non-strict, so two families whose best copy is the same (fused) target gene can both count.
- **[PROPOSED]** With fewer than two represented families, order, orientation and spacing are 0, and similarity is 0 when none is represented. Such candidates always fail `min_matches`.
- **orientation**: for each pair of neighbouring represented families in reference order, check whether their relative strand (same or opposite) in the target equals that in the reference. orientation = fraction of agreeing pairs. Reversal-invariant.
- **spacing**: `min(a, b) / max(a, b)`, where `a` is the unpadded candidate span and `b` is the reference span from the first to the last represented family. It penalizes large insertions or compression.

How missing genes, insertions, duplications, splits and fusions affect the score:

- Missing genes lower content (and similarity indirectly).
- Insertions lower spacing only.
- Duplicate copies are not penalized.
- Split and fusion events do not alter the score in v0.1. They are reported as flags. A split reference gene still counts as represented through combined coverage (7.7).

**[PROPOSED] Direction.** A candidate is *reversed* relative to the reference when, taking the represented families' best copies in target order, the longest decreasing run of reference positions is longer than the longest increasing run (the runs of the order component). On a tie, it is reversed when more best copies lie on the opposite strand to their best member than on the same strand. Otherwise it is *forward*. Direction never changes the score; order and orientation stay reversal-invariant.

### 8.1 Coverage metrics

Kept separate:

- `reference_coverage`: represented families / total families.
- `mean_protein_coverage`: mean 7.7 coverage across represented families.
- Reference genes without a protein are excluded from the denominator and reported in the run metadata as `reference_genes_without_protein`.

### 8.2 Missing and additional genes

- Missing reference genes are allowed and reported (as missing families, listing member genes).
- Additional target genes between matched genes are retained.
- Neither automatically disqualifies a candidate beyond the detection rules.

────────

## 9. Split proteins, fused proteins and duplications

### 9.1 Duplications

Retain all homologous copies. Never discard a copy because another scored better. Preserve all target genes and match relationships.

### 9.2 Potential splits

Detect and flag only. A potential split: multiple target proteins align to complementary portions of one reference protein.

- Compute the union of aligned reference positions across the target proteins.
- Count that union toward 7.7 coverage.
- Preserve individual target proteins and alignments.

```text
Reference protein B:  |------------------------|
Target protein B1:    |----------|
Target protein B2:               |-------------|
```

### 9.3 Potential fusions

Detect and flag only. A potential fusion: one target protein aligns to different portions of reference proteins from **different families** (7.6).

- Coverage is calculated independently per reference protein.
- Matching several members of the same family is not a fusion.
- Be careful to distinguish fusions from proteins sharing a common domain.

### 9.4 Heuristics [PROPOSED]

Thresholds need validation on real data. The user chose the starting values below as unvalidated defaults; all are configurable in `RegionDetection` (`min_split_contribution`, `max_fragment_overlap`, `min_fusion_part_coverage`, each 0.2).

- **Split flag**: at least 2 target genes in the candidate each contribute at least 20 % of the reference protein's length not already covered by the others (unique contribution), their pairwise overlap is at most 20 % of the shorter alignment, and the union reaches `min_reference_gene_coverage`.
- **Fusion flag**: one target gene has non-overlapping alignments (overlap at most 20 % of the shorter) to at least 2 different families, each covering at least 20 % of the target protein's length.
- Only GeneMatches from the candidate's evidence genes that pass the pair-level thresholds are used. Several members of one family never make a fusion.
- Fragments and fused parts are chosen greedily: largest first (unique contribution for splits, part length for fusions; ties in genomic order or protein position), skipping one that overlaps an already chosen one by more than the limit. A real fusion therefore stays visible when two other families only share a domain.
- Neighbouring orientation and adjacency of fragments are recorded but not required in v0.1. **[DECIDED]** A split event records its fragments in genomic order, `adjacent` (no other gene between them on the contig) and `same_strand`. A fusion event records the target gene and its families in order along the protein.

────────

## 10. Candidate boundaries and padding

- Boundaries are the outermost matching events of the chain or window, not predefined flanking genes. With zero padding, the candidate runs from the start of the first matched gene to the end of the last matched gene.
- All intervening genes, including unrelated insertions, are retained.
- Padding options: `padding_genes` and `padding_bp`, default 0. If both are given, the larger resulting extent is used **[PROPOSED]**. Padding is clipped to the contig ends.
- **[PROPOSED]** `padding_genes` adds whole genes from the contig's sorted gene list on each side of the outermost events; `padding_bp` adds base pairs.
- **[PROPOSED]** `max_region_bp` is checked against the **unpadded** span. Padding is added afterwards and affects the reconstructed/displayed interval only, not the evidence used for detection.
- A candidate's `touches_contig_edge` flag is set when its first or last matched gene is partial or lies within `edge_margin_bp` of the contig end (**[PROPOSED]** default 500). It signals possible truncation by assembly fragmentation.
  - **[PROPOSED]** "Within" is inclusive and applies to both contig ends: a boundary event starting at most `edge_margin_bp` from the contig start, or ending at most that far from the contig end.

────────

## 11. Visualization

Matplotlib only in v0.1.

### 11.1 Overview plot

- Reference region as the top line, each candidate on its own line.
- Gene arrows show strand.
- Genes in the same reference **family** share a colour across lines **[DECIDED]** (paralog members share a colour).
- Missing families are shown in neutral grey. Additional target genes use a neutral colour.
- Candidate score, reference coverage and flags (split, fusion, edge) are shown.
- Contig IDs and 1-based genomic coordinates are shown.
- Gene labels where practical.
- **[PROPOSED]** Every candidate line is drawn in the reference direction, and never silently: a reversed candidate (section 8) keeps its real 1-based coordinates on a descending x-axis, its title ends with `· reversed`, and its axis reads `<contig> position (1-based, reversed)`. The reference line is never flipped. This holds for the overview, the combined overview and the candidate plot.

- **[DECIDED]** Family colours: families take the eight validated colour-blind-safe palette hues in reference order (blue `#2a78d6`, orange `#eb6834`, aqua `#1baf7a`, yellow `#eda100`, magenta `#e87ba4`, green `#008300`, violet `#6250d6`, red `#e34948`); family 9 reuses hue 1, and so on. Colour is never the only cue: every coloured gene carries its family ID as a label, and a legend maps family IDs to colours.
- **[PROPOSED]** A matched target gene takes the colour of the family of its best-bit-score passing match (7.6); on equal bit scores, the family that comes first in reference order wins. Additional genes and reference genes without a protein are light grey (`#d9d8d2`) and unlabelled.
- **[PROPOSED]** A gene arrow is labelled with its family's `/gene` names (unique, in genomic order, joined with `/`), else with the family ID (`anmK`, `F002`); products never go on arrows. Neighbouring arrow labels alternate between two heights so short neighbours do not collide. Legend entries read `F001 anmK` for a named family and `F002` otherwise, with a name longer than 24 characters cut to 23 characters and `…`. The legend wraps into as many columns as fit the figure width.
- **[PROPOSED]** Missing families are not drawn as a note on plot lines; they are reported in `summary.tsv` (`missing_families`, 14) and the API. On the candidate plot they remain visible as grey reference genes.
- **[PROPOSED]** A matched target gene's opacity shows how conserved it is: `0.5 + 0.5 × (identity − min_identity) / (100 − min_identity)`, clamped to [0.5, 1] (1 when `min_identity` is 100), where identity is that of the passing match deciding its colour and `min_identity` is the search's pair-level threshold (6.2). The scale is fixed, so the same identity looks the same in every plot. Reference genes, grey genes and best-copy outlines stay fully opaque. When a plot draws matched target genes, its legend ends with an opacity key of three grey swatches labelled `identity 30%`, `identity 65%` and `identity 100%` (`min_identity`, the midpoint and 100).
- **[PROPOSED]** With several target genomes, a combined overview shows the reference line, then every candidate of every genome: genomes in input order, candidates in rank order within each. Each candidate title starts with its genome name (`a · #1 candidate_001 · …`). Genomes without candidates are listed in a note on the reference line (`no candidates: c`). Family colours, labels and the legend are shared, since every genome uses the same reference.
- **[PROPOSED]** Titles and the no-candidates note of the combined overview use genome labels, and its lines follow the genome order (12), so a `--names` file sets it.
- **[PROPOSED]** Every line has its own 1-based x-axis over its interval: the reference region, or a candidate's display interval (10).
- **[PROPOSED]** The overview draws at most `PlotOptions.candidates` lines per genome (default 3; 0 draws all), the best-ranked ones. The reference line's note says what was left out: `2 more candidates in summary.tsv` for one genome, `more candidates in summary.tsv: a +2, b +1` in the combined overview. Tables and candidate plots always cover every candidate.
- **[PROPOSED]** An overview longer than `PlotOptions.lines_per_page` lines (default 30, at least 2) is split into pages. Each page starts with the reference line and carries the shared legend; with more than one page, the reference note ends with `page i/n`.

### 11.2 Individual candidate plot

Detailed plot for a single candidate, showing every copy and the score components.

- **[DECIDED]** A panel below the gene tracks shows one horizontal bar per score component on a shared 0–1 axis with its value printed; the weighted score is in the title.
- **[PROPOSED]** The reference track shows genes of missing families in mid grey (`#9a998f`) with their labels kept. The best copy of each represented family has a thicker outline.

### 11.3 Responsibilities

Plotting consumes `Cluster` and result models. It never reruns searches or detection.

### 11.4 API and files [PROPOSED]

- `results.plot()` returns the overview and `results.plot_candidate(rank)` a candidate plot, each as a `matplotlib.figure.Figure` that is not shown.
- `results.write()` adds `plots/regions.png` and `plots/candidate_NNN.png` (PNG, 150 dpi, no Software metadata, so files stay byte-identical). Stale `candidate_*.png` files are removed like the candidate TSVs. With no candidates, `regions.png` shows the reference line and the note "no candidates".
- **[PROPOSED]** `PlotOptions(candidates: int = 3, lines_per_page: int = 30)`, a frozen dataclass in `config.py`, is validated at construction. `write(output_dir, plots=PlotOptions())` takes it, and `plots=None` writes no plots. Page 1 is `plots/regions.png`, later pages `plots/regions_002.png`, `plots/regions_003.png`, …; stale `regions_*.png` are removed like stale candidate files. `plot()` returns page 1 and `plot_pages()` every page, for `RegionResults` and `GenomeSetResults`.
- **[PROPOSED]** With `plots=None`, plot files an earlier run wrote (`regions.png`, `regions_*.png`, `candidate_*.png`) are removed, so the directory matches the run (14), and no `plots/` directory is created.

────────

## 12. Python API

```python
from homolocus import Genome, ReferenceRegion, RegionFinder

reference = ReferenceRegion.from_genbank(
    "ecoli.gb", contig="NC_000913.3", start=1_900_000, end=1_920_000,
)
genome = Genome.from_file("target.fna", annotation="auto")

finder = RegionFinder(reference=reference)       # accepts search, detection, family, scoring configs
results = finder.find(genome)

results.candidates
results.best_candidate
results.to_dataframe()
results.plot()
```

The API exposes intermediate models (alignments, matches, families, candidates, score components). The exact signatures stay provisional. No save/load of results in v0.1.

- **[DECIDED]** `results.to_dataframe()` has one row per candidate (the `summary.tsv` columns); `results.genes_dataframe()` has the gene-level rows (the candidate TSV columns plus `candidate_id`).
- **[PROPOSED]** `finder.find_all(genomes)` returns `GenomeSetResults`: the `RegionResults` of each genome in input order, iterable and indexable by genome name. Its `to_dataframe()` holds every candidate of every genome with a leading `genome` column, `plot()` returns the combined overview (11), and `write(output_dir)` writes the multi-genome layout (14). Like `find`, it has no file side effects.
- **[PROPOSED]** `read_genome_names(path)` reads a mapping file (13) into file name → presentation name, in file order. `apply_genome_names(genomes, names)` orders genomes by it, listed genomes first in file order and the others after in input order, sets their `display_name`, and warns about each line that matches no genome. A genome matches a line by its source file name, ignoring the path; a compressed genome is listed as named on disk (`strain.fna.gz`). `GenomeSetResults` keeps input order, so a run follows the mapping.
- **[DECIDED]** `finder.find(genome)` has no file side effects: DIAMOND runs in a temporary directory and the raw tables and token mapping are kept in the results. `results.write(output_dir)` writes every output file of section 14, the raw DIAMOND tables byte for byte.
- **[PROPOSED]** `RegionFinder(reference, search=DiamondSearch(), family_grouping=FamilyGrouping(), detection=RegionDetection(), scoring=ScoringWeights())`. `results.best_candidate` is `None` when nothing passes.

────────

## 13. CLI

A thin wrapper over the Python API with no duplicated logic.

```bash
homolocus find \
    --reference ecoli.gb --contig NC_000913.3 --start 1900000 --end 1920000 \
    --genome target.fna --output results/
```

```bash
homolocus find --reference ecoli_pai.gbk --genome target.gbk --output results/
```

`--start` and `--end` are 1-based inclusive. The CLI exposes DIAMOND thresholds, `max_target_seqs`, detection thresholds, family grouping, scoring weights and padding. A configuration file is not required in v0.1.

- **[PROPOSED]** Options are generated from the configuration dataclasses: one flag per field of `DiamondSearch`, `FamilyGrouping` and `RegionDetection`, named after the field (`min_identity` → `--min-identity`), and one per `ScoringWeights` field prefixed `--weight-` (`--weight-order`). `--no-family-grouping` disables grouping (`family_min_identity=None`) and `--annotation {auto,pyrodigal}` selects target annotation. Unset flags keep the dataclass defaults, which `--help` shows. `--start` and `--end` are given together or not at all.
- **[PROPOSED]** `--genome` takes one or more paths (`--genome a.fna b.fna`, or repeated). All genomes are read before any search, so an unreadable one exits with 1 and no output. With several genomes each stdout candidate line starts with its genome name, and a genome without candidates prints `<genome>: no candidates`.
- **[PROPOSED]** `--names FILE` gives genomes presentation names and an order. The file has one genome per line: the genome's file name (path ignored), a tab, and its presentation name; blank lines and lines starting with `#` are ignored and there is no header. Genomes follow the file's order, unlisted ones after it in command-line order. A line that matches no `--genome` file is a warning. A malformed line, or a file name or presentation name listed twice, is an input error (exit 1). Stdout lines use labels.
- **[PROPOSED]** Exit codes: 0 on success; 2 for a usage error or an invalid setting; 1 when input cannot be read or is invalid, or DIAMOND is missing or fails; 130 when the run is interrupted (Ctrl-C), after printing `homolocus: interrupted`. Errors are one line on stderr without a traceback.
- **[PROPOSED]** After writing, `find` prints one line per candidate to stdout (rank, candidate ID, `contig:start-end` 1-based, score, reference coverage, flags), or `no candidates`, then the output directory. Warnings still go to stderr through the `warnings` module.
- **[PROPOSED]** The package reports progress through the `logging` logger `homolocus` and never adds handlers itself, so API users see it only when they configure logging (`logging.basicConfig(level=logging.INFO)`). `INFO` reports the reference (its genes, and its families after the self-search), each genome read, `searching genome i/n: <label>` followed by `<label>: k candidates`, and the output directory written. `DEBUG` adds each DIAMOND command line.
- **[PROPOSED]** `find` prints progress to stderr as `homolocus: <message>`, at `INFO` by default. `-q/--quiet` shows only warnings and errors; `-v/--verbose` adds `DEBUG`. Warnings still go through the `warnings` module, and stdout is unchanged.
- **[PROPOSED]** `homolocus --version` prints `homolocus <version> (DIAMOND <version>)`, or `(DIAMOND not found)` when no `diamond` is on `PATH`, and exits 0.
- **[PROPOSED]** A `plots` option group is generated from `PlotOptions`, prefixed `--plot-` (`--plot-candidates`, `--plot-lines-per-page`). `--no-plots` writes tables and `analysis.json` only.

────────

## 14. Outputs

```text
results/
├── summary.tsv
├── candidates/
│   ├── candidate_001.tsv
│   └── ...
├── plots/
│   ├── regions.png
│   ├── candidate_001.png
│   └── ...
├── alignments/
│   ├── diamond.tsv             # raw main search
│   ├── reference_self.tsv      # raw self search (families)
│   └── protein_ids.tsv         # DIAMOND token -> canonical gene ID
└── analysis.json
```

- All coordinates in files are 1-based inclusive.
- `analysis.json` holds run metadata only: parameters, tool versions, thresholds, warnings, skipped features (pseudogenes, compound CDS, untranslatable genes), `reference_genes_without_protein`, and the reference families. It is **not** a reloadable results object.
- **[PROPOSED]** `summary.tsv` columns: `rank, candidate_id, contig, start, end, span_bp, direction (forward|reversed), score, content, similarity, order, orientation, spacing, reference_coverage, mean_protein_coverage, n_target_genes, n_matched_genes, n_missing_families, missing_families, flags`. `missing_families` lists the missing family IDs joined with `;` (`F002;F007`), empty when none.
- **[PROPOSED]** Candidate TSV columns: `target_gene_id, locus_tag, product, start, end, strand, role (matched|additional), reference_gene_id, family_id, family_name, identity, query_coverage, subject_coverage, bit_score, is_best_copy, flags`.
- **[PROPOSED]** Writing results is idempotent. Re-running into an existing output directory overwrites the files the package manages and leaves other files untouched. Identical input and parameters give identical output files: rows and candidates are ordered deterministically, and no output file carries a timestamp except the run metadata in `analysis.json`.
- **[DECIDED]** `analysis.json` carries no timestamp either, so identical input and parameters give byte-identical output directories.
- **[PROPOSED]** Output details:
  - Candidate IDs are `candidate_001`, `candidate_002`, … in rank order.
  - `summary.tsv` `start`/`end` are the unpadded span; `display_start`/`display_end` follow them with the padded interval (10).
  - The candidate TSV has one row per gene and GeneMatch (genes in genomic order, then reference gene ID), and one row with empty match columns for a gene without matches.
  - Flag tokens are joined with `;` and empty when none: `split`, `fusion`, `partial`, `edge` in the summary; `partial`, `split_fragment`, `fused` on gene rows.
  - Floats are written with 6 significant digits.
  - `analysis.json` is sorted, indented JSON with `homolocus_version`, `tool_versions` (DIAMOND, Pyrodigal, Biopython), `parameters` (including `plots`: the `PlotOptions` values, or `null` without plots), `reference`, `target`, `families` (each with `id`, `name` and `gene_ids`), `reference_genes_without_protein`, `skipped_features`, `untranslatable_genes` (coordinates 1-based) and `warnings`. `warnings` lists the reference's parse warnings prefixed `reference: `, then the target's prefixed `target: `, then the warnings raised during the run, each group in the order raised.
  - Managed files are `summary.tsv`, `candidates/candidate_*.tsv`, `plots/regions.png`, `plots/regions_*.png`, `plots/candidate_*.png`, `alignments/*` and `analysis.json`. Rewriting removes stale `candidate_*.tsv` and `candidate_*.png` files from an earlier run and leaves other files alone.
- **[PROPOSED]** With several target genomes, `results/<genome>/` holds the layout above for each genome, `results/summary.tsv` lists every candidate of every genome with a leading `genome` column (genomes in input order, then rank), and `results/plots/regions.png` is the combined overview (11). A single genome keeps the flat layout above. Subdirectories of genomes from an earlier run are left alone.
- **[PROPOSED]** The combined summary's `genome` column holds each genome's label; folders are named after `name`. `analysis.json`'s `target` gains `display_name` (`null` when none).
  - With family grouping disabled the self-search is not run and `alignments/reference_self.tsv` is empty.

────────

## 15. Dependencies and packaging

```toml
[build-system]
requires = ["hatchling>=1.27"]
build-backend = "hatchling.build"

[project]
name = "homolocus"
version = "0.1.0"
description = "Find homologous genomic regions in bacterial genomes from the proteins of a reference region, tolerant of gene loss, insertions, duplications and rearrangements."
requires-python = ">=3.10"
dependencies = ["biopython>=1.80", "pandas>=2.0", "matplotlib>=3.8", "pyrodigal>=3.0"]

[project.optional-dependencies]
dev = ["pytest", "pytest-cov", "hypothesis", "ruff", "mypy", "pandas-stubs"]

[project.scripts]
homolocus = "homolocus.cli:main"
```

- The block above is an excerpt. `pyproject.toml` is authoritative and also sets the readme, the MIT license (`license`, `license-files`), the wheel package path and the ruff, mypy, pytest and coverage configuration. `hatchling>=1.27` is required for the `license-files` metadata.
- Mandatory: Biopython, pandas, Matplotlib, Pyrodigal.
- DIAMOND is an external executable checked at runtime.
- Tests needing DIAMOND are marked `@pytest.mark.diamond` and skipped with a clear message when it is absent.
- **[DECIDED]** Continuous integration (GitHub Actions, `.github/workflows/tests.yml`) runs `ruff check`, `ruff format --check`, `mypy src` and the full test suite on Python 3.10–3.13 for every pull request and push to `main`. It installs DIAMOND 2.2.9 from the official release, pinned by SHA-256. It sets `HOMOLOCUS_REQUIRE_DIAMOND=1`, which makes a missing DIAMOND a test failure instead of a skip.
- **[PROPOSED]** The full suite keeps at least 98 % branch coverage: `fail_under` in the coverage configuration, checked in every CI test job.
- **[PROPOSED]** A scheduled workflow (`.github/workflows/prerelease.yml`) runs the full suite weekly, and on demand, on Python 3.13 against the newest pre-releases of the dependencies (`pip install --pre`), to catch upstream changes early. It does not block pull requests.
- **[PROPOSED]** Runtime dependencies carry tested lower bounds and no upper bounds: Biopython 1.80, pandas 2.0, Matplotlib 3.8, Pyrodigal 3.0. A CI job (`minimum versions`) installs exactly these versions on Python 3.10, with NumPy below 2 because pandas 2.0 is built against NumPy 1, and runs the full suite. A floor that cannot be proven is raised, not guessed.
- **[PROPOSED]** Releases: a CI job (`build`) builds the sdist and wheel on every pull request, checks them with `twine check --strict`, checks the wheel contents and installs the wheel into a clean environment to run `homolocus --help`. `.github/workflows/release.yml` builds again on a `v*` tag, refuses a tag that differs from the package version, and publishes to PyPI through Trusted Publishing from the `pypi` environment. Tagging and enabling the publisher are the maintainer's steps. `CHANGELOG.md` follows Keep a Changelog.
- mypy runs in strict mode. Missing type information is ignored only for Biopython, which ships none. pandas is typed through `pandas-stubs`.
- The quality-gate tools of section 17 ([crapper](https://github.com/unclebob/crapper), [dryer](https://github.com/unclebob/dryer), [mutator](https://github.com/unclebob/mutator)) are external, optional development tools checked out next to the repository. They are not Python dependencies of the package. Their output (`.metrics/`, `target/`) is not committed.

────────

## 16. Package structure

```text
homolocus/
├── pyproject.toml
├── README.md
├── LICENSE
├── SPECIFICATION.md
├── CLAUDE.md
├── CHANGELOG.md
├── .github/workflows/
│   ├── tests.yml                  # lint, types, tests, minimum versions, build
│   └── release.yml                # build and publish on a v* tag
├── src/homolocus/
│   ├── __init__.py
│   ├── cli.py
│   ├── finder.py
│   ├── results.py                 # RegionResults, summary and gene rows, data frames
│   ├── genome_names.py            # presentation names and order from a mapping file (13)
│   ├── output.py                  # writes the output directory (14), I/O layer
│   ├── config.py                  # DiamondSearch, RegionDetection, FamilyGrouping, ScoringWeights
│   ├── coordinates.py             # 0-based half-open <-> 1-based inclusive
│   ├── intervals.py               # interval unions for coverage (5.7, 7.7)
│   ├── models/
│   │   ├── __init__.py
│   │   ├── gene.py
│   │   ├── contig.py
│   │   ├── genome.py
│   │   ├── reference.py
│   │   ├── skipped_feature.py
│   │   ├── gene_family.py
│   │   ├── alignment.py
│   │   ├── gene_match.py
│   │   └── cluster.py
│   ├── annotation/
│   │   ├── __init__.py
│   │   ├── genbank.py
│   │   ├── genome_file.py         # format detection, Genome.from_file dispatch
│   │   ├── prodigal.py
│   │   └── translation.py         # CDS translation rules (3.2), pure
│   ├── search/
│   │   ├── __init__.py
│   │   ├── diamond.py
│   │   └── protein_index.py       # DIAMOND tokens <-> gene IDs, FASTA, mapping file
│   ├── detection/
│   │   ├── __init__.py
│   │   ├── families.py
│   │   ├── matching.py
│   │   ├── candidate_generation.py
│   │   ├── scoring.py
│   │   ├── flags.py               # split / fusion / edge flags
│   │   └── evaluation.py          # coverage, Cluster building, filter, rank, dedup
│   └── visualization/
│       ├── __init__.py
│       ├── layout.py              # tracks and gene glyphs from results (pure)
│       └── region_plot.py
└── tests/
    ├── test_models/
    ├── test_annotation/
    ├── test_search/
    ├── test_detection/
    └── test_visualization/
```

────────

## 17. Testing strategy

Cover biological edge cases, not only successful parsing. Use small synthetic genomes and proteins for deterministic unit tests, plus an end-to-end test with a small known reference and target.

**Models and coordinates:** 0-based half-open internally, 1-based inclusive conversion at I/O (round-trip), forward and reverse genes, sorting, empty contigs, coding-orientation sequences, gene ID format.

**GenBank parsing:** with and without `/translation`, `/transl_table`, `/codon_start`, reverse-strand translation, pseudogene exclusion, internal stop handling, compound CDS skipped with warning, boundary-overlapping CDS excluded, contained CDS included, partial genes flagged.

**FASTA and Pyrodigal:** one contig, multiple contigs, contigs without genes, short contigs in meta mode, coordinates and proteins, `auto` selection logic.

**DIAMOND and matching:** multiple hits per query, multiple HSPs per pair, coordinate conversion, union-based coverage without double counting, aggregated evidence, no per-HSP coverage filtering, `max_target_seqs` behaviour.

**Families:** paralog grouping at the identity/coverage thresholds, shared-domain genes not grouped, grouping disabled, family counted once in coverage.

**Region detection:** missing genes, inserted genes, duplicated homologues as separate events, rearranged order, reversed candidate, distance and span limits, long chains split into windows and deduplicated, multiple candidates on one contig, separate contigs never combined, `min_matches` and `min_reference_coverage`, padding after span check, contig-edge flag.

**Scoring:** each component on constructed cases, reversal invariance of order and orientation, weight normalization, tie-breaking.

**Split/fusion:** complementary fragments counted as one represented gene, one target matching distinct families, same-family matches not flagged as fusion, overlapping HSPs not double counted, flags carried into output.

**Visualization:** multiple candidate lines, consistent family colours, arrows and strand, missing and additional genes distinguishable.

**Test behaviour, not implementation.** Tests call the public API and module functions the way their user does and check the result that user would see. A test that would still pass if every function returned nothing is removed. Write the failing test before the code it covers.

**Property-based tests.** Use `hypothesis` for invariants such as coordinate round-trips (1-based inclusive ↔ 0-based half-open) and interval unions (coverage never double counts and never exceeds 1).

**Quality gates.** Run on `src/` before a change is considered done:

- Every function has a CRAP score below 10 (`crapper`). Reduce it by splitting the function, testing it, or both. Ruff's complexity limit of 9 is an early static guard for the same rule.
- No pair of functions has a duplication score of 0.8 or above (`dryer`).
- No surviving mutants in changed code (`mutator`): cover the uncovered, kill the survivors, delete dead code.

────────

## 18. Open questions and proposed defaults awaiting confirmation

**Closed in revision 2:** duplicate handling, coverage levels, user coordinates, candidate grouping, gene IDs, circular/compound CDS, `auto` annotation, translation failures, paralog families, scoring approach.

**Proposed defaults to confirm (implement as written, revisit in review):**

1. Matching-event rule (OR of query and subject coverage at 0.5), 7.2.
2. Family thresholds (90 % identity, 0.8 coverage, connected components), 7.6.
3. `min_matches` AND `min_reference_coverage`, with `min_matches >= 2`, 7.5.
4. Candidate overlap suppression (`max_candidate_overlap=0.5`), 7.3.
5. Padding applied after the span check and larger-of-two combination, section 10.
6. Component formulas for similarity, order, orientation and spacing, section 8.
7. DIAMOND fields and `--max-hsps 0`, 6.2 and 6.4.
8. Units: identity in percent, coverage as fractions, 0.2.
9. Output columns, section 14.
10. Pyrodigal training strategy for assemblies of 20 kb or more, 2.2.
11. `edge_margin_bp` default of 500, section 10.
12. Idempotent, deterministic output writing, section 14.
13. Invalid supplied `/translation` definition and one summary warning per contig, 3.2.
14. Pyrodigal meta mode restricted to table-11 models, and training only when a contig is at least 20 kb, 2.2.
15. Inclusive pair-level threshold comparisons, 5.7.
16. DIAMOND token prefixes `r`/`t`, minimum DIAMOND 2.0, and deterministic sorting of parsed alignments, 4.1, 6.2, 6.4.
17. Family linking from aggregated self-search evidence with inclusive thresholds, and family ID numbering, 7.6.
18. Literal intervening-gene count, single over-long event windows, window ordering, and gene-based padding, 7.3, 10.
19. Best member and best copy, non-strict order, components with fewer than two families, full ranking key, greedy overlap deduplication, inclusive edge margin, 7.7, 7.9, 8, 10.
20. Split and fusion thresholds (20 % each, unvalidated) and greedy choice of fragments and parts, 9.4.
21. 12-decimal score precision when ranking, 7.9.
22. RegionFinder constructor, candidate IDs, summary and gene-row layout, flag tokens, float precision, analysis.json keys, managed files, no self-search without grouping, 12, 14.
23. Matched-gene colour and its tie rule, grey shades, per-line axes, best-copy outline, plotting API and plot files, 11.
24. CLI options generated from the configuration dataclasses, exit codes and the stdout summary, 13.
25. Tested dependency floors without caps, the build and minimum-versions CI jobs, and the tag-triggered PyPI release, 15.
26. Family names from `/gene`, else `/product`, joined with `/`; `F001 sitA` labels truncated at 24 characters in plots; `family_name` in candidate TSVs, 5.5, 7.6, 11, 14.
27. Several target genomes per run: shared families, per-genome searches, genome names from file names, `find_all` and `GenomeSetResults`, repeatable `--genome`, per-genome subdirectories with a combined summary and overview plot, 2.2, 5.3, 11, 12, 13, 14.
28. Candidate direction (order runs, then strands), `Cluster.reversed`, the `direction` summary column, and reversed candidates drawn on a descending axis, 5.8, 8, 11, 14.
29. Uninformative products do not name families; arrows labelled by `/gene` name else family ID, staggered in two rows; a wrapping legend; missing families reported in `summary.tsv` instead of on plots, 7.6, 11, 14.
30. Matched-gene opacity from the identity of its colour match on a fixed scale (50 % at `min_identity`, 100 % at 100 % identity) with an opacity key, 11.
31. Genome presentation names and order from a `--names` mapping file (file name, tab, name; path ignored), shown as labels while folders keep `name`; unmatched lines warn, malformed or duplicated ones exit 1, 5.3, 11, 12, 13, 14, 16.
32. gzip-compressed reference and target files, detected from content, with `.gz` removed from genome names, 2.1, 2.2, 5.3, 12.
33. Short model reprs without sequences or raw output; parse warnings kept on `Genome` and `ReferenceRegion` and written to `analysis.json` with `reference: `/`target: ` prefixes; family names in `analysis.json`, 5, 5.3, 5.4, 14.
34. Progress through the `homolocus` logger with no handlers in the package; CLI progress on stderr at INFO with `-q`/`-v`; `--version` with the DIAMOND version; exit code 130 on interrupt, 13.
35. Overview capped at `PlotOptions.candidates` per genome with a note naming what was left out; pages of at most `lines_per_page` lines (`regions.png`, `regions_002.png`, …); `plot_pages()`; `--no-plots`; `parameters.plots` in `analysis.json`, 11, 11.4, 13, 14.
36. A 98 % branch-coverage floor and a weekly pre-release CI job, 15.

**Still open:**

- Distinguishing true fusions from shared domains beyond the family rule (the 9.4 overlap limit covers the simple case).
- Weights and component definitions need validation on real data (v0.1 ships them as unvalidated defaults).

Do not silently resolve these during implementation. Record the chosen default and ask when behaviour is materially affected.

────────

## 19. Implementation principles

- Keep parsing, annotation, search, matching, candidate generation, scoring and plotting in separate modules.
- Keep DIAMOND-specific logic isolated.
- Preserve raw evidence wherever practical.
- Avoid irreversible one-to-one matching assumptions.
- Keep all significant homologous copies and individual HSPs.
- Keep protein-level confidence separate from region-level confidence.
- Treat gene order as evidence, not a requirement.
- Analyse contigs independently.
- Make thresholds configurable and document provisional defaults.
- Convert coordinates only at I/O boundaries.
- Keep the Python API primary. The CLI is a wrapper.
- Add tests alongside each component.
- Do not implement features excluded from v0.1.
- Parse and validate outside data once, where it enters (GenBank, FASTA, DIAMOND output, CLI arguments). Inside the package, trust it.
- Model the domain so illegal states cannot be built: strand is +1 or −1, `start < end`, coverage lies in `[0, 1]`, identity in `[0, 100]`. Enforce this when objects are constructed.
- Never silence the type checker (no blanket `# type: ignore`, no `cast` to hide an error).
- Make operations idempotent (section 14).
- Use the terms of section 20 in code, tests and documentation.

────────

## 20. Glossary

Short definitions. The rules live in the referenced sections.

| Term | Meaning |
|------|---------|
| Reference region | The interval on one reference contig whose proteins define what is searched for (2.1, 5.4). |
| Target genome | The single assembly searched for the region. It has one or more contigs (2.2, 5.3). |
| Contig | One sequence of an assembly, treated as linear and analysed independently (2.3, 5.2). |
| Gene | A protein-coding CDS. No other feature types are modelled in v0.1 (3.1, 5.1). |
| Canonical gene ID | `{contig_id}\|{start_1based}-{end_1based}\|{strand_char}`, the only key for genes (4.1). |
| Partial gene | A gene truncated at a contig edge or marked partial. Kept and flagged (2.2). |
| Pseudogene | A CDS annotated as a pseudogene or with an internal stop. Excluded entirely (3.1, 3.2). |
| Untranslatable gene | A gene whose protein could not be obtained. Kept but not searchable (3.2, 8.1). |
| Skipped feature | A CDS excluded during parsing (pseudogene, internal stop or compound location), recorded on the model (3.2, 5.10). |
| HSP / Alignment | One DIAMOND high-scoring segment pair between a reference and a target protein (5.6). |
| GeneMatch | All HSPs between one reference protein and one target protein, aggregated (5.7). |
| Pair-level thresholds | E-value and identity filters applied to a GeneMatch. No coverage filter at this level (5.7, 6.2). |
| Gene family | Reference genes grouped as paralogs by the self-search. Counted once in coverage and scoring (5.5, 7.6). |
| Matching event | A target gene with a GeneMatch that passes pair-level thresholds and the event coverage rule (7.2). |
| Intervening genes | Genes between two consecutive matching events in the contig's sorted gene list (7.3). |
| Chain | Matching events linked by single linkage under the distance limits (7.3). |
| Window | A part of a chain longer than `max_region_bp`, anchored on one event (7.3). |
| Span | Distance from the start of the first matching event to the end of the last, without padding (7.3, 10). |
| Candidate | One evaluated chain or window on one contig. The model is `Cluster` (5.8, 7.3). |
| Matched gene | A target gene in a candidate that is a matching event (5.8, 7.7). Genes with only weak, sub-event matches are additional. |
| Additional gene | A target gene in a candidate that is not matched (5.8, 8.2). |
| Reference-gene coverage | Union of aligned reference positions across a candidate's target genes, as a fraction of the reference protein (7.7). |
| Represented family | A family whose best member reaches `min_reference_gene_coverage` in a candidate (7.7). |
| Missing family | A family that is not represented in a candidate (8.2). |
| Reference coverage | Represented families divided by families with at least one searchable gene (8.1). |
| Duplicated homologue | One of several target copies matching the same family. Each copy is its own matching event (7.4, 9.1). |
| Best copy | The copy of a family in a candidate with the highest bit score. Used for order and orientation (7.4). |
| Potential split | Several target proteins covering complementary parts of one reference protein. Flagged only (9.2, 9.4). |
| Potential fusion | One target protein covering parts of reference proteins from different families. Flagged only (9.3, 9.4). |
| Padding | Genes or base pairs added around a candidate for display, after the span check (10). |
| Contig edge | A candidate touches it when its first or last matched gene is partial or within `edge_margin_bp` of a contig end (10). |

End of v0.1 working specification (revision 2)
