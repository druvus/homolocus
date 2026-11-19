# homolocus

Find homologous genomic regions in bacterial genomes from the proteins of a reference region, tolerant of gene loss, insertions, duplications and rearrangements.

You give homolocus a reference region, such as a pathogenicity island in a GenBank genome, and a target assembly. It searches the reference proteins against the target proteins with [DIAMOND](https://github.com/bbuchfink/diamond), then finds and ranks the places in the target that hold a related collection of those proteins. It doesn't need predefined flanking genes. Protein content is the main evidence; gene order and orientation support it but are not required.

It tolerates:

- missing reference genes and inserted unrelated genes;
- duplicated homologues (every copy is kept);
- rearranged gene order and regions in reverse orientation;
- split and fused proteins (detected and flagged);
- repeated reference families such as IS-element transposases (paralogs are grouped into one family).

> **Status: v0.1, pre-alpha.** The thresholds and score weights are starting values that haven't been validated on real data. They are all configurable, and the [specification](SPECIFICATION.md) lists those still awaiting review (section 18).

## Install

homolocus needs Python 3.10 or later and the DIAMOND executable.

1. Install DIAMOND 2.0 or later from its [releases](https://github.com/bbuchfink/diamond/releases), and put `diamond` on your `PATH`. CI uses DIAMOND 2.2.9. You can check with `diamond version`.
2. Install homolocus from a checkout of this repository (it isn't on PyPI yet):

```bash
git clone https://github.com/emilhaegglund/homolocus.git
cd homolocus
pip install .
```

Biopython, pandas, Matplotlib and Pyrodigal are installed with it.

## Quick start: command line

Search for a region of a reference genome, given in 1-based inclusive coordinates:

```bash
homolocus find --reference ecoli.gbk --contig NC_000913.3 --start 1900000 --end 1920000 \
    --genome target.fna --output results/
```

Or give a GenBank file that holds only the region:

```bash
homolocus find --reference ecoli_pai.gbk --genome target.gbk --output results/
```

`find` writes the output directory and prints one line per candidate, then where the results are:

```text
#1 candidate_001  contig_7:3061-4548  score 0.900  coverage 0.75
#2 candidate_002  contig_2:88120-101400  score 0.612  coverage 0.50  split;edge
results written to results/
```

Search several target genomes in one run by giving `--genome` more than one file:

```bash
homolocus find --reference ecoli_pai.gbk --genome strain_a.fna strain_b.gbk --output results/
```

The reference families and its self-search are computed once. Each genome is then searched on its own, so its results are the same as a run with that genome alone. Genomes are named after their files (`strain_a`, `strain_b`), and names must be unique. Each candidate line starts with its genome, and a genome without candidates prints `strain_b: no candidates`.

Run `homolocus find --help` to see every option and its default. The exit code is:

- 0 on success, including when nothing is found;
- 2 for a usage error or an invalid setting;
- 1 when an input can't be read or is invalid, or DIAMOND is missing or fails.

Errors are one line on stderr.

## Quick start: Python

```python
from homolocus import Genome, ReferenceRegion, RegionDetection, RegionFinder

reference = ReferenceRegion.from_genbank("island.gbk")  # the whole record is the region
genome = Genome.from_file("target.fna")  # FASTA: genes predicted with Pyrodigal

finder = RegionFinder(reference, detection=RegionDetection(padding_genes=2))
results = finder.find(genome)  # writes no files

best = results.best_candidate  # None when nothing passes
summary = results.to_dataframe()  # one row per candidate, like summary.tsv
genes = results.genes_dataframe()  # one row per gene and match, with candidate_id
results.plot().savefig("overview.png")  # reference and every candidate
results.plot_candidate(1).savefig("candidate_001.png")  # one candidate in detail
results.write("results")  # every output file listed below
```

For several genomes, use `find_all`. It returns the results of every genome, in order and by name:

```python
genomes = [Genome.from_file(path) for path in ("strain_a.fna", "strain_b.gbk")]
found = finder.find_all(genomes)  # families and self-search computed once
found["strain_a"].best_candidate
found.to_dataframe()  # every candidate, with a leading genome column
found.plot()  # the reference and every genome's candidates in one figure
found.write("results")  # one folder per genome plus combined files
```

To take a region from a whole genome, name the contig and the 1-based inclusive coordinates:

```python
reference = ReferenceRegion.from_genbank(
    "ecoli.gbk", contig="NC_000913.3", start=1_900_000, end=1_920_000
)
```

`RegionFinder` takes four settings objects: `search=DiamondSearch(...)`, `family_grouping=FamilyGrouping(...)`, `detection=RegionDetection(...)` and `scoring=ScoringWeights(...)`. Their fields are the options listed below.

Each candidate (`results.candidates`, ranked best first) keeps its evidence:

- `components` holds the five score components;
- `families` holds the per-family coverage, best member and best copy;
- `matches` holds every DIAMOND HSP behind each gene pair;
- `role(gene_id)` says whether a gene is `matched` or `additional`;
- the split, fusion and edge flags.

## Inputs

**Reference.** The reference is a GenBank file: either a whole genome with `--contig`, `--start` and `--end`, or a file holding just the region (the record is the region; `--contig` is needed only when the file has several records).

- Only CDS features lying entirely inside the region are used.
- Pseudogenes, CDS with internal stops and compound (`join`) CDS are skipped and listed in `analysis.json`.
- A supplied `/translation` is used when it is valid. Otherwise the CDS is translated with its `/transl_table` and `/codon_start`.

**Target.** The target is one assembly with one or more contigs, as GenBank or FASTA. The format is detected from the file content, not the extension.

- `--annotation auto` (the default) uses the CDS of a GenBank file when it has any, and predicts genes with Pyrodigal otherwise, which includes every FASTA file.
- `--annotation pyrodigal` always predicts genes.
- Contigs are linear and analysed independently. A candidate never spans two contigs.

## How it works

1. **Families.** Reference proteins are searched against each other, and close paralogs are grouped into families `F001`, `F002`, … (at least 90 % identity and 0.8 coverage by default). A family counts once in coverage and scoring. A family is named from its genes' `/gene` qualifiers in the reference GenBank, or their `/product` when none has a `/gene` (paralogs with different names are joined: `tnpA/tnpB`). Products such as "hypothetical protein" or "conserved hypothetical protein" don't name a family.
2. **Search.** Reference proteins are searched against target proteins with DIAMOND. Every hit and every HSP is kept, then aggregated per gene pair (union-based coverage, never double counted).
3. **Matching events.** A target gene is an event when it has a match that passes the E-value and identity thresholds and covers at least half of either protein. That admits both split fragments and fused proteins.
4. **Candidates.** On each contig, nearby events are chained (limited by intervening genes and the base-pair gap). Chains longer than `--max-region-bp` are split into windows. A candidate needs at least `--min-matches` represented families and `--min-reference-coverage`.
5. **Score.** The score is a weighted sum of five components, each between 0 and 1:
    - **content:** represented families out of all families;
    - **similarity:** identity times coverage of the best copies;
    - **order:** longest increasing or decreasing run of families;
    - **orientation:** agreement of neighbouring strands;
    - **spacing:** candidate span versus reference span.

    Reversal doesn't change order or orientation. A candidate is called *reversed* when its families mostly run backwards compared with the reference (by gene order, then by strand on a tie).

6. **Flags.** Candidates are flagged `split` (one reference protein found as several fragments), `fusion` (one target protein matching several families), `partial` (a partial gene) or `edge` (close to a contig end, so possibly truncated). Flags never change the score.

Candidates are ranked by score, and an overlapping, weaker window of the same locus is dropped.

## Outputs

```text
results/
├── summary.tsv               # one row per candidate, ranked
├── candidates/
│   ├── candidate_001.tsv     # genes and matches of each candidate
│   └── ...
├── plots/
│   ├── regions.png           # reference and every candidate, one line each
│   ├── candidate_001.png     # one candidate with its score components
│   └── ...
├── alignments/
│   ├── diamond.tsv           # raw DIAMOND output of the main search
│   ├── reference_self.tsv    # raw DIAMOND output of the family self-search
│   └── protein_ids.tsv       # DIAMOND sequence IDs to gene IDs
└── analysis.json             # parameters, tool versions, families, warnings, skipped features
```

- **Coordinates:** all coordinates in files are 1-based inclusive. Gene IDs are `contig|start-end|strand`, for example `contig_7|3061-3519|+`.
- **`summary.tsv`:**
    - rank and candidate ID, the contig, and `start`/`end` (the matched span);
    - `display_start`/`display_end` (with padding) and `span_bp`;
    - `direction`: `forward`, or `reversed` when the candidate runs opposite to the reference;
    - the score and its five components;
    - `reference_coverage` and `mean_protein_coverage`;
    - gene, match and missing-family counts, and `missing_families` (the missing family IDs, joined with `;`);
    - `flags`.
- **Candidate TSVs:** one row per gene and match (one empty-match row for a gene without matches).
    - The `role` is `matched` or `additional`.
    - `family_name` is the family's name from the reference annotation, empty when it has none.
    - `is_best_copy` marks the copy used for scoring.
    - Gene `flags` are `partial`, `split_fragment` and `fused`.
- **Plots:** every candidate reads in the reference direction. A reversed candidate keeps its real coordinates on a descending axis, with `reversed` in its title and axis label. Families keep one colour across all plots. Each gene arrow is labelled with its family's `/gene` name, or the family ID when there is none (`anmK`, `F002`), and neighbouring labels alternate between two rows. The legend maps every family to its colour and name (`F001 anmK`, names over 24 characters shortened) and wraps to fit. A matched gene's opacity shows how conserved it is: solid at 100 % identity, fading to half opacity at the search's `--min-identity`, with a small key in the legend; exact identities are in the candidate TSVs. Unrelated genes are grey, and the best copy of each family is outlined. Missing families aren't noted on the plots; `summary.tsv` lists them.
- **Several genomes:** `results/<genome>/` holds the layout above for each genome. `results/summary.tsv` lists every candidate of every genome with a leading `genome` column. `results/plots/regions.png` shows the reference and every genome's candidates, titled by genome, and notes genomes without candidates. A single genome keeps the layout above.
- **Deterministic:** the same input and settings give byte-identical files, with no timestamps. Writing again into a directory replaces the files homolocus manages, removes stale candidate files, and leaves other files alone.

## Options

All settings have defaults, so only `--reference`, `--genome` and `--output` are required. Fractions run from 0 to 1, and identities are percentages.

### DIAMOND search

| Option | Default | Meaning |
|---|---|---|
| `--evalue` | `1e-05` | E-value a gene pair must reach |
| `--min-identity` | `30.0` | identity a gene pair must reach, percent |
| `--sensitivity` | `sensitive` | DIAMOND sensitivity mode, `fast` to `ultra-sensitive` |
| `--max-target-seqs` | `100` | targets per query; 0 is unlimited |
| `--max-hsps` | `0` | HSPs per pair; 0 reports all |
| `--hsp-evalue` | `0.001` | per-HSP E-value passed to DIAMOND |
| `--threads` | `unset` | DIAMOND threads; unset leaves DIAMOND's default |

### Family grouping

| Option | Default | Meaning |
|---|---|---|
| `--family-min-identity` | `90.0` | identity linking two reference genes into a family, percent |
| `--family-min-coverage` | `0.8` | coverage both genes must reach |
| `--no-family-grouping` | `off` | treat every reference gene as its own family |

### Region detection and padding

| Option | Default | Meaning |
|---|---|---|
| `--min-matches` | `3` | represented families a candidate needs (at least 2) |
| `--min-reference-coverage` | `0.3` | share of reference families a candidate needs |
| `--min-reference-gene-coverage` | `0.5` | coverage at which a reference gene counts as represented |
| `--min-event-coverage` | `0.5` | query or subject coverage that makes a target gene an event |
| `--max-intervening-genes` | `10` | most genes between two neighbouring events |
| `--max-inter-match-bp-gap` | `20000` | largest gap between two neighbouring events, bp |
| `--max-region-bp` | `50000` | largest candidate span before padding, bp |
| `--max-candidate-overlap` | `0.5` | overlap above which a weaker candidate is dropped; 1 keeps all |
| `--padding-genes` | `0` | genes added on each side of a candidate |
| `--padding-bp` | `0` | base pairs added on each side; the larger padding wins |
| `--edge-margin-bp` | `500` | distance from a contig end that sets the `edge` flag |
| `--min-split-contribution` | `0.2` | share of a reference protein each split fragment must add |
| `--max-fragment-overlap` | `0.2` | largest overlap between split fragments or fused parts |
| `--min-fusion-part-coverage` | `0.2` | share of a target protein each fused part must cover |

### Scoring weights

Weights are normalised to sum to 1.

| Option | Default | Meaning |
|---|---|---|
| `--weight-content` | `0.4` | represented families out of all families |
| `--weight-similarity` | `0.3` | identity and coverage of the best copies |
| `--weight-order` | `0.15` | gene order |
| `--weight-orientation` | `0.1` | relative strands |
| `--weight-spacing` | `0.05` | candidate span versus reference span |

Other options are `--annotation auto|pyrodigal` and, for the reference, `--contig`, `--start` and `--end`.

## Limitations in v0.1

- No GFF input, and no RNA or other non-CDS features.
- Pseudogenes, compound (`join`) CDS and CDS wrapping the origin of a circular contig are skipped.
- Matches on different contigs are never combined.
- Split and fused proteins are flagged, not reconstructed.
- Plots are static PNGs.
- Results can't be saved and reloaded as an object. `analysis.json` holds run metadata only.

## Development

```bash
pip install -e ".[dev]"
pytest -m "not diamond"   # fast tests, no DIAMOND needed
pytest                    # full suite, needs diamond on PATH
ruff check --fix && ruff format
mypy src
```

Tests that need DIAMOND are marked `diamond` and skipped with a message when it is missing. Set `HOMOLOCUS_REQUIRE_DIAMOND=1` to make them fail instead, as CI does.

CI (`.github/workflows/tests.yml`) runs lint, types and the full suite on Python 3.10–3.13 with DIAMOND 2.2.9. The optional quality gates (CRAP score, duplication and mutation testing) are described in `CLAUDE.md` and section 17 of the specification.

The design, every decision and its status, and a glossary are in [SPECIFICATION.md](SPECIFICATION.md).

## Licence

MIT, see [LICENSE](LICENSE).
