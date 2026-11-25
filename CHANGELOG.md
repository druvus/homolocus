# Changelog

All notable changes to homolocus are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Families are named from the reference GenBank annotation: their genes' `/gene` qualifiers, else `/product`. Names appear in plot labels and legends (`F001 sitA`), in `analysis.json` `families`, in a new `family_name` column of the candidate TSVs, and on `GeneFamily.name` and `FamilyEvidence.family_name`.
- Several target genomes in one run:
  - `--genome` takes one or more files, and `RegionFinder.find_all()` returns `GenomeSetResults`;
  - families and the self-search are computed once, and each genome's results equal a run on it alone;
  - output has a folder per genome, a combined `summary.tsv` with a `genome` column, and a combined overview plot.
- Candidates are called `forward` or `reversed` relative to the reference (gene order, then strand). The call is in `Cluster.reversed` and a new `direction` column of `summary.tsv`. Plots draw every candidate in the reference direction: a reversed one sits on a descending axis, marked `reversed`.
- `summary.tsv` gains `missing_families`, listing each candidate's missing family IDs.
- Plots fade each matched gene by the identity of its match: full opacity at 100 %, half at `--min-identity`. An identity key is in the legend.
- `--names FILE` gives genomes presentation names and an order (file name, tab, name), used in the terminal summary, the combined summary, the plots and `analysis.json`. `Genome.display_name` and `Genome.label` carry it in the API, and output folders keep file names.
- README sections "Reading the results" (what the score does and does not say, how to read the components and flags) and "Troubleshooting" (each input error with its cause and fix), and DIAMOND installation through conda.
- Overview plots scale to many genomes: `regions.png` draws the best three candidates per genome (`--plot-candidates`, 0 for all) and notes how many more are in `summary.tsv`; overviews longer than 30 lines (`--plot-lines-per-page`) continue in `regions_002.png`, …. `--no-plots` skips plotting. In the API, `PlotOptions`, `write(output_dir, plots=...)` and `plot_pages()`; `analysis.json` records the plot options.
- Progress messages: `find` reports on stderr the reference and its families, each genome read and searched with its candidate count, and where results are written. `-q/--quiet` hides them and `-v/--verbose` adds each DIAMOND command. In Python they go to the `homolocus` logger, which prints nothing unless logging is configured.
- `homolocus --version` prints the homolocus and DIAMOND versions.
- Reference and target files may be gzip-compressed (`.fna.gz`, `.gbff.gz`), detected from the content. Genome names drop the `.gz` (`strain.fna.gz` is `strain`), and a damaged gzip file is an input error naming the file.

### Changed

- Plots label gene arrows by `/gene` name, else family ID, on two alternating rows. The legend wraps to fit the figure, and missing families are no longer noted on plot lines.
- Products such as "hypothetical protein" no longer name a family.
- Genomes read from files are named after the file (`Genome.name`), so `analysis.json` now records `target.name` instead of `null`.
- Settings that no candidate could meet now fail before any genome is searched: `min_matches` above the number of reference families, and a reference region without searchable protein-coding genes, which has its own message.
- Clearer errors: a GenBank file saved without its sequence names the file and record and says how to fix it; duplicate gene IDs are listed; duplicate genome names say they come from the file names.
- `GenomeSetResults[key]` looks up genome names before display names, so a display name can no longer hide another genome.
- `analysis.json` `warnings` now include the GenBank parse warnings of the reference and the target (prefixed `reference: ` and `target: `), before the warnings of the run. `Genome.parse_warnings` and `ReferenceRegion.parse_warnings` keep them in the API.
- Printing a genome, contig, reference or results shows a one-line summary instead of every sequence and the raw DIAMOND output, so they stay readable in notebooks. Genes are printed without their sequences and qualifiers. Equality is unchanged.
- An interrupted `homolocus find` (Ctrl-C) prints `homolocus: interrupted` and exits with 130 instead of a traceback.

## [0.1.0] - Unreleased

The first release. Thresholds and score weights are unvalidated starting values; see section 18 of the specification.

### Added

- **Inputs:**
  - reference regions from GenBank, as a whole genome with a contig and 1-based coordinates, or as a region-only file;
  - target assemblies from GenBank or FASTA, with the format detected from content;
  - Pyrodigal gene prediction for FASTA and for GenBank without CDS (`auto` or `pyrodigal` annotation);
  - CDS translation rules, with skipped pseudogenes, internal stops and compound CDS recorded in the run metadata.
- **Search:**
  - DIAMOND protein search (2.0 or later) that keeps every hit and HSP, aggregated per gene pair with union-based coverage;
  - paralogous reference genes grouped into families from a DIAMOND self-search.
- **Detection:**
  - candidate regions found per contig from chains and windows of matching events, tolerant of missing and inserted genes, duplicated homologues, rearranged order and reversed regions;
  - a weighted score of content, similarity, order, orientation and spacing, with every component exposed;
  - flags for split proteins, fused proteins, partial genes and candidates at a contig edge.
- **Results:** `RegionResults`, with ranked candidates, `to_dataframe()` and `genes_dataframe()`.
- **Output:** a deterministic output directory with `summary.tsv`, per-candidate TSVs, raw DIAMOND tables, `analysis.json` and PNG plots.
- **Plots:** Matplotlib overview and per-candidate plots, with consistent family colours and labels.
- **Command line:** `homolocus find`, with every setting as a flag, exit codes 0, 1 and 2, and a stdout summary.
- **Requirements:** Python 3.10–3.13. Tested minimum versions are Biopython 1.80, pandas 2.0, Matplotlib 3.8 and Pyrodigal 3.0.

[0.1.0]: https://github.com/emilhaegglund/homolocus/releases/tag/v0.1.0
