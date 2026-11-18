# Changelog

All notable changes to homolocus are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Families are named from the reference GenBank annotation: their genes' `/gene` qualifiers, else `/product`. Names appear in plot labels and legends (`F001 sitA`), in the missing-families note, in a new `family_name` column of the candidate TSVs, and on `GeneFamily.name` and `FamilyEvidence.family_name`.
- Several target genomes in one run:
  - `--genome` takes one or more files, and `RegionFinder.find_all()` returns `GenomeSetResults`;
  - families and the self-search are computed once, and each genome's results equal a run on it alone;
  - output has a folder per genome, a combined `summary.tsv` with a `genome` column, and a combined overview plot.
- Candidates are called `forward` or `reversed` relative to the reference (gene order, then strand). The call is in `Cluster.reversed` and a new `direction` column of `summary.tsv`. Plots draw every candidate in the reference direction: a reversed one sits on a descending axis, marked `reversed`.

### Changed

- Genomes read from files are named after the file (`Genome.name`), so `analysis.json` now records `target.name` instead of `null`.

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
