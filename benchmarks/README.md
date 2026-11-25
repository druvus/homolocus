# Benchmarks

Checks against real, public genomes. They live outside the test suite: they
need network access and DIAMOND, they take minutes rather than seconds, and
they report what they found instead of asserting it (spec 17). `pytest`
never runs them.

What they show informs section 18 of `SPECIFICATION.md`, where the
thresholds and score weights are still open. No default is changed on the
strength of a benchmark without discussing it first.

## lee.py

The locus of enterocyte effacement (LEE) is a pathogenicity island present
in enterohaemorrhagic and enteropathogenic *Escherichia coli* and absent
from commensal K-12 and uropathogenic CFT073, so the expected answers are
known.

```bash
python benchmarks/lee.py
```

It downloads about 45 MB from NCBI into `benchmarks/data/` (git-ignored,
reused on later runs) and writes `benchmarks/results/lee.md`.

The reference region is derived from the Sakai annotation, not hard-coded:
the densest cluster of CDS carrying LEE gene names defines it, and the
report says which genes those were. Targets are read as FASTA, so Pyrodigal
predicts their genes and no result depends on how each genome was
annotated. The genome the reference came from is searched too, as a
control: it should be found with a near-perfect score.

The committed report records one run, with the homolocus and DIAMOND
versions it used.
