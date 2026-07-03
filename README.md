# rblat

**In-process BLAT sequence alignment for R.**

`rblat` runs the [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html) aligner (Jim
Kent; the [djhshih/blat](https://github.com/djhshih/blat) fork of v35) entirely
inside R through a compiled C backend driven by `.Call`. There is **no
`system()` call to an external `blat` binary**, and — if you pass sequences in
memory — no temporary files at all. Results come back as the same 21-column PSL
`data.frame` (`blat_colnames`) that the command-line `blat` + `readBlat()`
produce.

It is intended as a drop-in, in-process replacement for the command-line blat
step used by the [`trellis`](https://github.com/cancer-genomics/trellis)
structural-variant pipeline.

## Why

The command-line workflow — write a query FASTA, spawn `blat`, read the `.psl`
back — requires `blat` on the `PATH` and pays process-spawn and disk-I/O
overhead on every call. `rblat` links the BLAT engine directly into the package
shared library, so the same alignment happens in-process and returns a
`data.frame` directly.

## Installation

Requires a C toolchain and `zlib` (both standard on macOS/Linux).

```r
# from a local checkout
R CMD INSTALL rblat
# or, from R:
# install.packages("rblat", repos = NULL, type = "source")
```

## Quick start

```r
library(rblat)

# --- FASTA files (drop-in for the command-line call) ---
aln <- rblat("reference.fa", "reads.fa", minScore = 0, stepSize = 1)
head(aln)
names(aln)           # == blat_colnames (21-column PSL layout)

# --- in-memory sequences (no files touched) ---
ref   <- c(chr1 = "ACGT...")                 # named character vector
reads <- c(read1 = "ACG...", read2 = "TTT...")
rblat(ref, reads, minScore = 30)

# --- Biostrings input ---
rblat(Biostrings::DNAStringSet(ref),
      Biostrings::DNAStringSet(reads))
```

`database` and `query` may each be **FASTA file path(s)**, a **named character
vector** (names = sequence identifiers), or a **`Biostrings::DNAStringSet`** —
both arguments must use the same style. A query with no acceptable alignment
returns a well-typed *empty* `data.frame` (not an error).

## Exported functions

| Function | Purpose |
|----------|---------|
| `rblat(database, query, ...)` | Align in-process; returns a PSL `data.frame`. |
| `readPsl(filename)` | Read a `.psl` file into the same `data.frame` (handles the empty-result case correctly). |
| `blat_colnames` | The 21 PSL column names, matching `trellis::blat_colnames`. |
| `rblat_version()` | Confirm the vendored BLAT engine is linked. |

### Options

`rblat()` exposes BLAT's common tuning options as named arguments
(`minScore`, `minMatch`, `tileSize`, `stepSize`, `minIdentity`, `maxGap`,
`maxIntron`, `repMatch`, `t`, `q`, `mask`, `qMask`, `ooc`, `oneOff`, `fastMap`,
plus raw pass-through `args`). Each numeric option defaults to `NULL`, meaning
"use BLAT's own *type-aware* default" (e.g. `tileSize` is 11 for DNA but 5 for
protein). See `?rblat`.

Current limits: in-memory input supports simple **dna/dna** and **prot/prot**;
translated modes (`dnax`/`rnax`) and file-based target masking require FASTA
file input. See `vignettes/rblat.Rmd` and `DEVELOPMENT.md`.

## The `trellis` use-case

`trellis` re-aligns candidate reads with blat to check whether they span a novel
sequence junction. That step currently does:

```r
writeToFasta(tags, "reads.fa")
system("blat -minScore=0 -stepSize=1 reference.fa reads.fa out.psl")
blat <- readBlat("out.psl")
blat <- annotateBlatRecords(blat, tags)   # downstream QC unchanged
```

With `rblat` the middle three steps collapse to one call returning the identical
`data.frame`:

```r
blat <- rblat("reference.fa", "reads.fa", minScore = 0, stepSize = 1)
```

The package bundles `trellis`/`svbams`'s own blat fixtures together with the
real command-line-blat goldens, and its test suite verifies that `rblat`
reproduces them **byte-for-byte** — with no dependency on `trellis`. The
vignette (`vignettes/rblat.Rmd`) walks through this and sketches a future
file-free integration.

## License

`rblat` links the BLAT aligner core directly and therefore inherits the
**original BLAT license**: free for personal, academic, and non-profit use;
**commercial use requires a separate license** from
[Kent Informatics](https://kentinformatics.com). This is the same licensing that
governs the command-line `blat` binary. See [`LICENSE`](LICENSE) for the
verbatim terms.

## Status

Engine and R interface are complete (in-memory I/O, options/modes, robustness,
and trellis-equivalence regression tests; 112-check `testthat` suite). See
[`DEVELOPMENT.md`](DEVELOPMENT.md) for phase-by-phase build notes.
