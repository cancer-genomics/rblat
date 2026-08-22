# trellis blat regression fixtures

These files reproduce the exact command-line `blat` step that the `trellis`
package (and its `svbams` test-data package) relies on. They let rblat's test
suite prove, without any dependency on or edits to `trellis`, that `rblat()` is
a drop-in replacement for that step.

## Provenance

Copied verbatim from `svbams` (Cancer Genomics lab, `cancer_genomics/svbams`):

| rblat file    | source in svbams                        |
|---------------|-----------------------------------------|
| `ref1.fa`     | `inst/extdata/fake_reference.fa`        |
| `query1.fa`   | `inst/extdata/fake_query.fa`            |
| `ref2.fa`     | `inst/extdata/fake_reference2.fa`       |
| `query2.fa`   | `inst/extdata/fake_query2.fa`           |
| `golden1.psl` | `data-raw/test_alignment.txt`           |
| `golden2.psl` | `data-raw/test_alignment2.txt`          |

The two goldens were produced by real, external command-line `blat` (v35,
sourceforge build) via `svbams/data-raw/blat_test.R`, invoked as:

```
blat -minScore=0 -stepSize=1 <ref> <query> <out>
```

read back into trellis with `readBlat()` (`read.delim(skip = 5)` + `blat_colnames`).

## What the regression test asserts

`rblat(ref, query, minScore = 0, stepSize = 1)` returns a data.frame **equal**
to `readPsl(golden)` — i.e. byte-for-byte identical to what command-line blat
produced and what trellis consumes today. It also checks the two invariants the
svbams script documents: query1 yields two single-block alignments; query2
yields `blockSizes == "50,50,"`.
