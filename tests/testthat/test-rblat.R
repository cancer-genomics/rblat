ed <- function(...) system.file("extdata", ..., package = "rblat")

test_that("engine is linked and callable", {
  expect_match(rblat_version(), "engine linked")
})

test_that("small case matches golden (in-process, no system call)", {
  df <- rblat(ed("small/target.fa"), ed("small/query.fa"),
              minScore = 20, minMatch = 1)
  golden <- readPsl(ed("small/golden.psl"))
  expect_equal(df, golden)
  expect_identical(names(df), blat_colnames)
})

test_that("v29skips ex2 matches golden with default options", {
  df <- rblat(ed("v29skips/ex2_database.fa"), ed("v29skips/ex2_query.fa"))
  golden <- readPsl(ed("v29skips/ex2_golden.psl"))
  expect_equal(nrow(df), 1L)
  expect_equal(df, golden)
})

test_that("results are repeatable across calls (option-state reset)", {
  a <- rblat(ed("small/target.fa"), ed("small/query.fa"), minScore = 20, minMatch = 1)
  b <- rblat(ed("v29skips/ex2_database.fa"), ed("v29skips/ex2_query.fa"))
  a2 <- rblat(ed("small/target.fa"), ed("small/query.fa"), minScore = 20, minMatch = 1)
  expect_equal(a, a2)
})

test_that("no-hit query yields an empty but well-typed data.frame", {
  df <- rblat(ed("small/target.fa"), ed("small/query.fa"))  # 25bp < default minScore 30
  expect_equal(nrow(df), 0L)
  expect_identical(names(df), blat_colnames)
  expect_type(df[["match"]], "integer")
})

test_that("malformed input is a clean R error, not a crashed session", {
  missing_query <- tempfile(fileext = ".fa")  # does not exist
  expect_error(rblat(ed("small/target.fa"), missing_query))
  # session survives:
  expect_equal(1 + 1, 2)
})

# ---- Phase 4: in-memory sequence input (no FASTA files) --------------------

# Minimal FASTA -> named character vector (uses first whitespace-delimited
# token of each header, matching how blat names sequences).
.read_fasta_vec <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(ln)]
  h <- grepl("^>", ln)
  nm <- sub("[[:space:]].*$", "", sub("^>", "", ln[h]))
  seqs <- tapply(ln[!h], cumsum(h)[!h], paste, collapse = "")
  stats::setNames(as.character(seqs), nm)
}

test_that("in-memory (named vector) input matches file input and golden", {
  db <- .read_fasta_vec(ed("small/target.fa"))
  q  <- .read_fasta_vec(ed("small/query.fa"))
  mem  <- rblat(db, q, minScore = 20, minMatch = 1)
  file <- rblat(ed("small/target.fa"), ed("small/query.fa"), minScore = 20, minMatch = 1)
  expect_equal(mem, file)
  expect_equal(mem, readPsl(ed("small/golden.psl")))
})

test_that("in-memory ex2 matches golden and is repeatable across calls", {
  db <- .read_fasta_vec(ed("v29skips/ex2_database.fa"))
  q  <- .read_fasta_vec(ed("v29skips/ex2_query.fa"))
  a <- rblat(db, q)
  # interleave a different in-memory call, then repeat
  smalldb <- .read_fasta_vec(ed("small/target.fa"))
  smallq  <- .read_fasta_vec(ed("small/query.fa"))
  invisible(rblat(smalldb, smallq, minScore = 20, minMatch = 1))
  a2 <- rblat(db, q)
  expect_equal(a, readPsl(ed("v29skips/ex2_golden.psl")))
  expect_equal(a, a2)
})

test_that("in-memory no-hit query yields an empty but well-typed data.frame", {
  db <- .read_fasta_vec(ed("small/target.fa"))
  q  <- .read_fasta_vec(ed("small/query.fa"))
  df <- rblat(db, q)  # 25bp < default minScore 30
  expect_equal(nrow(df), 0L)
  expect_identical(names(df), blat_colnames)
  expect_type(df[["match"]], "integer")
})

test_that("in-memory input works via a Biostrings DNAStringSet", {
  skip_if_not_installed("Biostrings")
  db <- Biostrings::DNAStringSet(.read_fasta_vec(ed("small/target.fa")))
  q  <- Biostrings::DNAStringSet(.read_fasta_vec(ed("small/query.fa")))
  df <- rblat(db, q, minScore = 20, minMatch = 1)
  expect_equal(df, readPsl(ed("small/golden.psl")))
})

test_that("input dispatch rejects mixed / unnamed / unsupported inputs cleanly", {
  db <- .read_fasta_vec(ed("small/target.fa"))
  q  <- .read_fasta_vec(ed("small/query.fa"))
  # one file path, one in-memory
  expect_error(rblat(db, ed("small/query.fa")), "both")
  # unnamed character that isn't a file path
  expect_error(rblat(c("ACGT", "TTTT"), q), "named")
  # translated modes not supported in memory (use valid protein tileSize so we
  # reach the guard rather than blat's tileSize check)
  expect_error(rblat(db, q, t = "dnax", q = "dnax", tileSize = 5), "translated")
  expect_equal(1 + 1, 2)  # session survives all of the above
})

# ---- Phase 5: options & modes + edge cases --------------------------------

.rnd_dna <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")

test_that("default options reproduce ex2 golden without any numeric flags", {
  # rblat's numeric options now default to NULL (BLAT's own type-aware defaults);
  # the DNA default run must still match the golden generated with blat defaults.
  df <- rblat(ed("v29skips/ex2_database.fa"), ed("v29skips/ex2_query.fa"))
  expect_equal(df, readPsl(ed("v29skips/ex2_golden.psl")))
})

test_that("a spliced query yields a correct multi-block (intron) alignment", {
  set.seed(1)
  exon1 <- .rnd_dna(80); intron <- .rnd_dna(2000); exon2 <- .rnd_dna(80)
  target <- c(chr = paste0(.rnd_dna(50), exon1, intron, exon2, .rnd_dna(50)))
  mrna   <- c(tx  = paste0(exon1, exon2))  # spliced transcript (no intron)
  df <- rblat(target, mrna, minScore = 30, minMatch = 1)
  expect_equal(nrow(df), 1L)
  expect_equal(df$blockcount, 2L)
  expect_equal(df$match, 160L)
  # each block field is comma-terminated with exactly blockcount entries
  split_ints <- function(s) as.integer(strsplit(s, ",", fixed = TRUE)[[1]])
  bs <- split_ints(df$blockSizes); qs <- split_ints(df$qStarts); ts <- split_ints(df$tStarts)
  expect_length(bs, 2); expect_length(qs, 2); expect_length(ts, 2)
  expect_equal(sum(bs), df$match)                 # blocks account for all matches
  expect_equal(qs, c(0L, 80L))                    # query is contiguous across the join
  expect_gt(ts[2] - (ts[1] + bs[1]), 1000L)       # target gap == the intron
})

test_that("minIdentity filters low-identity alignments", {
  set.seed(2)
  base <- .rnd_dna(100)
  target <- c(chr = paste0(.rnd_dna(30), base, .rnd_dna(30)))
  qs <- strsplit(base, "")[[1]]
  pos <- c(10, 30, 50, 70, 90)                    # 5 mismatches -> ~95% identity
  qs[pos] <- ifelse(qs[pos] == "A", "C", "A")
  query <- c(q = paste(qs, collapse = ""))
  expect_equal(nrow(rblat(target, query, minIdentity = 90, minScore = 30)), 1L)
  expect_equal(nrow(rblat(target, query, minIdentity = 98, minScore = 30)), 0L)
})

test_that("protein/protein search works in memory using type-aware defaults", {
  # No tileSize given: BLAT must pick its protein default (5), not DNA's 11.
  pep <- "MKVLADESGHIKWYACDEFGHIKLMNPQRSTVWYMKVLADESG"
  target <- c(prot1 = paste0("AAA", pep, "GGG"))
  query  <- c(pquery = pep)
  df <- rblat(target, query, t = "prot", q = "prot")
  expect_true(nrow(df) >= 1L)
  top <- df[which.max(df$match), ]
  expect_equal(top$match, nchar(pep))
  expect_equal(top$Qname, "pquery")
  expect_equal(top$strand, "+")
})

test_that("file-based / index options are wired through and validated", {
  set.seed(3)
  target <- c(chr = .rnd_dna(200)); query <- c(q = substring(target, 40, 120))
  # repMatch and tileSize are accepted and (here) don't change a clean hit
  base <- rblat(target, query, minScore = 30, minMatch = 1)
  expect_equal(rblat(target, query, minScore = 30, minMatch = 1, repMatch = 256), base)
  expect_equal(rblat(target, query, minScore = 30, minMatch = 1, tileSize = 8), base)
  # qMask (case-based) is allowed with in-memory input
  expect_s3_class(rblat(target, query, qMask = "lower"), "data.frame")
  # target mask needs files; missing ooc file -> clean R errors, session survives
  expect_error(rblat(target, query, mask = "lower"), "mask")
  expect_error(rblat(target, query, ooc = tempfile(fileext = ".ooc")), "ooc")
  expect_equal(1 + 1, 2)
})

# ---- Phase 6: robustness (exit/longjmp safety, warnings, reentrancy) -------

test_that("results stay correct after repeated caught errors (longjmp cleanup)", {
  set.seed(4)
  target <- c(chr = .rnd_dna(300)); query <- c(q = substring(target, 50, 180))
  good <- rblat(target, query, minScore = 30, minMatch = 1)
  expect_gt(nrow(good), 0L)
  # Hammer the engine's errAbort/longjmp path, which must disarm the collector
  # and in-memory inputs each time, then confirm a normal run is still correct.
  for (i in 1:25) {
    expect_error(rblat(target, query, t = "dnax", q = "dnax", tileSize = 5))  # translated -> errAbort
    expect_error(rblat(target, query, mask = "out"))                          # rejected in R
  }
  expect_equal(rblat(target, query, minScore = 30, minMatch = 1), good)
})

test_that("BLAT warn() is surfaced as an R warning, run still succeeds", {
  # A query longer than BLAT's qWarnSize (5,000,000) triggers a warn(); all-'N'
  # keeps the search itself trivial.
  db <- c(chr = paste(rep("N", 100L), collapse = ""))
  q  <- c(bigq = paste(rep("N", 5000001L), collapse = ""))
  expect_warning(df <- rblat(db, q, minScore = 30), "take a while")
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0L)
})

test_that("a normal run emits no spurious warnings", {
  set.seed(5)
  target <- c(chr = .rnd_dna(300)); query <- c(q = substring(target, 50, 180))
  expect_no_warning(rblat(target, query, minScore = 30, minMatch = 1))
})

# ---- Phase 7: trellis integration (drop-in for the command-line blat step) --
#
# The trellis package aligns unmapped reads with an EXTERNAL command-line blat
# and reads the .psl with readBlat() (== read.delim(skip = 5) + blat_colnames).
# svbams/data-raw/blat_test.R runs `blat -minScore=0 -stepSize=1 <ref> <query>`
# on two small fixtures and checks its output. We vendored those fixtures and
# their real command-line-blat goldens into inst/extdata/trellis (see the
# README there) so we can prove rblat() reproduces that step byte-for-byte --
# with NO dependency on, or edits to, trellis.

test_that("rblat reproduces command-line blat on trellis fixture 1 (two hits)", {
  # trellis invocation: blat -minScore=0 -stepSize=1 ref query out
  df  <- rblat(ed("trellis/ref1.fa"), ed("trellis/query1.fa"),
               minScore = 0, stepSize = 1)
  gld <- readPsl(ed("trellis/golden1.psl"))   # real command-line blat output
  expect_equal(df, gld)                        # byte-for-byte == what trellis reads
  expect_identical(names(df), blat_colnames)
  # invariant documented in svbams/data-raw/blat_test.R: two single-block hits
  expect_equal(nrow(df), 2L)
  expect_equal(df$blockcount, c(1L, 1L))
})

test_that("rblat reproduces command-line blat on trellis fixture 2 (split hit)", {
  df  <- rblat(ed("trellis/ref2.fa"), ed("trellis/query2.fa"),
               minScore = 0, stepSize = 1)
  gld <- readPsl(ed("trellis/golden2.psl"))
  expect_equal(df, gld)
  # invariant documented in svbams: a single two-block alignment, "50,50,"
  expect_equal(nrow(df), 1L)
  expect_equal(df$blockSizes, "50,50,")
})

test_that("in-memory input matches the trellis file-based result and golden", {
  # Same fixtures via named character vectors (whitespace stripped, as the
  # svbams fixtures wrap sequences with a leading space). Proves the in-memory
  # path is equally a drop-in for trellis's blat step.
  clean <- function(v) stats::setNames(gsub("[[:space:]]", "", v), names(v))
  ref1 <- clean(.read_fasta_vec(ed("trellis/ref1.fa")))
  qry1 <- clean(.read_fasta_vec(ed("trellis/query1.fa")))
  mem  <- rblat(ref1, qry1, minScore = 0, stepSize = 1)
  file <- rblat(ed("trellis/ref1.fa"), ed("trellis/query1.fa"),
                minScore = 0, stepSize = 1)
  expect_equal(mem, file)
  expect_equal(mem, readPsl(ed("trellis/golden1.psl")))
})

test_that("rblat output is column-compatible with trellis::readBlat", {
  # trellis::readBlat is read.delim(skip = 5, header = FALSE) + blat_colnames.
  # Emulate it locally (no trellis dependency) on the vendored golden and
  # confirm rblat() yields an identical data.frame -- so downstream trellis
  # code (annotateBlatRecords, blatScores, blatGRanges, ...) sees the same
  # object whether fed by command-line blat or rblat.
  readBlat_like <- function(path) {
    b <- utils::read.delim(path, skip = 5, header = FALSE,
                           sep = "\t", stringsAsFactors = FALSE)
    colnames(b) <- blat_colnames
    b
  }
  df  <- rblat(ed("trellis/ref2.fa"), ed("trellis/query2.fa"),
               minScore = 0, stepSize = 1)
  ref <- readBlat_like(ed("trellis/golden2.psl"))
  expect_identical(names(df), names(ref))
  expect_equal(df, ref)
})
