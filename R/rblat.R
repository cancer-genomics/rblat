#' Report that the vendored BLAT engine is linked and callable
#' @return A character string.
#' @export
rblat_version <- function() {
  .Call(C_rblat_version)
}

#' Standard 21 PSL column names (matches trellis::blat_colnames)
#' @export
blat_colnames <- c(
  "match", "mis-match", "rep.match", "N's",
  "Qgapcount", "Qgapbases", "Tgapcount", "Tgapbases",
  "strand", "Qname", "Qsize", "Qstart", "Qend",
  "Tname", "Tsize", "Tstart", "Tend",
  "blockcount", "blockSizes", "qStarts", "tStarts")

# Columns that are integer counts/coordinates (everything except these is character)
.psl_int_cols <- c(
  "match", "mis-match", "rep.match", "N's",
  "Qgapcount", "Qgapbases", "Tgapcount", "Tgapbases",
  "Qsize", "Qstart", "Qend", "Tsize", "Tstart", "Tend", "blockcount")

#' Read a PSL file into a data.frame with BLAT column names
#'
#' Unlike a naive `read.delim(skip = 5)`, this correctly returns an empty (but
#' well-typed) data.frame when BLAT found no alignments (the file then contains
#' only the 5-line psLayout header).
#'
#' @param filename Path to a `.psl` file.
#' @return A data.frame with columns [blat_colnames].
#' @export
readPsl <- function(filename) {
  lines <- readLines(filename, warn = FALSE)
  data_lines <- if (length(lines) > 5L) lines[-seq_len(5L)] else character(0)
  data_lines <- data_lines[nzchar(data_lines)]
  if (length(data_lines) == 0L) {
    df <- as.data.frame(
      lapply(blat_colnames, function(nm)
        if (nm %in% .psl_int_cols) integer(0) else character(0)),
      stringsAsFactors = FALSE)
    names(df) <- blat_colnames
    return(df)
  }
  df <- utils::read.delim(text = data_lines, header = FALSE, sep = "\t",
                          stringsAsFactors = FALSE, colClasses = "character")
  names(df) <- blat_colnames
  for (nm in .psl_int_cols) df[[nm]] <- as.integer(df[[nm]])
  df
}

# Interpret an rblat sequence argument. Returns a list(names, seqs) for
# in-memory input, or NULL when `x` is one or more existing FASTA file paths
# (file mode). Accepts a named character vector (names = identifiers, values =
# sequences) or a Biostrings DNAStringSet / XStringSet.
.rblat_seq_input <- function(x, what) {
  if (methods::is(x, "XStringSet")) {
    seqs <- as.character(x)
    nms <- names(seqs)
    if (is.null(nms)) nms <- names(x)
    if (is.null(nms)) nms <- paste0("seq", seq_along(seqs))
    return(list(names = as.character(nms), seqs = unname(as.character(seqs))))
  }
  if (is.character(x)) {
    nms <- names(x)
    if (!is.null(nms)) {
      if (anyNA(nms) || any(!nzchar(nms)))
        stop(sprintf("`%s`: every sequence must have a non-empty name", what))
      return(list(names = nms, seqs = unname(x)))
    }
    if (length(x) >= 1L && all(file.exists(x)))
      return(NULL)  # file mode
    stop(sprintf(paste0("`%s` must be either existing FASTA file path(s) or a ",
                        "*named* character vector of sequences (names = ",
                        "sequence identifiers)"), what))
  }
  stop(sprintf("`%s`: unsupported input type '%s'", what, class(x)[1]))
}

#' Run BLAT in-process and return a PSL data.frame
#'
#' Calls the BLAT engine compiled into this package via `.Call` -- no `system()`
#' call to an external `blat` binary, and no PSL file: alignments are collected
#' into memory by the C backend and returned as a data.frame directly.
#'
#' `database` and `query` may be given either as **in-memory sequences** (a
#' named character vector, names = sequence identifiers, or a Biostrings
#' `DNAStringSet`) or as **FASTA file path(s)**. Both must use the same style.
#' With in-memory sequences the entire run touches no files at all.
#'
#' @param database,query Target/database and query sequences: a named character
#'   vector, a `DNAStringSet`, or FASTA file path(s).
#' @param minScore,minMatch,tileSize,stepSize,minIdentity,maxGap,maxIntron,repMatch
#'   Numeric BLAT tuning options. Each defaults to `NULL`, meaning "use BLAT's
#'   own default" -- which is *type-aware*: e.g. `tileSize` is 11 for DNA but 5
#'   for protein, `minMatch` 2 vs 1, `minIdentity` 90 vs 25. Set a value to
#'   override. (Forcing DNA numbers on a protein/translated search is why these
#'   are not hard-coded here.)
#' @param t,q Target/query sequence type: one of "dna", "prot", "dnax", "rnax",
#'   "rna". In-memory input currently supports only the simple dna/dna and
#'   prot/prot combinations; translated modes (dnax/rnax, protein-vs-dna) require
#'   FASTA file input for now.
#' @param mask,qMask Repeat masking for target / query. Case-based values
#'   ("lower", "upper") work everywhere; file-based target masking
#'   ("out", "file.out") needs FASTA file input (`mask` is rejected with
#'   in-memory sequences).
#' @param ooc Path to an over-occurring N-mer (`.ooc`) file to speed up indexing.
#' @param oneOff,fastMap Logical BLAT flags.
#' @param args Character vector of extra raw BLAT flags, e.g. `c("-fine")`.
#' @return A data.frame of PSL records with columns [blat_colnames].
#' @export
rblat <- function(database, query,
                  minScore = NULL, minMatch = NULL, tileSize = NULL,
                  stepSize = NULL, minIdentity = NULL, maxGap = NULL,
                  maxIntron = NULL, repMatch = NULL,
                  t = NULL, q = NULL,
                  mask = NULL, qMask = NULL, ooc = NULL,
                  oneOff = FALSE, fastMap = FALSE,
                  args = character(0)) {
  dbIn <- .rblat_seq_input(database, "database")
  qIn  <- .rblat_seq_input(query, "query")
  if (is.null(dbIn) != is.null(qIn))
    stop("`database` and `query` must both be file paths or both be in-memory sequences")
  in_memory <- !is.null(dbIn)

  # Options that need real files can't run on in-memory sequences; fail early
  # with a clear message (the C engine would otherwise errAbort).
  if (in_memory && !is.null(mask))
    stop("`mask` (target masking) requires FASTA file input; not supported with ",
         "in-memory sequences")
  if (in_memory && !is.null(ooc) && !file.exists(ooc))
    stop("`ooc` file does not exist: ", ooc)

  # Only emit a flag when the user set the option, so BLAT applies its own
  # type-aware default otherwise.
  int_opt <- function(name, val)
    if (is.null(val)) character(0) else sprintf("-%s=%d", name, as.integer(val))
  num_opt <- function(name, val)
    if (is.null(val)) character(0) else sprintf("-%s=%g", name, val)
  str_opt <- function(name, val)
    if (is.null(val)) character(0) else sprintf("-%s=%s", name, val)

  flags <- c(
    "-verbose=0",
    int_opt("minScore", minScore),
    int_opt("minMatch", minMatch),
    int_opt("tileSize", tileSize),
    int_opt("stepSize", stepSize),
    num_opt("minIdentity", minIdentity),
    int_opt("maxGap", maxGap),
    int_opt("maxIntron", maxIntron),
    int_opt("repMatch", repMatch),
    str_opt("t", t),
    str_opt("q", q),
    str_opt("mask", mask),
    str_opt("qMask", qMask),
    str_opt("ooc", ooc),
    if (isTRUE(oneOff))  "-oneOff"  else character(0),
    if (isTRUE(fastMap)) "-fastMap" else character(0),
    args)

  if (!in_memory) {
    # File mode: database/query are FASTA paths.
    argv <- c("blat", flags, database, query)
    return(.Call(C_rblat_run_collect, argv))
  }
  # In-memory mode: placeholder positional args satisfy blat's argc check.
  argv <- c("blat", flags, "rblat_db", "rblat_query")
  .Call(C_rblat_run_collect_mem, argv,
        dbIn$names, dbIn$seqs, qIn$names, qIn$seqs)
}
