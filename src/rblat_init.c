/* rblat C entry points and .Call registration.
 *
 * This file includes ONLY R headers. The BLAT/Kent headers redefine
 * TRUE/FALSE/boolean and clash with Rinternals.h, so we reach the engine through
 * the plain-C prototypes in rblat.h (no Kent types cross this boundary).
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <string.h>
#include "rblat.h"

/* From rblat_blat.c: initializes nucleotide lookup tables. Declared here (rather
 * than via dnautil.h) to avoid pulling in Kent headers. */
extern void dnaUtilOpen(void);

/* Re-emit any warn() text the engine captured during a successful run as an R
 * warning, then free the malloc'd string. No-op on NULL. Trailing newline (the
 * engine appends one per warn()) is trimmed so the R warning reads cleanly. */
static void rblat_emit_warning(char *warnMsg)
{
    if (warnMsg == NULL)
        return;
    size_t n = strlen(warnMsg);
    while (n > 0 && warnMsg[n - 1] == '\n')
        warnMsg[--n] = '\0';
    if (n > 0)
        Rf_warning("blat: %s", warnMsg);
    free(warnMsg);
}

SEXP rblat_version(void)
{
    dnaUtilOpen();
    return Rf_mkString("rblat: BLAT v35 engine linked (Phase 1 stub)");
}

/* Run blat given a character vector argv (including argv[0] = "blat").
 * optionInit() mutates the argv strings and array in place, so we must pass
 * WRITABLE copies -- never the CHAR() pointers into R's string pool. */
SEXP rblat_run(SEXP r_argv)
{
    if (TYPEOF(r_argv) != STRSXP)
        Rf_error("rblat_run: argv must be a character vector");
    int argc = LENGTH(r_argv);
    if (argc < 4)
        Rf_error("rblat_run: need at least argv0 + database + query + output");

    /* Writable copies in R_alloc memory (freed automatically at .Call return). */
    char **argv = (char **) R_alloc((size_t) argc, sizeof(char *));
    for (int i = 0; i < argc; i++) {
        const char *s = CHAR(STRING_ELT(r_argv, i));
        size_t n = strlen(s) + 1;
        char *buf = R_alloc(n, 1);
        memcpy(buf, s, n);
        argv[i] = buf;
    }

    char *errMsg = NULL;
    int ok = rblat_run_argv(argc, argv, &errMsg);
    if (!ok) {
        char msg[1024];
        snprintf(msg, sizeof(msg), "blat error: %s",
                 errMsg ? errMsg : "(no message)");
        if (errMsg) free(errMsg);
        Rf_error("%s", msg);          /* does not return */
    }
    return R_NilValue;
}

/* The 21 PSL columns, in blat_colnames order. Kept in sync with R/rblat.R. */
static const char *rblat_colnames[] = {
    "match", "mis-match", "rep.match", "N's",
    "Qgapcount", "Qgapbases", "Tgapcount", "Tgapbases",
    "strand", "Qname", "Qsize", "Qstart", "Qend",
    "Tname", "Tsize", "Tstart", "Tend",
    "blockcount", "blockSizes", "qStarts", "tStarts"};
#define RBLAT_NCOL 21

static SEXP rblat_int_col(const int *src, int n)
/* Build an integer vector column from a C int array. */
{
    SEXP v = PROTECT(Rf_allocVector(INTSXP, n));
    if (n > 0)
        memcpy(INTEGER(v), src, (size_t) n * sizeof(int));
    UNPROTECT(1);
    return v;
}

static SEXP rblat_chr_col(char **src, int n)
/* Build a character vector column from a C char* array. */
{
    SEXP v = PROTECT(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; i++)
        SET_STRING_ELT(v, i, Rf_mkChar(src[i]));
    UNPROTECT(1);
    return v;
}

/* Build the 21-column PSL data.frame from a collected result. The column set,
 * order, and types match readPsl()'s output so goldens compare equal. */
static SEXP rblat_result_to_df(struct rblatResult *r)
{
    int n = r->nrow;
    SEXP df    = PROTECT(Rf_allocVector(VECSXP, RBLAT_NCOL));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, RBLAT_NCOL));
    for (int j = 0; j < RBLAT_NCOL; j++)
        SET_STRING_ELT(names, j, Rf_mkChar(rblat_colnames[j]));

    SET_VECTOR_ELT(df,  0, rblat_int_col(r->match,       n));
    SET_VECTOR_ELT(df,  1, rblat_int_col(r->misMatch,    n));
    SET_VECTOR_ELT(df,  2, rblat_int_col(r->repMatch,    n));
    SET_VECTOR_ELT(df,  3, rblat_int_col(r->nCount,      n));
    SET_VECTOR_ELT(df,  4, rblat_int_col(r->qNumInsert,  n));
    SET_VECTOR_ELT(df,  5, rblat_int_col(r->qBaseInsert, n));
    SET_VECTOR_ELT(df,  6, rblat_int_col(r->tNumInsert,  n));
    SET_VECTOR_ELT(df,  7, rblat_int_col(r->tBaseInsert, n));
    SET_VECTOR_ELT(df,  8, rblat_chr_col(r->strand,      n));
    SET_VECTOR_ELT(df,  9, rblat_chr_col(r->qName,       n));
    SET_VECTOR_ELT(df, 10, rblat_int_col(r->qSize,       n));
    SET_VECTOR_ELT(df, 11, rblat_int_col(r->qStart,      n));
    SET_VECTOR_ELT(df, 12, rblat_int_col(r->qEnd,        n));
    SET_VECTOR_ELT(df, 13, rblat_chr_col(r->tName,       n));
    SET_VECTOR_ELT(df, 14, rblat_int_col(r->tSize,       n));
    SET_VECTOR_ELT(df, 15, rblat_int_col(r->tStart,      n));
    SET_VECTOR_ELT(df, 16, rblat_int_col(r->tEnd,        n));
    SET_VECTOR_ELT(df, 17, rblat_int_col(r->blockCount,  n));
    SET_VECTOR_ELT(df, 18, rblat_chr_col(r->blockSizes,  n));
    SET_VECTOR_ELT(df, 19, rblat_chr_col(r->qStarts,     n));
    SET_VECTOR_ELT(df, 20, rblat_chr_col(r->tStarts,     n));

    Rf_setAttrib(df, R_NamesSymbol, names);

    /* Compact automatic row.names c(NA, -n), as produced by read.table/
     * as.data.frame -- keeps equality with readPsl()'s data.frame. */
    SEXP rn = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -n;
    Rf_setAttrib(df, R_RowNamesSymbol, rn);

    Rf_setAttrib(df, R_ClassSymbol, Rf_mkString("data.frame"));
    UNPROTECT(3);
    return df;
}

/* Run blat in-process and return the PSL records as a data.frame directly --
 * no PSL text and no output file (Phase 3). argv is argv0 + database + query
 * (no trailing output path). */
SEXP rblat_run_collect(SEXP r_argv)
{
    if (TYPEOF(r_argv) != STRSXP)
        Rf_error("rblat_run_collect: argv must be a character vector");
    int argc = LENGTH(r_argv);
    if (argc < 3)
        Rf_error("rblat_run_collect: need at least argv0 + database + query");

    /* Writable copies (optionInit mutates argv in place). */
    char **argv = (char **) R_alloc((size_t) argc, sizeof(char *));
    for (int i = 0; i < argc; i++) {
        const char *s = CHAR(STRING_ELT(r_argv, i));
        size_t n = strlen(s) + 1;
        char *buf = R_alloc(n, 1);
        memcpy(buf, s, n);
        argv[i] = buf;
    }

    struct rblatResult *result = NULL;
    char *errMsg = NULL, *warnMsg = NULL;
    int ok = rblat_run_argv_collect(argc, argv, &result, &errMsg, &warnMsg);
    if (!ok) {
        char msg[1024];
        snprintf(msg, sizeof(msg), "blat error: %s",
                 errMsg ? errMsg : "(no message)");
        if (errMsg) free(errMsg);
        Rf_error("%s", msg);          /* does not return */
    }
    rblat_emit_warning(warnMsg);

    SEXP df = PROTECT(rblat_result_to_df(result));
    rblat_result_free(result);
    UNPROTECT(1);
    return df;
}

/* View an R character vector as a C char* array (borrowing R's string pool --
 * valid for the duration of the .Call; the engine cloneString's what it keeps).
 * Returns a 1-slot dummy for length 0 so the pointer is never NULL. */
static char **rblat_strvec(SEXP v, int *n)
{
    *n = LENGTH(v);
    char **a = (char **) R_alloc((size_t)(*n > 0 ? *n : 1), sizeof(char *));
    for (int i = 0; i < *n; i++)
        a[i] = (char *) CHAR(STRING_ELT(v, i));
    return a;
}

/* Run blat in-process with in-memory database/query sequences, returning a PSL
 * data.frame directly (Phase 4) -- no input FASTA, no output file. r_argv holds
 * argv0 + option flags + two placeholder positional args. */
SEXP rblat_run_collect_mem(SEXP r_argv, SEXP r_dbNames, SEXP r_dbSeqs,
                           SEXP r_qNames, SEXP r_qSeqs)
{
    if (TYPEOF(r_argv) != STRSXP)
        Rf_error("rblat_run_collect_mem: argv must be a character vector");
    if (TYPEOF(r_dbNames) != STRSXP || TYPEOF(r_dbSeqs) != STRSXP ||
        TYPEOF(r_qNames)  != STRSXP || TYPEOF(r_qSeqs)  != STRSXP)
        Rf_error("rblat_run_collect_mem: names/sequences must be character vectors");
    if (LENGTH(r_dbNames) != LENGTH(r_dbSeqs) ||
        LENGTH(r_qNames)  != LENGTH(r_qSeqs))
        Rf_error("rblat_run_collect_mem: names and sequences differ in length");

    int argc = LENGTH(r_argv);
    char **argv = (char **) R_alloc((size_t) argc, sizeof(char *));
    for (int i = 0; i < argc; i++) {
        const char *s = CHAR(STRING_ELT(r_argv, i));
        size_t n = strlen(s) + 1;
        char *buf = R_alloc(n, 1);
        memcpy(buf, s, n);
        argv[i] = buf;
    }

    int nDb, nDbS, nQ, nQS;
    char **dbNames = rblat_strvec(r_dbNames, &nDb);
    char **dbSeqs  = rblat_strvec(r_dbSeqs,  &nDbS);
    char **qNames  = rblat_strvec(r_qNames,  &nQ);
    char **qSeqs   = rblat_strvec(r_qSeqs,   &nQS);

    struct rblatResult *result = NULL;
    char *errMsg = NULL, *warnMsg = NULL;
    int ok = rblat_run_argv_collect_mem(argc, argv, dbNames, dbSeqs, nDb,
                                        qNames, qSeqs, nQ, &result,
                                        &errMsg, &warnMsg);
    if (!ok) {
        char msg[1024];
        snprintf(msg, sizeof(msg), "blat error: %s",
                 errMsg ? errMsg : "(no message)");
        if (errMsg) free(errMsg);
        Rf_error("%s", msg);          /* does not return */
    }
    rblat_emit_warning(warnMsg);

    SEXP df = PROTECT(rblat_result_to_df(result));
    rblat_result_free(result);
    UNPROTECT(1);
    return df;
}

static const R_CallMethodDef CallEntries[] = {
    {"C_rblat_version",         (DL_FUNC) &rblat_version,         0},
    {"C_rblat_run",             (DL_FUNC) &rblat_run,             1},
    {"C_rblat_run_collect",     (DL_FUNC) &rblat_run_collect,     1},
    {"C_rblat_run_collect_mem", (DL_FUNC) &rblat_run_collect_mem, 5},
    {NULL, NULL, 0}
};

void R_init_rblat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
