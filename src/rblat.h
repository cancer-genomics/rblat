/* rblat internal C API shared between the R-facing layer and the engine adapter.
 * Deliberately free of both R and Kent headers so it can be included by either. */
#ifndef RBLAT_H
#define RBLAT_H

/* Run blat with the given argv (argv[0] is the program name, e.g. "blat").
 * Returns 1 on success, 0 if an errAbort/exit was caught. On failure,
 * *errMsgOut receives a malloc'd message string (caller frees), else NULL. */
int rblat_run_argv(int argc, char **argv, char **errMsgOut);

/* Suppress blat's status printf's ("Loaded ...", "Searched ...") when nonzero. */
extern int rblat_quiet;

/* --- Phase 3: in-memory PSL collection ----------------------------------- *
 * Alignments are appended directly to growable C arrays instead of being
 * fprintf'd to a .psl file, then handed to R as a data.frame -- no PSL text,
 * no output file.  This struct holds only standard C types so it can be read
 * on the R-facing side (rblat_init.c) without pulling in the Kent headers that
 * clash with Rinternals.h.  One entry per array slot; fields are in the same
 * order as blat_colnames (the 21-column PSL layout). */
struct rblatResult
    {
    int nrow;                   /* Rows filled so far. */
    int cap;                    /* Allocated capacity of each array. */
    /* Integer columns. */
    int *match;                 /* Matching bases (non-repeat). */
    int *misMatch;              /* Mismatching bases. */
    int *repMatch;              /* Matching bases in repeat regions. */
    int *nCount;                /* Count of 'N' bases. */
    int *qNumInsert;            /* Query gap count. */
    int *qBaseInsert;           /* Query gap bases. */
    int *tNumInsert;            /* Target gap count. */
    int *tBaseInsert;           /* Target gap bases. */
    int *qSize;                 /* Query sequence size. */
    int *qStart;                /* Query alignment start. */
    int *qEnd;                  /* Query alignment end. */
    int *tSize;                 /* Target sequence size. */
    int *tStart;                /* Target alignment start. */
    int *tEnd;                  /* Target alignment end. */
    int *blockCount;            /* Number of aligned blocks. */
    /* Character columns (each element malloc'd; freed by rblat_result_free). */
    char **strand;              /* "+"/"-" (or two chars if target strand shown). */
    char **qName;               /* Query sequence name. */
    char **tName;               /* Target sequence name. */
    char **blockSizes;          /* Comma-terminated block sizes. */
    char **qStarts;             /* Comma-terminated query block starts. */
    char **tStarts;             /* Comma-terminated target block starts. */
    };

/* Arm/disarm in-memory collection for the next rblat_blat_main() call. When set
 * to a non-NULL collector, blat() appends to it and opens no output file; pass
 * NULL to restore normal file-based output. Defined in rblat_blat.c. */
void rblat_blat_set_collector(struct rblatResult *r);

/* Allocate an empty result collector. */
struct rblatResult *rblat_result_new(void);

/* Free a result collector and every string it owns. Safe on NULL. */
void rblat_result_free(struct rblatResult *r);

/* Run blat with in-memory collection instead of file output. Returns 1 on
 * success (with *outResult pointing to a collector the caller frees via
 * rblat_result_free), 0 if an errAbort/exit was caught (then *outResult is NULL
 * and *errMsgOut is a malloc'd message the caller frees with free()).
 * On success, any blat warn() output captured during the run is returned in
 * *warnMsgOut (malloc'd, caller frees) or NULL if there were no warnings. */
int rblat_run_argv_collect(int argc, char **argv,
                           struct rblatResult **outResult,
                           char **errMsgOut, char **warnMsgOut);

/* Like rblat_run_argv_collect, but the database and query sequences are passed
 * in memory as parallel (name, sequence) C-string arrays -- no input FASTA
 * files. argv still carries the option flags (and two placeholder positional
 * args). Same return/ownership contract as rblat_run_argv_collect. */
int rblat_run_argv_collect_mem(int argc, char **argv,
                               char **dbNames, char **dbSeqs, int nDb,
                               char **qNames, char **qSeqs, int nQ,
                               struct rblatResult **outResult,
                               char **errMsgOut, char **warnMsgOut);

#endif /* RBLAT_H */
