/* rblat engine adapter.
 *
 * Bridges the R-facing layer (rblat_init.c, which includes R headers) and the
 * BLAT engine (rblat_blat.c). This translation unit includes the Kent headers
 * but NOT R's headers, because common.h redefines TRUE/FALSE/boolean and clashes
 * with Rinternals.h. The R layer talks to us only through the plain-C prototype
 * in rblat.h.
 */
#include "common.h"
#include "errabort.h"
#include "errCatch.h"
#include "dnaseq.h"
#include "rblat.h"

extern int rblat_blat_main(int argc, char *argv[]);

/* From rblat_input.c / rblat_blat.c (Kent-typed, so not in rblat.h). */
extern struct dnaSeq *rblat_seqs_from_arrays(char **names, char **seqs, int n);
extern void rblat_blat_set_inputs(struct dnaSeq *db, struct dnaSeq *query,
                                  boolean enable);

static char *rblat_msg_dup(struct dyString *msg)
/* Return a plain malloc'd copy of a captured errCatch message (so the R side
 * can free() it after errCatchFree releases the dyString), or NULL if there is
 * nothing to report. */
{
if (msg == NULL || msg->string == NULL || msg->string[0] == '\0')
    return NULL;
size_t n = strlen(msg->string) + 1;
char *dst = malloc(n);
if (dst != NULL)
    memcpy(dst, msg->string, n);
return dst;
}

int rblat_run_argv(int argc, char **argv, char **errMsgOut)
/* Run blat's main() under an errCatch guard. Returns 1 on success, 0 on a
 * caught errAbort/exit; on failure *errMsgOut points to a malloc'd copy of the
 * message (caller frees with free()), NULL otherwise. */
{
struct errCatch *errCatch = errCatchNew();
*errMsgOut = NULL;
if (errCatchStart(errCatch))
    {
    rblat_blat_main(argc, argv);
    }
errCatchEnd(errCatch);
boolean gotError = errCatch->gotError;
if (gotError && errCatch->message != NULL && errCatch->message->string != NULL)
    {
    /* Copy out with plain malloc so the R side can free() it after errCatchFree
     * releases the dyString. */
    char *src = errCatch->message->string;
    size_t n = strlen(src) + 1;
    char *dst = malloc(n);
    if (dst != NULL)
        memcpy(dst, src, n);
    *errMsgOut = dst;
    }
errCatchFree(&errCatch);
return gotError ? 0 : 1;
}

int rblat_run_argv_collect(int argc, char **argv,
                           struct rblatResult **outResult,
                           char **errMsgOut, char **warnMsgOut)
/* Run blat's main() with in-memory collection under an errCatch guard. On
 * success returns 1 and *outResult owns the collected rows (and *warnMsgOut
 * holds any captured warn() text or NULL); on a caught errAbort/exit returns 0,
 * *outResult is NULL, and *errMsgOut is a malloc'd message (caller frees). */
{
struct rblatResult *result = rblat_result_new();
struct errCatch *errCatch = errCatchNew();
*errMsgOut = NULL;
*warnMsgOut = NULL;
*outResult = NULL;
rblat_blat_set_collector(result);
if (errCatchStart(errCatch))
    {
    rblat_blat_main(argc, argv);
    }
errCatchEnd(errCatch);
rblat_blat_set_collector(NULL);   /* disarm even if we longjmp'd out */
boolean gotError = errCatch->gotError;
if (gotError)
    {
    *errMsgOut = rblat_msg_dup(errCatch->message);
    rblat_result_free(result);
    }
else
    {
    /* On success the message holds any warn() output captured during the run. */
    *warnMsgOut = rblat_msg_dup(errCatch->message);
    *outResult = result;
    }
errCatchFree(&errCatch);
return gotError ? 0 : 1;
}

int rblat_run_argv_collect_mem(int argc, char **argv,
                               char **dbNames, char **dbSeqs, int nDb,
                               char **qNames, char **qSeqs, int nQ,
                               struct rblatResult **outResult,
                               char **errMsgOut, char **warnMsgOut)
/* Run blat with in-memory database/query sequences and in-memory collection.
 * Builds Kent dnaSeq lists from the (name, sequence) arrays, runs the engine
 * under an errCatch guard, and frees the lists. Same contract as
 * rblat_run_argv_collect. */
{
struct rblatResult *result = rblat_result_new();
struct dnaSeq *dbList = rblat_seqs_from_arrays(dbNames, dbSeqs, nDb);
struct dnaSeq *qList = rblat_seqs_from_arrays(qNames, qSeqs, nQ);
struct errCatch *errCatch = errCatchNew();
*errMsgOut = NULL;
*warnMsgOut = NULL;
*outResult = NULL;
rblat_blat_set_collector(result);
rblat_blat_set_inputs(dbList, qList, TRUE);
if (errCatchStart(errCatch))
    {
    rblat_blat_main(argc, argv);
    }
errCatchEnd(errCatch);
rblat_blat_set_inputs(NULL, NULL, FALSE);   /* disarm even on longjmp */
rblat_blat_set_collector(NULL);
boolean gotError = errCatch->gotError;
if (gotError)
    {
    *errMsgOut = rblat_msg_dup(errCatch->message);
    rblat_result_free(result);
    }
else
    {
    *warnMsgOut = rblat_msg_dup(errCatch->message);
    *outResult = result;
    }
errCatchFree(&errCatch);
/* We own these lists (blat() does not free them in memory mode). */
freeDnaSeqList(&dbList);
freeDnaSeqList(&qList);
return gotError ? 0 : 1;
}
