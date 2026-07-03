/* rblat_collect.c - in-memory PSL collection for rblat (Phase 3).
 *
 * Provides a gfOutput handler that, instead of fprintf'ing PSL records to a
 * file, appends the 21 PSL fields to growable C arrays (struct rblatResult).
 * The R-facing layer then turns those arrays into a data.frame directly, so a
 * blat run produces no PSL text and no output file.
 *
 * The alignment-analysis math here is a faithful copy of savePslx() in gfOut.c
 * (match / mismatch / repMatch / N / insert counting, identity filtering, and
 * strand coordinate flipping); only the final "emit a record" step differs --
 * where savePslx does fprintf, we append to arrays.  Keeping the arithmetic
 * identical is what guarantees rblat's in-memory output matches the .psl a
 * file-based blat run would have produced.
 *
 * This translation unit includes the Kent headers but NOT R's headers.
 */
#include "common.h"
#include "hash.h"
#include "dystring.h"
#include "dnautil.h"
#include "trans3.h"
#include "fuzzyFind.h"
#include "genoFind.h"
#include "bits.h"
#include "rblat.h"

/* gfOutput.data payload: just a back-pointer to the collector. Allocated with
 * a Kent allocator so gfOutputFree()'s freeMem(out->data) can release it. */
struct rblatCollectData
    {
    struct rblatResult *result;
    };

struct rblatResult *rblat_result_new(void)
/* Allocate an empty result collector (all arrays NULL until first append). */
{
struct rblatResult *r;
AllocVar(r);
r->nrow = 0;
r->cap = 0;
return r;
}

void rblat_result_free(struct rblatResult *r)
/* Free a result collector and every string it owns. */
{
int i;
if (r == NULL)
    return;
for (i = 0; i < r->nrow; ++i)
    {
    freeMem(r->strand[i]);
    freeMem(r->qName[i]);
    freeMem(r->tName[i]);
    freeMem(r->blockSizes[i]);
    freeMem(r->qStarts[i]);
    freeMem(r->tStarts[i]);
    }
freeMem(r->match);       freeMem(r->misMatch);    freeMem(r->repMatch);
freeMem(r->nCount);      freeMem(r->qNumInsert);  freeMem(r->qBaseInsert);
freeMem(r->tNumInsert);  freeMem(r->tBaseInsert); freeMem(r->qSize);
freeMem(r->qStart);      freeMem(r->qEnd);        freeMem(r->tSize);
freeMem(r->tStart);      freeMem(r->tEnd);        freeMem(r->blockCount);
freeMem(r->strand);      freeMem(r->qName);       freeMem(r->tName);
freeMem(r->blockSizes);  freeMem(r->qStarts);     freeMem(r->tStarts);
freeMem(r);
}

static void rblatResultGrow(struct rblatResult *r)
/* Ensure room for at least one more row, growing all arrays together. */
{
if (r->nrow < r->cap)
    return;
int newCap = (r->cap == 0) ? 16 : r->cap * 2;
r->match       = needMoreMem(r->match,       r->cap * sizeof(int),  newCap * sizeof(int));
r->misMatch    = needMoreMem(r->misMatch,    r->cap * sizeof(int),  newCap * sizeof(int));
r->repMatch    = needMoreMem(r->repMatch,    r->cap * sizeof(int),  newCap * sizeof(int));
r->nCount      = needMoreMem(r->nCount,      r->cap * sizeof(int),  newCap * sizeof(int));
r->qNumInsert  = needMoreMem(r->qNumInsert,  r->cap * sizeof(int),  newCap * sizeof(int));
r->qBaseInsert = needMoreMem(r->qBaseInsert, r->cap * sizeof(int),  newCap * sizeof(int));
r->tNumInsert  = needMoreMem(r->tNumInsert,  r->cap * sizeof(int),  newCap * sizeof(int));
r->tBaseInsert = needMoreMem(r->tBaseInsert, r->cap * sizeof(int),  newCap * sizeof(int));
r->qSize       = needMoreMem(r->qSize,       r->cap * sizeof(int),  newCap * sizeof(int));
r->qStart      = needMoreMem(r->qStart,      r->cap * sizeof(int),  newCap * sizeof(int));
r->qEnd        = needMoreMem(r->qEnd,        r->cap * sizeof(int),  newCap * sizeof(int));
r->tSize       = needMoreMem(r->tSize,       r->cap * sizeof(int),  newCap * sizeof(int));
r->tStart      = needMoreMem(r->tStart,      r->cap * sizeof(int),  newCap * sizeof(int));
r->tEnd        = needMoreMem(r->tEnd,        r->cap * sizeof(int),  newCap * sizeof(int));
r->blockCount  = needMoreMem(r->blockCount,  r->cap * sizeof(int),  newCap * sizeof(int));
r->strand      = needMoreMem(r->strand,      r->cap * sizeof(char *), newCap * sizeof(char *));
r->qName       = needMoreMem(r->qName,       r->cap * sizeof(char *), newCap * sizeof(char *));
r->tName       = needMoreMem(r->tName,       r->cap * sizeof(char *), newCap * sizeof(char *));
r->blockSizes  = needMoreMem(r->blockSizes,  r->cap * sizeof(char *), newCap * sizeof(char *));
r->qStarts     = needMoreMem(r->qStarts,     r->cap * sizeof(char *), newCap * sizeof(char *));
r->tStarts     = needMoreMem(r->tStarts,     r->cap * sizeof(char *), newCap * sizeof(char *));
r->cap = newCap;
}

static void collectPslx(char *chromName, int chromSize, int chromOffset,
	struct ffAli *ali, struct dnaSeq *tSeq, struct dnaSeq *qSeq,
	boolean isRc, enum ffStringency stringency, int minMatch,
	struct rblatResult *result, struct hash *t3Hash,
	boolean reportTargetStrand, boolean targetIsRc,
	struct hash *maskHash, int minIdentity)
/* Analyse one alignment and, if it passes the identity filter, append its PSL
 * fields to the collector. Counting logic mirrors savePslx() in gfOut.c. */
{
struct ffAli *ff, *nextFf;
struct ffAli *right = ffRightmost(ali);
DNA *needle = qSeq->dna;
DNA *hay = tSeq->dna;
int nStart = ali->nStart - needle;
int nEnd = right->nEnd - needle;
int hStart, hEnd;
int nInsertBaseCount = 0;
int nInsertCount = 0;
int hInsertBaseCount = 0;
int hInsertCount = 0;
int matchCount = 0;
int mismatchCount = 0;
int repMatch = 0;
int countNs = 0;
DNA *np, *hp, n, h;
int blockSize;
int i;
struct trans3 *t3List = NULL;
Bits *maskBits = NULL;

if (maskHash != NULL)
    maskBits = hashMustFindVal(maskHash, tSeq->name);
if (t3Hash != NULL)
    t3List = hashMustFindVal(t3Hash, tSeq->name);
hStart = trans3GenoPos(ali->hStart, tSeq, t3List, FALSE) + chromOffset;
hEnd = trans3GenoPos(right->hEnd, tSeq, t3List, TRUE) + chromOffset;

/* Count up matches, mismatches, inserts, etc. */
for (ff = ali; ff != NULL; ff = nextFf)
    {
    nextFf = ff->right;
    blockSize = ff->nEnd - ff->nStart;
    np = ff->nStart;
    hp = ff->hStart;
    for (i=0; i<blockSize; ++i)
	{
	n = np[i];
	h = hp[i];
	if (n == 'n' || h == 'n')
	    ++countNs;
	else
	    {
	    if (n == h)
		{
		if (maskBits != NULL)
		    {
		    int seqOff = hp + i - hay;
		    if (bitReadOne(maskBits, seqOff))
			++repMatch;
		    else
			++matchCount;
		    }
		else
		    ++matchCount;
		}
	    else
		++mismatchCount;
	    }
	}
    if (nextFf != NULL)
	{
	int nhStart = trans3GenoPos(nextFf->hStart, tSeq, t3List, FALSE) + chromOffset;
	int ohEnd = trans3GenoPos(ff->hEnd, tSeq, t3List, TRUE) + chromOffset;
	int hGap = nhStart - ohEnd;
	int nGap = nextFf->nStart - ff->nEnd;

	if (nGap != 0)
	    {
	    ++nInsertCount;
	    nInsertBaseCount += nGap;
	    }
	if (hGap != 0)
	    {
	    ++hInsertCount;
	    hInsertBaseCount += hGap;
	    }
	}
    }

/* See if it looks good enough to output, and if so append it. */
    {
    int gaps = nInsertCount + (stringency == ffCdna ? 0: hInsertCount);
    if ((matchCount + repMatch + mismatchCount) > 0)
	{
	int id = roundingScale(1000, matchCount + repMatch - 2*gaps, matchCount + repMatch + mismatchCount);
	if (id >= minIdentity)
	    {
	    struct dyString *blockSizes = newDyString(256);
	    struct dyString *qStarts = newDyString(256);
	    struct dyString *tStarts = newDyString(256);
	    char strandBuf[3];
	    int si = 0;
	    int row;

	    if (isRc)
		{
		int temp;
		int oSize = qSeq->size;
		temp = nStart;
		nStart = oSize - nEnd;
		nEnd = oSize - temp;
		}
	    if (targetIsRc)
		{
		int temp;
		temp = hStart;
		hStart = chromSize - hEnd;
		hEnd = chromSize - temp;
		}

	    strandBuf[si++] = (isRc ? '-' : '+');
	    if (reportTargetStrand)
		strandBuf[si++] = (targetIsRc ? '-' : '+');
	    strandBuf[si] = '\0';

	    /* Build the three comma-terminated block strings exactly as
	     * savePslx would have fprintf'd them ("%ld," / "%d,"). */
	    for (ff = ali; ff != NULL; ff = ff->right)
		dyStringPrintf(blockSizes, "%ld,", (long)(ff->nEnd - ff->nStart));
	    for (ff = ali; ff != NULL; ff = ff->right)
		dyStringPrintf(qStarts, "%ld,", (long)(ff->nStart - needle));
	    for (ff = ali; ff != NULL; ff = ff->right)
		dyStringPrintf(tStarts, "%d,", trans3GenoPos(ff->hStart, tSeq, t3List, FALSE) + chromOffset);

	    rblatResultGrow(result);
	    row = result->nrow;
	    result->match[row]       = matchCount;
	    result->misMatch[row]    = mismatchCount;
	    result->repMatch[row]    = repMatch;
	    result->nCount[row]      = countNs;
	    result->qNumInsert[row]  = nInsertCount;
	    result->qBaseInsert[row] = nInsertBaseCount;
	    result->tNumInsert[row]  = hInsertCount;
	    result->tBaseInsert[row] = hInsertBaseCount;
	    result->qSize[row]       = qSeq->size;
	    result->qStart[row]      = nStart;
	    result->qEnd[row]        = nEnd;
	    result->tSize[row]       = chromSize;
	    result->tStart[row]      = hStart;
	    result->tEnd[row]        = hEnd;
	    result->blockCount[row]  = ffAliCount(ali);
	    result->strand[row]      = cloneString(strandBuf);
	    result->qName[row]       = cloneString(qSeq->name);
	    result->tName[row]       = cloneString(chromName);
	    result->blockSizes[row]  = cloneString(blockSizes->string);
	    result->qStarts[row]     = cloneString(qStarts->string);
	    result->tStarts[row]     = cloneString(tStarts->string);
	    result->nrow += 1;

	    dyStringFree(&blockSizes);
	    dyStringFree(&qStarts);
	    dyStringFree(&tStarts);
	    }
	}
    }
}

static void rblatCollectOut(char *chromName, int chromSize, int chromOffset,
	struct ffAli *ali, struct dnaSeq *tSeq, struct hash *t3Hash,
	struct dnaSeq *qSeq, boolean qIsRc, boolean tIsRc,
	enum ffStringency stringency, int minMatch, struct gfOutput *out)
/* gfOutput.out callback: collect one alignment into the result arrays. */
{
struct rblatCollectData *cd = out->data;
collectPslx(chromName, chromSize, chromOffset, ali, tSeq, qSeq,
    qIsRc, stringency, minMatch, cd->result, t3Hash,
    out->reportTargetStrand, tIsRc, out->maskHash, out->minGood);
}

struct gfOutput *rblatOutputCollector(int goodPpt, boolean qIsProt,
	boolean tIsProt, struct rblatResult *result)
/* Set up a gfOutput that appends alignments to an in-memory collector instead
 * of writing a .psl file. queryOut/fileHead stay NULL (no per-query or header
 * output). Mirrors gfOutputPsl() setup minus the file. */
{
struct gfOutput *out;
struct rblatCollectData *cd;
AllocVar(out);
out->minGood = goodPpt;
out->qIsProt = qIsProt;
out->tIsProt = tIsProt;
AllocVar(cd);
cd->result = result;
out->out = rblatCollectOut;
out->data = cd;
return out;
}
