/* rblat_input.c - build in-memory sequence lists for rblat (Phase 4).
 *
 * Turns R-supplied (name, sequence) pairs into a Kent dnaSeq list so a blat run
 * needs no input FASTA file. The list is fed to blat() as the database, and to
 * searchOneIndex() as the query, in place of the file-reading paths.
 *
 * Sequences are copied verbatim (case preserved) -- exactly what
 * faMixedSpeedReadNext() would have handed back from a .fa file. Per-run case
 * conversion (faToDna/faToProtein) then happens downstream, identical to the
 * file path.
 *
 * Kent headers only; no R headers (the char** arrays cross the boundary from
 * rblat_init.c as plain C).
 */
#include "common.h"
#include "dnaseq.h"
#include "rblat.h"

struct dnaSeq *rblat_seqs_from_arrays(char **names, char **seqs, int n)
/* Build a dnaSeq list from parallel name/sequence arrays, preserving order.
 * Each element owns cloned copies of its name and sequence bytes; free the
 * whole list with freeDnaSeqList(). */
{
struct dnaSeq *list = NULL;
int i;
for (i = 0; i < n; ++i)
    {
    char *dna = cloneString(seqs[i]);   /* mutable; faToDna edits in place */
    struct dnaSeq *seq = newDnaSeq(dna, (int) strlen(dna), names[i]);
    slAddHead(&list, seq);
    }
slReverse(&list);
return list;
}
