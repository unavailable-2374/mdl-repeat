/*
 * discover_mask.c — Genome masking subsystem for the discover engine.
 *
 * After a candidate family is extended, its instances must be masked
 * out of the genome to prevent the next seed iteration from re-finding
 * the same family.  Masking uses a 1-vs-1 banded DP alignment between
 * each putative occurrence and the consensus "master" sequence to
 * determine the precise extent to mask.
 *
 * This module is a faithful port of the RepeatScout masking routines
 * (lines 1387-1578 + 1207-1364 in the original RepeatScout source),
 * extracted from discover.c by M3#8 to keep that file under 2000
 * lines.  All shared state lives in DiscoverContext (see
 * discover_internal.h).
 *
 * Public entry point:
 *   mask_headptr(C, headptr) — mask all unmasked occurrences of the
 *                              family currently in masters[R-1].
 */

#include <stdio.h>
#include <stdlib.h>

#include "discover_internal.h"

/* ================================================================
 * Mask scoring: 1-vs-1 banded DP cell
 * (Exact port of RepeatScout lines 2071-2161)
 * ================================================================ */

static int compute_maskscore_right(DiscoverContext *C, char mbase, int repeaty,
                                   gpos_t seqpos, int offset)
{
    int oldoffset, tempscore, ans, ismatch, x;

    ans = -1000000000;

    /* Case A: gap in sequence */
    if (offset < C->MAXOFFSET) {
        oldoffset = offset + 1;
        tempscore = C->maskscore[(repeaty + 1) % 2][oldoffset + C->MAXOFFSET] + C->GAP;
        if (tempscore > ans) ans = tempscore;
    }

    /* Case B: diagonal */
    oldoffset = offset;
    tempscore = C->maskscore[(repeaty + 1) % 2][oldoffset + C->MAXOFFSET];
    if (mbase == C->sequence[seqpos + C->l + repeaty + offset])
        tempscore += C->MATCH;
    else
        tempscore += C->MISMATCH;
    if (tempscore > ans) ans = tempscore;

    /* Case C: multiple gaps */
    for (oldoffset = -C->MAXOFFSET; oldoffset < offset; oldoffset++) {
        ismatch = 0;
        for (x = oldoffset; x <= offset; x++) {
            if (mbase == C->sequence[seqpos + C->l + repeaty + x])
                ismatch = 1;
        }
        tempscore = C->maskscore[(repeaty + 1) % 2][oldoffset + C->MAXOFFSET];
        tempscore += (offset - oldoffset) * C->GAP;
        if (ismatch) tempscore += C->MATCH;
        else tempscore += C->MISMATCH;
        if (tempscore > ans) ans = tempscore;
    }

    return ans;
}

static int compute_maskscore_left(DiscoverContext *C, char mbase, int repeatw,
                                  gpos_t seqpos, int offset)
{
    int oldoffset, tempscore, ans, ismatch, x;

    ans = -1000000000;

    /* Case A: gap in sequence */
    if (offset > -C->MAXOFFSET) {
        oldoffset = offset - 1;
        tempscore = C->maskscore[(repeatw + 1) % 2][oldoffset + C->MAXOFFSET] + C->GAP;
        if (tempscore > ans) ans = tempscore;
    }

    /* Case B: diagonal */
    oldoffset = offset;
    tempscore = C->maskscore[(repeatw + 1) % 2][oldoffset + C->MAXOFFSET];
    if (mbase == C->sequence[seqpos - repeatw - 1 + offset])
        tempscore += C->MATCH;
    else
        tempscore += C->MISMATCH;
    if (tempscore > ans) ans = tempscore;

    /* Case C: multiple gaps */
    for (oldoffset = offset + 1; oldoffset <= C->MAXOFFSET; oldoffset++) {
        ismatch = 0;
        for (x = offset; x <= oldoffset; x++) {
            if (mbase == C->sequence[seqpos - repeatw - 1 + x])
                ismatch = 1;
        }
        tempscore = C->maskscore[(repeatw + 1) % 2][oldoffset + C->MAXOFFSET];
        tempscore += (oldoffset - offset) * C->GAP;
        if (ismatch) tempscore += C->MATCH;
        else tempscore += C->MISMATCH;
        if (tempscore > ans) ans = tempscore;
    }

    return ans;
}

/* ================================================================
 * Per-occurrence banded DP extension to determine mask extent
 * (Port of RepeatScout lines 1387-1578)
 * ================================================================ */

static gpos_t maskextend_right(DiscoverContext *C, int rc, gpos_t seqpos,
                               int modelpos)
{
    int offset, bestmaskscore_val, bestbestmaskscore, bestoffset, bestbestoffset;
    gpos_t bEnd;
    int x, bestx, maxext;
    char mbase;

    /* Initialize maskscore */
    for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
        C->maskscore[1][offset + C->MAXOFFSET] = 0;
        if (offset < 0) C->maskscore[1][offset + C->MAXOFFSET] += -offset * C->GAP;
        if (offset > 0) C->maskscore[1][offset + C->MAXOFFSET] += offset * C->GAP;
    }

    bestbestmaskscore = 0;
    bestbestoffset = 0;
    bestoffset = 0;
    bestx = 0;

    /* Find sequence boundary */
    bEnd = 0;
    x = 0;
    while (C->disc_boundaries[x] > 0) {
        if (seqpos < C->disc_boundaries[x] + PADLENGTH) {
            bEnd = C->disc_boundaries[x] + PADLENGTH;
            break;
        }
        x++;
    }

    /* Max extension distance in model */
    maxext = C->masterend[C->R - 1] - (modelpos + C->l + 1);
    if (rc)
        maxext = modelpos - (C->masterstart[C->R - 1] + 1);

    for (x = 0; x < maxext; x++) {
        if (bEnd > 0 && seqpos + C->l + x == bEnd)
            break;

        /* Get consensus base (4 combinations of direction x strand) */
        if (rc)
            mbase = compl_base(C->masters[C->R - 1][modelpos - x - 1]);
        else
            mbase = C->masters[C->R - 1][modelpos + C->l + x];

        bestmaskscore_val = -1000000000;
        for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
            if (bEnd > 0 && seqpos + C->l + x + offset >= bEnd) {
                C->maskscore[(x + 2) % 2][offset + C->MAXOFFSET] = 0;
            } else {
                C->maskscore[(x + 2) % 2][offset + C->MAXOFFSET] =
                    compute_maskscore_right(C, mbase, x, seqpos, offset);
            }
            if (C->maskscore[(x + 2) % 2][offset + C->MAXOFFSET] > bestmaskscore_val) {
                bestmaskscore_val = C->maskscore[(x + 2) % 2][offset + C->MAXOFFSET];
                bestoffset = offset;
            }
        }
        if (bestmaskscore_val > bestbestmaskscore) {
            bestbestmaskscore = bestmaskscore_val;
            bestbestoffset = bestoffset;
            bestx = x;
        }
        if (x - bestx >= C->WHEN_TO_STOP) break;
    }

    return seqpos + C->l + bestx + bestbestoffset;
}

static gpos_t maskextend_left(DiscoverContext *C, int rc, gpos_t seqpos,
                              int modelpos)
{
    int offset, bestmaskscore_val, bestbestmaskscore, bestoffset, bestbestoffset;
    int maxext;
    gpos_t bStart;
    int x, bestx;
    char mbase;

    /* Initialize maskscore */
    for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
        C->maskscore[1][offset + C->MAXOFFSET] = 0;
        if (offset < 0) C->maskscore[1][offset + C->MAXOFFSET] += -offset * C->GAP;
        if (offset > 0) C->maskscore[1][offset + C->MAXOFFSET] += offset * C->GAP;
    }

    bestbestoffset = 0;
    bestoffset = 0;
    bestbestmaskscore = 0;
    bestx = 0;

    /* Find sequence boundary */
    bStart = 0;
    x = 0;
    while (C->disc_boundaries[x] > 0) {
        if (seqpos < C->disc_boundaries[x] + PADLENGTH) {
            if (x > 0)
                bStart = C->disc_boundaries[x - 1] + PADLENGTH;
            break;
        }
        x++;
    }

    /* Max extension distance in model */
    maxext = modelpos - C->masterstart[C->R - 1];
    if (rc)
        maxext = C->masterend[C->R - 1] - (modelpos + C->l - 1);

    for (x = 0; x < maxext; x++) {
        if (seqpos - x - 1 == bStart)
            break;

        /* Get consensus base */
        if (rc)
            mbase = compl_base(C->masters[C->R - 1][modelpos + C->l + x]);
        else
            mbase = C->masters[C->R - 1][modelpos - x - 1];

        bestmaskscore_val = -1000000000;
        for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
            if (seqpos - x - 1 + offset < bStart) {
                C->maskscore[(x + 2) % 2][offset + C->MAXOFFSET] = 0;
            } else {
                C->maskscore[(x + 2) % 2][offset + C->MAXOFFSET] =
                    compute_maskscore_left(C, mbase, x, seqpos, offset);
            }
            if (C->maskscore[(x + 2) % 2][offset + C->MAXOFFSET] > bestmaskscore_val) {
                bestmaskscore_val = C->maskscore[(x + 2) % 2][offset + C->MAXOFFSET];
                bestoffset = offset;
            }
        }
        if (bestmaskscore_val > bestbestmaskscore) {
            bestbestmaskscore = bestmaskscore_val;
            bestbestoffset = bestoffset;
            bestx = x;
        }
        if (x - bestx >= C->WHEN_TO_STOP) break;
    }

    return seqpos - bestx + bestbestoffset - 1;
}

/* ================================================================
 * mask_headptr — public entry: mask the family in masters[R-1]
 * (Port of RepeatScout lines 1207-1364, with dynamic value[])
 * ================================================================ */

void mask_headptr(DiscoverContext *C, struct llist **headptr)
{
    int smallh, h, h2;
    gpos_t x, y;
    gpos_t startmask, endmask;
    gpos_t bestsequencey_val, bestsequencew_val;
    struct repeatllist **repeatheadptr;
    struct repeatllist *repeattmp, *nextrepeattmp;
    struct llist *tmp, *tmp2;
    struct posllist *postmp;
    int rrev;

    /* Build small hash table of l-mers from consensus */
    repeatheadptr = malloc(SMALLHASH_SIZE * sizeof(*repeatheadptr));
    if (repeatheadptr == NULL) {
        fprintf(stderr, "discover: out of memory for repeatheadptr\n"); exit(1);
    }
    for (smallh = 0; smallh < SMALLHASH_SIZE; smallh++)
        repeatheadptr[smallh] = NULL;

    for (x = C->masterstart[C->R - 1]; x <= C->masterend[C->R - 1] - C->l + 1; x++) {
        smallh = smallhash_function(C->masters[C->R - 1] + x);
        if (smallh < 0) continue;

        /* Check if already present */
        repeattmp = repeatheadptr[smallh];
        while (repeattmp != NULL) {
            if (lmermatch(C, repeattmp->value, C->masters[C->R - 1] + x))
                break;
            repeattmp = repeattmp->next;
        }
        if (repeattmp != NULL) continue;

        /* Add new entry (FIXED: dynamic value allocation) */
        repeattmp = malloc(sizeof(*repeattmp));
        if (repeattmp == NULL) {
            fprintf(stderr, "discover: out of memory\n"); exit(1);
        }
        repeattmp->value = malloc(C->l * sizeof(char));
        if (repeattmp->value == NULL) {
            fprintf(stderr, "discover: out of memory for value\n"); exit(1);
        }
        for (y = 0; y < C->l; y++)
            repeattmp->value[y] = C->masters[C->R - 1][x + y];
        repeattmp->repeatpos = x;
        repeattmp->next = repeatheadptr[smallh];
        repeatheadptr[smallh] = repeattmp;
    }

    /* Find and extend hits */
    for (smallh = 0; smallh < SMALLHASH_SIZE; smallh++) {
        repeattmp = repeatheadptr[smallh];
        while (repeattmp != NULL) {
            /* Find matching entry in main headptr */
            h = hash_function(C, repeattmp->value);
            tmp = headptr[h];
            while (tmp != NULL) {
                if (lmermatch(C, C->sequence + tmp->lastocc, repeattmp->value))
                    break;
                else if (lmermatchrc(C, C->sequence + tmp->lastocc, repeattmp->value))
                    break;
                tmp = tmp->next;
            }
            if (tmp == NULL) { repeattmp = repeattmp->next; continue; }

            /* Extend each available hit (coverage-aware: position may
             * be re-extended if claim_count < MAX_FAMILY_CLAIMS). */
            postmp = tmp->pos;
            while (postmp != NULL) {
                x = postmp->this;
                if (C->removed[x] >= MAX_FAMILY_CLAIMS) {
                    postmp = postmp->next;
                    continue;
                }

                /* Strand determination */
                rrev = 1;
                if (lmermatch(C, C->sequence + x, repeattmp->value))
                    rrev = 0;

                /* 1-vs-1 banded DP extension */
                bestsequencey_val = maskextend_right(C, rrev, x, repeattmp->repeatpos);
                bestsequencew_val = maskextend_left(C, rrev, x, repeattmp->repeatpos);

                /* Determine mask range */
                if (EXTRAMASK) {
                    startmask = bestsequencew_val - C->l + 1;
                    if (startmask < 0) startmask = 0;
                    endmask = bestsequencey_val;
                } else {
                    startmask = bestsequencew_val;
                    endmask = bestsequencey_val - C->l + 1;
                }

                /* Mask and update frequencies — restored original
                 * decrement behavior.  Empirical test (chr4 RM-remap):
                 * preserving bystander l-mer freq HURT recall (0.389 →
                 * 0.362), not helped.  Reason: bystander l-mers whose
                 * occurrences are MOSTLY in mask shadow keep high freq
                 * → find_besttmp re-picks them → build_pos finds all
                 * positions filtered → wasted MAXR slot.  The original
                 * freq decrement naturally drained these bystander
                 * l-mers, which is the desired behavior. */
                for (y = startmask; y <= endmask; y++) {
                    uint8_t prev_claim = (uint8_t)C->removed[y];
                    if (prev_claim >= MAX_FAMILY_CLAIMS) continue;

                    if (prev_claim == 0) {
                        h2 = hash_function(C, C->sequence + y);
                        if (h2 < 0) {
                            C->removed[y] = (char)(prev_claim + 1);
                            continue;
                        }
                        tmp2 = headptr[h2];
                        while (tmp2 != NULL) {
                            if (lmermatcheither(C, C->sequence + tmp2->lastocc, C->sequence + y)) {
                                tmp2->freq -= 1;
                                break;
                            }
                            tmp2 = tmp2->next;
                        }
                    }
                    C->removed[y] = (char)(prev_claim + 1);
                }

                postmp = postmp->next;
            }
            repeattmp = repeattmp->next;
        }
    }

    /* Free repeatheadptr */
    for (smallh = 0; smallh < SMALLHASH_SIZE; smallh++) {
        repeattmp = repeatheadptr[smallh];
        while (repeattmp != NULL) {
            nextrepeattmp = repeattmp->next;
            free(repeattmp->value);
            free(repeattmp);
            repeattmp = nextrepeattmp;
        }
    }
    free(repeatheadptr);
}
