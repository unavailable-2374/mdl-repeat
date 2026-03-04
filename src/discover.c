/*
 * discover.c — RepeatScout-style seed-and-extend discovery engine
 *
 * Faithfully replicates RepeatScout v1.0.7 algorithms:
 *   - N-sequence simultaneous banded DP extension
 *   - Mask-and-extend loop for family discovery
 *   - 1-vs-1 banded DP masking
 *
 * Bug fixes vs RepeatScout:
 *   - Sequence boundary check: uses correctly computed bStart/bEnd
 *     (RepeatScout hardcodes bStart=0; bEnd=50000000)
 *   - repeatllist.value: dynamically allocated (fixes overflow when l > 19)
 *
 * The discovery loop is single-threaded (masking global state requires it).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "discover.h"

/* ================================================================
 * Constants
 * ================================================================ */

#define HASH_SIZE       16000057   /* prime, same as RepeatScout */
#define SMALLHASH_SIZE  5003       /* prime, for mask consensus lookup */
#define SMALLL          6          /* length of small hash l-mers */
#define EXTRAMASK       1          /* extend mask region by l-1 on each side */

/* ================================================================
 * Internal hash table types (RepeatScout convention)
 * ================================================================ */

struct posllist {
    int this;
    struct posllist *next;
};

struct llist {
    int freq;
    int lastocc;
    int lastplusocc;       /* last forward-strand occurrence position */
    int lastminusocc;      /* last reverse-complement occurrence position */
    struct llist *next;
    struct posllist *pos;
};

/* Fixed: dynamic value[] instead of char value[20] */
struct repeatllist {
    char *value;           /* dynamically allocated, length = l */
    int   repeatpos;
    struct repeatllist *next;
};

/* ================================================================
 * Static workspace (file-scoped, single-threaded)
 * ================================================================ */

/* Parameters */
static int l;                  /* l-mer length */
static int L;                  /* max extension per side */
static int MAXOFFSET;
static int MAXN;
static int MAXR;
static int MATCH;
static int MISMATCH;
static int GAP;
static int CAPPENALTY;
static int MINIMPROVEMENT;
static int WHEN_TO_STOP;
static float MAXENTROPY;
static int GOODLENGTH;
static int MINTHRESH;
static int TANDEMDIST;
static int VERBOSE;

/* Genome data */
static char *sequence;          /* doubled genome: forward + PAD + RC */
static char *sequence_owned;    /* if non-NULL, we allocated this (for doubling) */
static char *removed;           /* masking bitmap, length = doubled genome length */
static int   length;            /* total length (doubled: 2*orig + PADLENGTH) */
static int   orig_length;       /* original genome length (before doubling) */
static int  *disc_boundaries;   /* pre-padding boundary positions (int copy) */
static int   disc_num_sequences;

/* Extension workspace */
static char  *master;           /* [2*L+l] working consensus */
static char **masters;          /* [MAXR] saved consensus arrays */
static int   *masters_allocated;/* [MAXR] allocation flags */
static int   *masterstart;      /* [MAXR] consensus start in master coords */
static int   *masterend;        /* [MAXR] consensus end in master coords */

static int   *pos;              /* [MAXN] occurrence positions */
static char  *rev;              /* [MAXN] strand flags (0=fwd, 1=rev) */
static int   *upperBoundI;      /* [MAXN] boundary index per occurrence */
static int    N;                /* current number of occurrences */

static int ***score;            /* [2][MAXN][2*MAXOFFSET+1] DP scores */
static int  **score_of_besty;   /* [MAXN][2*MAXOFFSET+1] checkpoint scores */
static int  **maskscore;        /* [2][2*MAXOFFSET+1] mask DP scores */

static int    totalbestscore;
static int    besttotalbestscore;
static int   *bestbestscore;    /* [MAXN] per-occurrence best scores */
static int   *savebestscore;    /* [MAXN] checkpoint of bestbestscore */

static int   *best_left_offset; /* [MAXN] */
static int   *save_left_offset; /* [MAXN] */
static int   *best_right_offset;/* [MAXN] */
static int   *save_right_offset;/* [MAXN] */

static int    besty, bestw;
static int    nrepeatocc, nactiverepeatocc;
static int    bestnrepeatocc, bestnactiverepeatocc;
static int    R;
static int    prevbestfreq, prevbesthash;

/* ================================================================
 * DNA utility functions
 * ================================================================ */

static inline char compl_base(char c)
{
    if (c == DNA_N) return DNA_N;
    return 3 - c;
}

/* ================================================================
 * Hash functions (symmetric w.r.t. reverse complement)
 * ================================================================ */

static int hash_function(const char *lmer)
{
    int x, ans, ans2;

    for (x = 0; x < l; x++)
        if (lmer[x] == DNA_N) return -1;

    ans = 0;
    for (x = 0; x < l; x++)
        ans = (4 * ans + (lmer[x] % 4)) % HASH_SIZE;
    ans2 = 0;
    for (x = 0; x < l; x++)
        ans2 = (4 * ans2 + ((3 - lmer[l - 1 - x]) % 4)) % HASH_SIZE;
    if (ans2 > ans) ans = ans2;
    return ans;
}

static int smallhash_function(const char *lmer)
{
    int x, ans;

    for (x = 0; x < SMALLL; x++)
        if (lmer[x] == DNA_N) return -1;

    ans = lmer[0];
    for (x = 1; x < SMALLL; x++)
        ans = (4 * ans + (lmer[x] % 4)) % SMALLHASH_SIZE;
    return ans;
}

/* ================================================================
 * L-mer matching functions
 * ================================================================ */

/* Exact forward match */
static int lmermatch(const char *lmer1, const char *lmer2)
{
    for (int x = 0; x < l; x++)
        if (lmer1[x] != lmer2[x]) return 0;
    return 1;
}

/* Reverse complement match */
static int lmermatchrc(const char *lmer1, const char *lmer2)
{
    for (int x = 0; x < l; x++)
        if (lmer1[x] + lmer2[l - 1 - x] != 3) return 0;
    return 1;
}

/* Forward or reverse complement match */
static int lmermatcheither(const char *lmer1, const char *lmer2)
{
    for (int x = 0; x < l; x++)
        if (lmer1[x] != lmer2[x]) return lmermatchrc(lmer1, lmer2);
    return 1;
}

/* ================================================================
 * Entropy (Shannon, natural log — returns negative value)
 * ================================================================ */

static double compute_entropy(const char *lmer)
{
    int count[4] = {0, 0, 0, 0};
    double answer = 0.0;

    for (int x = 0; x < l; x++)
        count[(int)lmer[x]] += 1;

    for (int x = 0; x < 4; x++) {
        if (count[x] == 0) continue;
        double y = (double)count[x] / (double)l;
        answer += y * log(y);
    }
    return answer;
}

/* ================================================================
 * Sequence boundary helper (FIXED — no hardcoded override)
 * ================================================================ */

static void get_boundaries(int n, int *bStart_out, int *bEnd_out)
{
    int bStart = PADLENGTH;
    int bEnd;

    if (upperBoundI[n] == -1) {
        /* Defensive: position past all boundaries */
        bStart = PADLENGTH;
        bEnd = length;
    } else {
        if (upperBoundI[n] > 0)
            bStart = disc_boundaries[upperBoundI[n] - 1] + PADLENGTH;
        bEnd = disc_boundaries[upperBoundI[n]] + PADLENGTH;
    }

    /* RepeatScout bug: bStart=0; bEnd=50000000; — NOT applied here */

    *bStart_out = bStart;
    *bEnd_out = bEnd;
}

/* ================================================================
 * Find sequence index for a padded position
 * ================================================================ */

static int find_seq_index(int padded_pos)
{
    int raw = padded_pos - PADLENGTH;
    if (raw < 0) return -1;

    for (int i = 0; i < disc_num_sequences; i++) {
        if (i == 0) {
            if (raw < disc_boundaries[0]) return 0;
        } else {
            if (raw >= disc_boundaries[i - 1] && raw < disc_boundaries[i])
                return i;
        }
    }
    return disc_num_sequences - 1;
}

/* ================================================================
 * Hash table: read from .freq file (RepeatScout format)
 * ================================================================ */

static void build_headptr_from_freq(struct llist **headptr, const char *freq_file)
{
    FILE *fp;
    char string[1000];
    int thisfreq, thisocc, x, h;
    struct llist *tmp;

    fp = fopen(freq_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "discover: could not open freq file %s\n", freq_file);
        exit(1);
    }

    for (h = 0; h < HASH_SIZE; h++)
        headptr[h] = NULL;

    while (1) {
        if (fscanf(fp, "%s", string) != 1) break;
        if (fscanf(fp, "%d", &thisfreq) != 1) {
            fprintf(stderr, "discover: error reading frequency from %s\n", freq_file);
            exit(1);
        }
        if (fscanf(fp, "%d", &thisocc) != 1) {
            fprintf(stderr, "discover: error reading occurrence from %s\n", freq_file);
            exit(1);
        }

        if ((int)strlen(string) != l) {
            fprintf(stderr, "discover: l-mer length mismatch: expected %d, got %d\n",
                    l, (int)strlen(string));
            exit(1);
        }

        for (x = 0; x < l; x++)
            string[x] = char_to_num(string[x]);

        h = hash_function(string);
        if (h < 0) continue;

        /* Check if already present (duplicate = end of useful data) */
        tmp = headptr[h];
        while (tmp != NULL) {
            if (lmermatcheither(sequence + tmp->lastocc, string) == 1)
                break;
            tmp = tmp->next;
        }
        if (tmp != NULL) break; /* already in table — done */

        /* Entropy filter */
        if (compute_entropy(sequence + thisocc) > MAXENTROPY) continue;

        tmp = malloc(sizeof(*tmp));
        if (tmp == NULL) { fprintf(stderr, "discover: out of memory\n"); exit(1); }
        tmp->freq = thisfreq;
        tmp->lastocc = thisocc;
        tmp->lastplusocc = thisocc;
        tmp->lastminusocc = -1000000;
        tmp->next = headptr[h];
        tmp->pos = NULL;
        headptr[h] = tmp;
    }

    fclose(fp);
}

/* ================================================================
 * Hash table: internal counting (no .freq file)
 * ================================================================ */

static void build_headptr_internal(struct llist **headptr)
{
    int h, x;
    struct llist *tmp;

    for (h = 0; h < HASH_SIZE; h++)
        headptr[h] = NULL;

    /*
     * Count l-mer frequencies with per-strand TANDEMDIST filtering,
     * matching RepeatScout's build_lmer_table behavior.
     * Only scan the forward copy (not the RC copy in doubled genome).
     */
    int count_end = orig_length - l;
    for (x = PADLENGTH; x <= count_end; x++) {
        h = hash_function(sequence + x);
        if (h < 0) continue;

        /* Check if already in table */
        tmp = headptr[h];
        while (tmp != NULL) {
            if (lmermatch(sequence + tmp->lastplusocc, sequence + x)) {
                /* Forward match: TANDEMDIST check against lastplusocc */
                if (x - tmp->lastplusocc >= TANDEMDIST)
                    tmp->freq++;
                tmp->lastplusocc = x;
                tmp->lastocc = x;
                break;
            } else if (lmermatchrc(sequence + tmp->lastplusocc, sequence + x)) {
                /* Reverse complement match: TANDEMDIST check against lastminusocc */
                if (x - tmp->lastminusocc >= TANDEMDIST)
                    tmp->freq++;
                tmp->lastminusocc = x;
                tmp->lastocc = x;
                break;
            }
            tmp = tmp->next;
        }

        if (tmp == NULL) {
            /* New entry — apply entropy filter */
            if (compute_entropy(sequence + x) > MAXENTROPY) continue;

            tmp = malloc(sizeof(*tmp));
            if (tmp == NULL) { fprintf(stderr, "discover: out of memory\n"); exit(1); }
            tmp->freq = 1;
            tmp->lastocc = x;
            tmp->lastplusocc = x;
            tmp->lastminusocc = -1000000;
            tmp->next = headptr[h];
            tmp->pos = NULL;
            headptr[h] = tmp;
        }
    }
}

/* ================================================================
 * Trim: remove l-mers with freq < MINTHRESH, reset freq to 0
 * ================================================================ */

static void trim_headptr(struct llist **headptr)
{
    int h;
    struct llist *tmp, *prevtmp, *nexttmp;

    for (h = 0; h < HASH_SIZE; h++) {
        prevtmp = NULL;
        tmp = headptr[h];
        while (tmp != NULL) {
            if (tmp->freq >= MINTHRESH) {
                /* Keep, but reset freq for rebuild by build_all_pos */
                tmp->freq = 0;
                prevtmp = tmp;
                tmp = tmp->next;
                continue;
            }
            /* Remove */
            nexttmp = tmp->next;
            free(tmp);
            tmp = nexttmp;
            if (prevtmp == NULL)
                headptr[h] = tmp;
            else
                prevtmp->next = tmp;
        }
    }
}

/* ================================================================
 * Build all positions: rescan genome, add positions, TANDEMDIST filter
 * (Exact port of RepeatScout build_all_pos, lines 613-743)
 * ================================================================ */

static void build_all_pos(struct llist **headptr)
{
    int x, h, pos1, pos2;
    struct llist *tmp;
    struct posllist *postmp, *postmp2, *prevpostmp, *nextpostmp;
    int currBoundary = 0;

    if (VERBOSE) fprintf(stderr, "discover: building all positions...\n");

    /* Pass 1: scan genome, add positions to matching l-mers */
    for (x = l - 1; x < length - l + 1; x++) {
        /* Check sequence boundaries */
        if (disc_boundaries[currBoundary] > 0) {
            if (x == disc_boundaries[currBoundary] + PADLENGTH)
                currBoundary++;
            if (x + l > disc_boundaries[currBoundary] + PADLENGTH)
                continue;
        }

        h = hash_function(sequence + x);
        if (h < 0) continue;

        tmp = headptr[h];
        while (tmp != NULL) {
            if (lmermatcheither(sequence + tmp->lastocc, sequence + x)) {
                /* Hit: add position */
                postmp = malloc(sizeof(*postmp));
                if (postmp == NULL) {
                    fprintf(stderr, "discover: out of memory\n"); exit(1);
                }
                postmp->this = x;
                postmp->next = tmp->pos;
                tmp->pos = postmp;
                tmp->freq++;
            }
            tmp = tmp->next;
        }
    }

    if (VERBOSE) fprintf(stderr, "discover: TANDEMDIST filtering...\n");

    /* Pass 2: remove position pairs within TANDEMDIST (same strand only) */
    for (h = 0; h < HASH_SIZE; h++) {
        tmp = headptr[h];
        while (tmp != NULL) {
            /* Mark within-TANDEMDIST same-strand pairs as negative */
            postmp = tmp->pos;
            while (postmp != NULL && postmp->next != NULL) {
                pos1 = postmp->this;
                postmp2 = postmp;
                while ((postmp2 = postmp2->next) != NULL) {
                    pos2 = postmp2->this;
                    if (pos1 - pos2 >= TANDEMDIST)
                        break;
                    if (lmermatch(sequence + pos1, sequence + pos2)) {
                        postmp->this = -pos1;
                        break;
                    }
                }
                postmp = postmp->next;
            }

            /* Remove marked (negative) positions */
            postmp = tmp->pos;
            prevpostmp = NULL;
            while (postmp != NULL) {
                if (postmp->this >= 0) {
                    prevpostmp = postmp;
                    postmp = postmp->next;
                    continue;
                }
                /* Remove */
                nextpostmp = postmp->next;
                tmp->freq--;
                removed[-postmp->this] = 1;
                free(postmp);
                postmp = nextpostmp;
                if (prevpostmp == NULL)
                    tmp->pos = postmp;
                else
                    prevpostmp->next = postmp;
            }

            tmp = tmp->next;
        }
    }

    if (VERBOSE) fprintf(stderr, "discover: position building complete\n");
}

/* ================================================================
 * Seed selection: find highest-frequency l-mer
 * (Exact port of RepeatScout find_besttmp, lines 748-788)
 * ================================================================ */

static struct llist *find_besttmp(struct llist **headptr)
{
    int h;
    struct llist *tmp, *besttmp;
    int bestfreq;

    /* Try to match prevbestfreq first (locality optimization) */
    for (h = prevbesthash; h < HASH_SIZE; h++) {
        tmp = headptr[h];
        while (tmp != NULL) {
            if (tmp->freq == prevbestfreq) {
                prevbesthash = h;
                return tmp;
            }
            tmp = tmp->next;
        }
    }

    /* Global search for max freq */
    besttmp = NULL;
    bestfreq = 0;
    for (h = 0; h < HASH_SIZE; h++) {
        tmp = headptr[h];
        while (tmp != NULL) {
            if (tmp->freq > bestfreq) {
                besttmp = tmp;
                bestfreq = tmp->freq;
                prevbesthash = h;
            }
            tmp = tmp->next;
        }
    }
    prevbestfreq = bestfreq;
    return besttmp;
}

/* ================================================================
 * Build occurrence list for a seed
 * (Exact port of RepeatScout build_pos, lines 996-1032)
 * ================================================================ */

static void build_pos(struct llist *besttmp)
{
    struct posllist *postmp;
    int x;

    postmp = besttmp->pos;
    N = 0;
    while (postmp != NULL) {
        if (removed[postmp->this] == 0) {
            pos[N] = postmp->this;

            /* Strand determination: lmermatch = forward, else reverse */
            if (lmermatch(sequence + besttmp->lastocc, sequence + postmp->this))
                rev[N] = 0;
            else
                rev[N] = 1;

            /* Find upper boundary index */
            upperBoundI[N] = -1;
            x = 0;
            while (disc_boundaries[x] != 0) {
                if (disc_boundaries[x] + PADLENGTH > pos[N]) {
                    upperBoundI[N] = x;
                    break;
                }
                x++;
            }

            N++;
            if (N == MAXN) break;
        }
        postmp = postmp->next;
    }
}

/* ================================================================
 * DP scoring: compute_score_right
 * (Exact port of RepeatScout lines 1877-1969, with boundary fix)
 * ================================================================ */

static int compute_score_right(int y, int n, int offset, char a)
{
    int oldoffset, tempscore, ans, ismatch, x;
    int bStart, bEnd;

    get_boundaries(n, &bStart, &bEnd);

    /* Boundary check */
    if (rev[n]) {
        if (pos[n] - (offset + y - L - l) - 1 < bStart)
            return 0;
    } else {
        if (pos[n] + offset + y - L >= bEnd)
            return 0;
    }

    ans = -1000000000;

    /* Case A: gap in sequence (oldoffset = offset+1) */
    if (offset < MAXOFFSET) {
        oldoffset = offset + 1;
        tempscore = score[(y - 1) % 2][n][oldoffset + MAXOFFSET] + GAP;
        if (tempscore > ans) ans = tempscore;
    }

    /* Case B: diagonal match/mismatch (oldoffset = offset) */
    oldoffset = offset;
    tempscore = score[(y - 1) % 2][n][oldoffset + MAXOFFSET];
    if (rev[n]) {
        if (a == compl_base(sequence[pos[n] - (offset + y - L - l) - 1]))
            tempscore += MATCH;
        else
            tempscore += MISMATCH;
    } else {
        if (a == sequence[pos[n] + offset + y - L])
            tempscore += MATCH;
        else
            tempscore += MISMATCH;
    }
    if (tempscore > ans) ans = tempscore;

    /* Case C: multiple gaps in consensus (oldoffset < offset) */
    for (oldoffset = -MAXOFFSET; oldoffset < offset; oldoffset++) {
        ismatch = 0;
        for (x = oldoffset; x <= offset; x++) {
            if (rev[n]) {
                if (a == compl_base(sequence[pos[n] - (x + y - L - l) - 1]))
                    ismatch = 1;
            } else {
                if (a == sequence[pos[n] + x + y - L])
                    ismatch = 1;
            }
        }
        tempscore = score[(y - 1) % 2][n][oldoffset + MAXOFFSET];
        tempscore += (offset - oldoffset) * GAP;
        if (ismatch) tempscore += MATCH;
        else tempscore += MISMATCH;
        if (tempscore > ans) ans = tempscore;
    }

    return ans;
}

/* ================================================================
 * DP scoring: compute_score_left
 * (Exact port of RepeatScout lines 1976-2069, with boundary fix)
 * ================================================================ */

static int compute_score_left(int w, int n, int offset, char a)
{
    int oldoffset, tempscore, ans, ismatch, x;
    int bStart, bEnd;

    get_boundaries(n, &bStart, &bEnd);

    /* Boundary check (note: directions are swapped vs right) */
    if (rev[n]) {
        if (pos[n] - (offset + w - L - l) - 1 >= bEnd)
            return 0;
    } else {
        if (pos[n] + offset + w - L < bStart)
            return 0;
    }

    ans = -1000000000;

    /* Case A: gap in sequence (oldoffset = offset-1) */
    if (offset > -MAXOFFSET) {
        oldoffset = offset - 1;
        tempscore = score[(w + 1) % 2][n][oldoffset + MAXOFFSET] + GAP;
        if (tempscore > ans) ans = tempscore;
    }

    /* Case B: diagonal (oldoffset = offset) */
    oldoffset = offset;
    tempscore = score[(w + 1) % 2][n][oldoffset + MAXOFFSET];
    if (rev[n]) {
        if (a == compl_base(sequence[pos[n] - (offset + w - L - l) - 1]))
            tempscore += MATCH;
        else
            tempscore += MISMATCH;
    } else {
        if (a == sequence[pos[n] + offset + w - L])
            tempscore += MATCH;
        else
            tempscore += MISMATCH;
    }
    if (tempscore > ans) ans = tempscore;

    /* Case C: multiple gaps in consensus (oldoffset > offset) */
    for (oldoffset = offset + 1; oldoffset <= MAXOFFSET; oldoffset++) {
        ismatch = 0;
        for (x = offset; x <= oldoffset; x++) {
            if (rev[n]) {
                if (a == compl_base(sequence[pos[n] - (x + w - L - l) - 1]))
                    ismatch = 1;
            } else {
                if (a == sequence[pos[n] + x + w - L])
                    ismatch = 1;
            }
        }
        tempscore = score[(w + 1) % 2][n][oldoffset + MAXOFFSET];
        tempscore += (oldoffset - offset) * GAP;
        if (ismatch) tempscore += MATCH;
        else tempscore += MISMATCH;
        if (tempscore > ans) ans = tempscore;
    }

    return ans;
}

/* ================================================================
 * Checkpoint: compute_totalbestscore_right
 * (Exact port of RepeatScout lines 1722-1776)
 * ================================================================ */

static void compute_totalbestscore_right(int y)
{
    int n, bestscore, offset;

    nrepeatocc = 0;
    nactiverepeatocc = 0;
    totalbestscore = 0;

    for (n = 0; n < N; n++) {
        bestscore = bestbestscore[n] + CAPPENALTY;
        if (bestscore < 0) bestscore = 0;

        int bestscore_offset = -100000;
        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
            if (score[y % 2][n][offset + MAXOFFSET] > bestscore) {
                bestscore = score[y % 2][n][offset + MAXOFFSET];
                bestscore_offset = offset;
            }
        }

        if (bestscore > 0) nrepeatocc++;
        if (bestscore > bestbestscore[n] + CAPPENALTY) nactiverepeatocc++;
        if (bestscore > bestbestscore[n]) {
            if (bestscore_offset >= -MAXOFFSET)
                best_right_offset[n] = bestscore_offset + y;
            bestbestscore[n] = bestscore;
        }
        totalbestscore += bestscore;
    }

    /* MINIMPROVEMENT checkpoint */
    if ((totalbestscore >= besttotalbestscore + (y - besty) * MINIMPROVEMENT)
        && (totalbestscore > besttotalbestscore)) {
        besty = y;
        besttotalbestscore = totalbestscore;
        bestnrepeatocc = nrepeatocc;
        bestnactiverepeatocc = nactiverepeatocc;
        for (n = 0; n < N; n++) {
            save_right_offset[n] = best_right_offset[n];
            savebestscore[n] = bestbestscore[n];
            for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
                score_of_besty[n][offset + MAXOFFSET] = score[y % 2][n][offset + MAXOFFSET];
        }
    }
}

/* ================================================================
 * Checkpoint: compute_totalbestscore_left
 * (Exact port of RepeatScout lines 1795-1841)
 *
 * KEY ASYMMETRY: does NOT update savebestscore or score_of_besty.
 * Only updates bestw, besttotalbestscore, save_left_offset.
 * This is intentional — see plan review #1.
 * ================================================================ */

static void compute_totalbestscore_left(int w)
{
    int n, bestscore, offset;

    nrepeatocc = 0;
    nactiverepeatocc = 0;
    totalbestscore = 0;

    for (n = 0; n < N; n++) {
        int bestscore_offset = -100000;
        bestscore = bestbestscore[n] + CAPPENALTY;
        if (bestscore < 0) bestscore = 0;

        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
            if (score[w % 2][n][offset + MAXOFFSET] > bestscore) {
                bestscore = score[w % 2][n][offset + MAXOFFSET];
                bestscore_offset = offset;
            }
        }

        if (bestscore > 0) nrepeatocc++;
        if (bestscore > bestbestscore[n] + CAPPENALTY) nactiverepeatocc++;
        if (bestscore > bestbestscore[n]) {
            if (bestscore_offset >= -MAXOFFSET)
                best_left_offset[n] = bestscore_offset + w;
            bestbestscore[n] = bestscore;
        }
        totalbestscore += bestscore;
    }

    /* MINIMPROVEMENT checkpoint — only save offsets, NOT scores/score_of_besty */
    if ((totalbestscore >= besttotalbestscore + (bestw - w) * MINIMPROVEMENT)
        && (totalbestscore > besttotalbestscore)) {
        bestw = w;
        besttotalbestscore = totalbestscore;
        bestnrepeatocc = nrepeatocc;
        bestnactiverepeatocc = nactiverepeatocc;
        for (n = 0; n < N; n++)
            save_left_offset[n] = best_left_offset[n];
    }
}

/* ================================================================
 * extend_right
 * (Exact port of RepeatScout lines 1034-1117)
 * ================================================================ */

static void extend_right(void)
{
    int y, n, bestscore, tempscore, offset;
    char a, besta;
    int newtotalbestscore, newtotalbestscore_a;
    besta = 0;

    /* Initialize score at seed's last base */
    y = L + l - 1;
    for (n = 0; n < N; n++) {
        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
            score[y % 2][n][offset + MAXOFFSET] = l * MATCH;
            if (offset < 0)
                score[y % 2][n][offset + MAXOFFSET] += -offset * GAP;
            if (offset > 0)
                score[y % 2][n][offset + MAXOFFSET] += offset * GAP;
        }
        bestbestscore[n] = l * MATCH;
    }

    /* Initialize checkpoint */
    besty = L + l - 1;
    for (n = 0; n < N; n++) {
        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
            score_of_besty[n][offset + MAXOFFSET] = score[y % 2][n][offset + MAXOFFSET];
    }
    besttotalbestscore = 0;
    compute_totalbestscore_right(y);

    /* Extend right, one base at a time */
    for (y = L + l; y < 2 * L + l; y++) {
        newtotalbestscore = 0;
        for (a = 0; a < 4; a++) {
            newtotalbestscore_a = 0;
            for (n = 0; n < N; n++) {
                bestscore = bestbestscore[n] + CAPPENALTY;
                if (bestscore < 0) bestscore = 0;
                for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
                    tempscore = compute_score_right(y, n, offset, a);
                    if (tempscore > bestscore)
                        bestscore = tempscore;
                }
                newtotalbestscore_a += bestscore;
            }
            if (newtotalbestscore_a > newtotalbestscore) {
                newtotalbestscore = newtotalbestscore_a;
                besta = a;
            }
        }

        master[y] = besta;
        for (n = 0; n < N; n++) {
            for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
                score[y % 2][n][offset + MAXOFFSET] = compute_score_right(y, n, offset, besta);
        }
        compute_totalbestscore_right(y);

        if (y - besty >= WHEN_TO_STOP) break;
    }

    if (y == 2 * L + l && VERBOSE)
        fprintf(stderr, "Warning: extended right all the way to %d\n", y);

    y = besty;
    totalbestscore = besttotalbestscore;
    nrepeatocc = bestnrepeatocc;
    nactiverepeatocc = bestnactiverepeatocc;
}

/* ================================================================
 * extend_left
 * (Exact port of RepeatScout lines 1119-1205)
 * ================================================================ */

static void extend_left(void)
{
    int w, n, bestscore, tempscore, offset;
    char a, besta;
    int newtotalbestscore, newtotalbestscore_a;
    besta = 0;

    /* Initialize from right extension checkpoint */
    w = L;
    for (n = 0; n < N; n++) {
        bestscore = savebestscore[n] + CAPPENALTY;
        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
            if (score_of_besty[n][offset + MAXOFFSET] > bestscore)
                bestscore = score_of_besty[n][offset + MAXOFFSET];
        }
        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
            score[w % 2][n][offset + MAXOFFSET] = bestscore;
            if (offset < 0)
                score[w % 2][n][offset + MAXOFFSET] += -offset * GAP;
            if (offset > 0)
                score[w % 2][n][offset + MAXOFFSET] += offset * GAP;
        }
        bestbestscore[n] = bestscore;
    }

    /* Initialize checkpoint */
    bestw = L;
    compute_totalbestscore_left(w);

    /* Extend left */
    for (w = L - 1; w >= 0; w--) {
        newtotalbestscore = 0;
        for (a = 0; a < 4; a++) {
            newtotalbestscore_a = 0;
            for (n = 0; n < N; n++) {
                bestscore = bestbestscore[n] + CAPPENALTY;
                if (bestscore < 0) bestscore = 0;
                for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
                    tempscore = compute_score_left(w, n, offset, a);
                    if (tempscore > bestscore)
                        bestscore = tempscore;
                }
                newtotalbestscore_a += bestscore;
            }
            if (newtotalbestscore_a > newtotalbestscore) {
                newtotalbestscore = newtotalbestscore_a;
                besta = a;
            }
        }

        master[w] = besta;
        for (n = 0; n < N; n++) {
            for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++)
                score[w % 2][n][offset + MAXOFFSET] = compute_score_left(w, n, offset, besta);
        }
        compute_totalbestscore_left(w);

        if (w - bestw <= -WHEN_TO_STOP) break;
    }

    if (w == -1 && VERBOSE)
        fprintf(stderr, "Warning: extended left all the way to %d\n", w);

    w = bestw;
    totalbestscore = besttotalbestscore;
    nrepeatocc = bestnrepeatocc;
    nactiverepeatocc = bestnactiverepeatocc;

    masterstart[R] = bestw;
    masterend[R] = besty;
    for (w = bestw; w <= besty; w++)
        masters[R][w] = master[w];
}

/* ================================================================
 * Mask scoring: compute_maskscore_right
 * (Exact port of RepeatScout lines 2071-2113)
 * ================================================================ */

static int compute_maskscore_right(char mbase, int repeaty, int seqpos, int offset)
{
    int oldoffset, tempscore, ans, ismatch, x;

    ans = -1000000000;

    /* Case A: gap in sequence */
    if (offset < MAXOFFSET) {
        oldoffset = offset + 1;
        tempscore = maskscore[(repeaty + 1) % 2][oldoffset + MAXOFFSET] + GAP;
        if (tempscore > ans) ans = tempscore;
    }

    /* Case B: diagonal */
    oldoffset = offset;
    tempscore = maskscore[(repeaty + 1) % 2][oldoffset + MAXOFFSET];
    if (mbase == sequence[seqpos + l + repeaty + offset])
        tempscore += MATCH;
    else
        tempscore += MISMATCH;
    if (tempscore > ans) ans = tempscore;

    /* Case C: multiple gaps */
    for (oldoffset = -MAXOFFSET; oldoffset < offset; oldoffset++) {
        ismatch = 0;
        for (x = oldoffset; x <= offset; x++) {
            if (mbase == sequence[seqpos + l + repeaty + x])
                ismatch = 1;
        }
        tempscore = maskscore[(repeaty + 1) % 2][oldoffset + MAXOFFSET];
        tempscore += (offset - oldoffset) * GAP;
        if (ismatch) tempscore += MATCH;
        else tempscore += MISMATCH;
        if (tempscore > ans) ans = tempscore;
    }

    return ans;
}

/* ================================================================
 * Mask scoring: compute_maskscore_left
 * (Exact port of RepeatScout lines 2116-2161)
 * ================================================================ */

static int compute_maskscore_left(char mbase, int repeatw, int seqpos, int offset)
{
    int oldoffset, tempscore, ans, ismatch, x;

    ans = -1000000000;

    /* Case A: gap in sequence */
    if (offset > -MAXOFFSET) {
        oldoffset = offset - 1;
        tempscore = maskscore[(repeatw + 1) % 2][oldoffset + MAXOFFSET] + GAP;
        if (tempscore > ans) ans = tempscore;
    }

    /* Case B: diagonal */
    oldoffset = offset;
    tempscore = maskscore[(repeatw + 1) % 2][oldoffset + MAXOFFSET];
    if (mbase == sequence[seqpos - repeatw - 1 + offset])
        tempscore += MATCH;
    else
        tempscore += MISMATCH;
    if (tempscore > ans) ans = tempscore;

    /* Case C: multiple gaps */
    for (oldoffset = offset + 1; oldoffset <= MAXOFFSET; oldoffset++) {
        ismatch = 0;
        for (x = offset; x <= oldoffset; x++) {
            if (mbase == sequence[seqpos - repeatw - 1 + x])
                ismatch = 1;
        }
        tempscore = maskscore[(repeatw + 1) % 2][oldoffset + MAXOFFSET];
        tempscore += (oldoffset - offset) * GAP;
        if (ismatch) tempscore += MATCH;
        else tempscore += MISMATCH;
        if (tempscore > ans) ans = tempscore;
    }

    return ans;
}

/* ================================================================
 * maskextend_right
 * (Exact port of RepeatScout lines 1387-1472)
 * ================================================================ */

static int maskextend_right(int rc, int seqpos, int modelpos)
{
    int offset, bestmaskscore_val, bestbestmaskscore, bestoffset, bestbestoffset;
    int bEnd, x, bestx, maxext;
    char mbase;

    /* Initialize maskscore */
    for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
        maskscore[1][offset + MAXOFFSET] = 0;
        if (offset < 0) maskscore[1][offset + MAXOFFSET] += -offset * GAP;
        if (offset > 0) maskscore[1][offset + MAXOFFSET] += offset * GAP;
    }

    bestbestmaskscore = 0;
    bestbestoffset = 0;
    bestoffset = 0;
    bestx = 0;

    /* Find sequence boundary */
    bEnd = 0;
    x = 0;
    while (disc_boundaries[x] > 0) {
        if (seqpos < disc_boundaries[x] + PADLENGTH) {
            bEnd = disc_boundaries[x] + PADLENGTH;
            break;
        }
        x++;
    }

    /* Max extension distance in model */
    maxext = masterend[R - 1] - (modelpos + l + 1);
    if (rc)
        maxext = modelpos - (masterstart[R - 1] + 1);

    for (x = 0; x < maxext; x++) {
        if (bEnd > 0 && seqpos + l + x == bEnd)
            break;

        /* Get consensus base (4 combinations of direction × strand) */
        if (rc)
            mbase = compl_base(masters[R - 1][modelpos - x - 1]);
        else
            mbase = masters[R - 1][modelpos + l + x];

        bestmaskscore_val = -1000000000;
        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
            if (bEnd > 0 && seqpos + l + x + offset >= bEnd) {
                maskscore[(x + 2) % 2][offset + MAXOFFSET] = 0;
            } else {
                maskscore[(x + 2) % 2][offset + MAXOFFSET] =
                    compute_maskscore_right(mbase, x, seqpos, offset);
            }
            if (maskscore[(x + 2) % 2][offset + MAXOFFSET] > bestmaskscore_val) {
                bestmaskscore_val = maskscore[(x + 2) % 2][offset + MAXOFFSET];
                bestoffset = offset;
            }
        }
        if (bestmaskscore_val > bestbestmaskscore) {
            bestbestmaskscore = bestmaskscore_val;
            bestbestoffset = bestoffset;
            bestx = x;
        }
        if (x - bestx >= WHEN_TO_STOP) break;
    }

    return seqpos + l + bestx + bestbestoffset;
}

/* ================================================================
 * maskextend_left
 * (Exact port of RepeatScout lines 1496-1578)
 * ================================================================ */

static int maskextend_left(int rc, int seqpos, int modelpos)
{
    int offset, bestmaskscore_val, bestbestmaskscore, bestoffset, bestbestoffset;
    int maxext, bStart, x, bestx;
    char mbase;

    /* Initialize maskscore */
    for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
        maskscore[1][offset + MAXOFFSET] = 0;
        if (offset < 0) maskscore[1][offset + MAXOFFSET] += -offset * GAP;
        if (offset > 0) maskscore[1][offset + MAXOFFSET] += offset * GAP;
    }

    bestbestoffset = 0;
    bestoffset = 0;
    bestbestmaskscore = 0;
    bestx = 0;

    /* Find sequence boundary */
    bStart = 0;
    x = 0;
    while (disc_boundaries[x] > 0) {
        if (seqpos < disc_boundaries[x] + PADLENGTH) {
            if (x > 0)
                bStart = disc_boundaries[x - 1] + PADLENGTH;
            break;
        }
        x++;
    }

    /* Max extension distance in model */
    maxext = modelpos - masterstart[R - 1];
    if (rc)
        maxext = masterend[R - 1] - (modelpos + l - 1);

    for (x = 0; x < maxext; x++) {
        if (seqpos - x - 1 == bStart)
            break;

        /* Get consensus base */
        if (rc)
            mbase = compl_base(masters[R - 1][modelpos + l + x]);
        else
            mbase = masters[R - 1][modelpos - x - 1];

        bestmaskscore_val = -1000000000;
        for (offset = -MAXOFFSET; offset <= MAXOFFSET; offset++) {
            if (seqpos - x - 1 + offset < bStart) {
                maskscore[(x + 2) % 2][offset + MAXOFFSET] = 0;
            } else {
                maskscore[(x + 2) % 2][offset + MAXOFFSET] =
                    compute_maskscore_left(mbase, x, seqpos, offset);
            }
            if (maskscore[(x + 2) % 2][offset + MAXOFFSET] > bestmaskscore_val) {
                bestmaskscore_val = maskscore[(x + 2) % 2][offset + MAXOFFSET];
                bestoffset = offset;
            }
        }
        if (bestmaskscore_val > bestbestmaskscore) {
            bestbestmaskscore = bestmaskscore_val;
            bestbestoffset = bestoffset;
            bestx = x;
        }
        if (x - bestx >= WHEN_TO_STOP) break;
    }

    return seqpos - bestx + bestbestoffset - 1;
}

/* ================================================================
 * mask_headptr: mask found family from genome
 * (Exact port of RepeatScout lines 1207-1364, with dynamic value[])
 * ================================================================ */

static void mask_headptr(struct llist **headptr)
{
    int x, y, smallh, h, h2;
    int startmask, endmask;
    int bestsequencey_val, bestsequencew_val;
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

    for (x = masterstart[R - 1]; x <= masterend[R - 1] - l + 1; x++) {
        smallh = smallhash_function(masters[R - 1] + x);
        if (smallh < 0) continue;

        /* Check if already present */
        repeattmp = repeatheadptr[smallh];
        while (repeattmp != NULL) {
            if (lmermatch(repeattmp->value, masters[R - 1] + x))
                break;
            repeattmp = repeattmp->next;
        }
        if (repeattmp != NULL) continue;

        /* Add new entry (FIXED: dynamic value allocation) */
        repeattmp = malloc(sizeof(*repeattmp));
        if (repeattmp == NULL) {
            fprintf(stderr, "discover: out of memory\n"); exit(1);
        }
        repeattmp->value = malloc(l * sizeof(char));
        if (repeattmp->value == NULL) {
            fprintf(stderr, "discover: out of memory for value\n"); exit(1);
        }
        for (y = 0; y < l; y++)
            repeattmp->value[y] = masters[R - 1][x + y];
        repeattmp->repeatpos = x;
        repeattmp->next = repeatheadptr[smallh];
        repeatheadptr[smallh] = repeattmp;
    }

    /* Find and extend hits */
    for (smallh = 0; smallh < SMALLHASH_SIZE; smallh++) {
        repeattmp = repeatheadptr[smallh];
        while (repeattmp != NULL) {
            /* Find matching entry in main headptr */
            h = hash_function(repeattmp->value);
            tmp = headptr[h];
            while (tmp != NULL) {
                if (lmermatch(sequence + tmp->lastocc, repeattmp->value))
                    break;
                else if (lmermatchrc(sequence + tmp->lastocc, repeattmp->value))
                    break;
                tmp = tmp->next;
            }
            if (tmp == NULL) { repeattmp = repeattmp->next; continue; }

            /* Extend each unmasked hit */
            postmp = tmp->pos;
            while (postmp != NULL) {
                x = postmp->this;
                if (removed[x] == 1) { postmp = postmp->next; continue; }

                /* Strand determination */
                rrev = 1;
                if (lmermatch(sequence + x, repeattmp->value))
                    rrev = 0;

                /* 1-vs-1 banded DP extension */
                bestsequencey_val = maskextend_right(rrev, x, repeattmp->repeatpos);
                bestsequencew_val = maskextend_left(rrev, x, repeattmp->repeatpos);

                /* Determine mask range */
                if (EXTRAMASK) {
                    startmask = bestsequencew_val - l + 1;
                    if (startmask < 0) startmask = 0;
                    endmask = bestsequencey_val;
                } else {
                    startmask = bestsequencew_val;
                    endmask = bestsequencey_val - l + 1;
                }

                /* Mask and update frequencies */
                for (y = startmask; y <= endmask; y++) {
                    if (removed[y] == 1) continue;
                    h2 = hash_function(sequence + y);
                    if (h2 < 0) continue;
                    tmp2 = headptr[h2];
                    while (tmp2 != NULL) {
                        if (lmermatcheither(sequence + tmp2->lastocc, sequence + y)) {
                            tmp2->freq -= 1;
                            break;
                        }
                        tmp2 = tmp2->next;
                    }
                    removed[y] = 1;
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

/* ================================================================
 * Workspace allocation and deallocation
 * ================================================================ */

static int allocate_workspace(void)
{
    int x, n;

    master = malloc((2 * L + l + 1) * sizeof(char));
    if (!master) return -1;
    master[2 * L + l] = '\0';

    masters_allocated = calloc(MAXR, sizeof(int));
    masterstart = malloc(MAXR * sizeof(int));
    masterend = malloc(MAXR * sizeof(int));
    if (!masters_allocated || !masterstart || !masterend) return -1;

    masters = malloc(MAXR * sizeof(*masters));
    if (!masters) return -1;

    rev = malloc(MAXN * sizeof(char));
    upperBoundI = malloc(MAXN * sizeof(int));
    pos = malloc(MAXN * sizeof(int));
    bestbestscore = malloc(MAXN * sizeof(int));
    savebestscore = malloc(MAXN * sizeof(int));
    best_left_offset = malloc(MAXN * sizeof(int));
    save_left_offset = malloc(MAXN * sizeof(int));
    best_right_offset = malloc(MAXN * sizeof(int));
    save_right_offset = malloc(MAXN * sizeof(int));
    if (!rev || !upperBoundI || !pos || !bestbestscore || !savebestscore ||
        !best_left_offset || !save_left_offset ||
        !best_right_offset || !save_right_offset) return -1;

    /* score[2][MAXN][2*MAXOFFSET+1] */
    score = malloc(2 * sizeof(*score));
    if (!score) return -1;
    for (x = 0; x < 2; x++) {
        score[x] = malloc(MAXN * sizeof(*score[x]));
        if (!score[x]) return -1;
        for (n = 0; n < MAXN; n++) {
            score[x][n] = malloc((2 * MAXOFFSET + 1) * sizeof(*score[x][n]));
            if (!score[x][n]) return -1;
        }
    }

    /* score_of_besty[MAXN][2*MAXOFFSET+1] */
    score_of_besty = malloc(MAXN * sizeof(*score_of_besty));
    if (!score_of_besty) return -1;
    for (n = 0; n < MAXN; n++) {
        score_of_besty[n] = malloc((2 * MAXOFFSET + 1) * sizeof(*score_of_besty[n]));
        if (!score_of_besty[n]) return -1;
    }

    /* maskscore[2][2*MAXOFFSET+1] */
    maskscore = malloc(2 * sizeof(*maskscore));
    if (!maskscore) return -1;
    for (x = 0; x < 2; x++) {
        maskscore[x] = malloc((2 * MAXOFFSET + 1) * sizeof(*maskscore[x]));
        if (!maskscore[x]) return -1;
    }

    return 0;
}

static void free_workspace(void)
{
    int x, n;

    free(master); master = NULL;
    if (masters) {
        for (int r = 0; r < MAXR; r++) {
            if (masters_allocated && masters_allocated[r])
                free(masters[r]);
        }
        free(masters); masters = NULL;
    }
    free(masters_allocated); masters_allocated = NULL;
    free(masterstart); masterstart = NULL;
    free(masterend); masterend = NULL;

    free(rev); rev = NULL;
    free(upperBoundI); upperBoundI = NULL;
    free(pos); pos = NULL;
    free(bestbestscore); bestbestscore = NULL;
    free(savebestscore); savebestscore = NULL;
    free(best_left_offset); best_left_offset = NULL;
    free(save_left_offset); save_left_offset = NULL;
    free(best_right_offset); best_right_offset = NULL;
    free(save_right_offset); save_right_offset = NULL;

    if (score) {
        for (x = 0; x < 2; x++) {
            if (score[x]) {
                for (n = 0; n < MAXN; n++)
                    free(score[x][n]);
                free(score[x]);
            }
        }
        free(score); score = NULL;
    }

    if (score_of_besty) {
        for (n = 0; n < MAXN; n++)
            free(score_of_besty[n]);
        free(score_of_besty); score_of_besty = NULL;
    }

    if (maskscore) {
        for (x = 0; x < 2; x++)
            free(maskscore[x]);
        free(maskscore); maskscore = NULL;
    }

    free(removed); removed = NULL;
    free(disc_boundaries); disc_boundaries = NULL;
    free(sequence_owned); sequence_owned = NULL;
    sequence = NULL;
}

static void free_headptr(struct llist **headptr)
{
    int h;
    struct llist *tmp, *nexttmp;
    struct posllist *postmp, *nextpostmp;

    if (!headptr) return;

    for (h = 0; h < HASH_SIZE; h++) {
        tmp = headptr[h];
        while (tmp != NULL) {
            nexttmp = tmp->next;
            /* Free position chain */
            postmp = tmp->pos;
            while (postmp != NULL) {
                nextpostmp = postmp->next;
                free(postmp);
                postmp = nextpostmp;
            }
            free(tmp);
            tmp = nexttmp;
        }
    }
    free(headptr);
}

/* ================================================================
 * Collect instances from extension occurrences
 *
 * After extend_right + extend_left, pos[0..N-1] holds the seed
 * positions and rev[0..N-1] holds strand flags. We derive
 * approximate instance locations and edit counts.
 * ================================================================ */

static void collect_instances_from_extension(CandidateFamily *fam)
{
    int cons_len = masterend[R] - masterstart[R] + 1;
    int n, i;

    fam->instances = malloc(N * sizeof(Instance));
    if (!fam->instances) {
        fprintf(stderr, "discover: out of memory for instances\n");
        fam->num_instances = 0;
        fam->cap_instances = 0;
        return;
    }
    fam->cap_instances = N;
    fam->num_instances = 0;

    for (n = 0; n < N; n++) {
        gpos_t genome_start, genome_end;

        if (rev[n] == 0) {
            /* Forward strand */
            genome_start = pos[n] + bestw - L;
            genome_end = pos[n] + besty - L + 1;
        } else {
            /* Reverse strand: extension reads genome in opposite direction */
            genome_start = pos[n] - (besty - L - l) - 1;
            genome_end = pos[n] - (bestw - L - l);
        }

        /* Bounds check — only keep instances from forward copy */
        if (genome_start < 0 || genome_end > orig_length)
            continue;

        glen_t aligned_length = genome_end - genome_start;
        if (aligned_length <= 0 || aligned_length > 100000)
            continue;

        /* Count mismatches by direct comparison */
        int edits = 0;
        int compare_len = (int)aligned_length < cons_len ? (int)aligned_length : cons_len;
        for (i = 0; i < compare_len; i++) {
            char seq_base;
            if (rev[n] == 0) {
                seq_base = sequence[genome_start + i];
            } else {
                seq_base = compl_base(sequence[(int)(genome_end - 1 - i)]);
            }
            char cons_base = masters[R][masterstart[R] + i];
            if (seq_base == DNA_N || seq_base != cons_base)
                edits++;
        }

        /* Populate instance */
        Instance *inst = &fam->instances[fam->num_instances];
        inst->position = genome_start;
        inst->aligned_length = aligned_length;
        inst->cons_start = 0;
        inst->cons_end = cons_len;
        inst->num_edits = edits;
        inst->divergence = (compare_len > 0) ? (float)edits / (float)compare_len : 0.0f;
        inst->score = bestbestscore[n];
        inst->strand = rev[n] ? -1 : 1;
        inst->seq_index = find_seq_index((int)genome_start);
        fam->num_instances++;
    }
}

/* ================================================================
 * Genome doubling: append reverse complement (matching RepeatScout add_rc)
 *
 * Creates: [forward_copy | PADLENGTH_padding | RC_copy]
 * Sets length = 2*original_length + PADLENGTH
 * This doubles the occurrences per l-mer, giving stronger extension support.
 * ================================================================ */

static char *build_doubled_genome(const char *orig_seq, int orig_len, int *new_len)
{
    int doubled_len = 2 * orig_len + PADLENGTH;
    char *doubled = malloc(doubled_len * sizeof(char));
    if (!doubled) return NULL;

    /* Copy forward strand */
    memcpy(doubled, orig_seq, orig_len);

    /* PADLENGTH spacer (DNA_N = 99) */
    memset(doubled + orig_len, 99, PADLENGTH);

    /* Reverse complement */
    for (int x = 0; x < orig_len; x++) {
        if (orig_seq[orig_len - 1 - x] == 99)
            doubled[orig_len + PADLENGTH + x] = 99;
        else
            doubled[orig_len + PADLENGTH + x] = 3 - orig_seq[orig_len - 1 - x];
    }

    *new_len = doubled_len;
    return doubled;
}

/* ================================================================
 * Public API: discover_default_params
 * ================================================================ */

void discover_default_params(DiscoverParams *params)
{
    params->l              = 0;        /* auto: ceil(1 + log4(N)) */
    params->L              = 10000;
    params->MAXOFFSET      = 5;
    params->MAXN           = 10000;
    params->MAXR           = 100000;
    params->MATCH          = 1;
    params->MISMATCH       = -1;
    params->GAP            = -5;
    params->CAPPENALTY     = -20;
    params->MINIMPROVEMENT = 3;
    params->WHEN_TO_STOP   = 100;
    params->MAXENTROPY     = -0.7f;
    params->GOODLENGTH     = 30;       /* plan: 50→30, MDL decides */
    params->MINTHRESH      = 2;        /* plan: 3→2, let MDL decide */
    params->TANDEMDIST     = 500;
    params->VERBOSE        = 0;
    params->freq_file      = NULL;
    params->freq_output    = NULL;
}

/* ================================================================
 * Dump frequency table to file (RepeatScout build_lmer_table format)
 *
 * Format: <lmer_string>\t<frequency>\t<position>\n
 * L-mers are canonicalized (lexicographically smaller of fwd/rc).
 * Sorted by frequency descending.
 * ================================================================ */

static void reverse_if_necessary(char *lmer_buf, int lmer_len)
{
    char rc;
    for (int x = 0; x < lmer_len; x++) {
        rc = 3 - lmer_buf[lmer_len - 1 - x];
        if (lmer_buf[x] < rc) return;
        if (lmer_buf[x] > rc) {
            /* Swap to RC */
            char tmp_buf[256];
            for (int y = 0; y < lmer_len; y++)
                tmp_buf[y] = 3 - lmer_buf[lmer_len - 1 - y];
            for (int y = 0; y < lmer_len; y++)
                lmer_buf[y] = tmp_buf[y];
            return;
        }
    }
    /* Palindrome — no change needed */
}

static int cmp_freq_desc(const void *a, const void *b)
{
    const int *fa = (const int *)a;
    const int *fb = (const int *)b;
    /* Each entry is 2 ints: [freq, occ]. Sort by freq descending. */
    if (fa[0] > fb[0]) return -1;
    if (fa[0] < fb[0]) return  1;
    return 0;
}

static void dump_freq_table(struct llist **hp, const char *path)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "discover: could not open %s for writing\n", path);
        return;
    }

    /* Count entries */
    int count = 0;
    for (int h = 0; h < HASH_SIZE; h++) {
        struct llist *tmp = hp[h];
        while (tmp) { count++; tmp = tmp->next; }
    }

    /* Collect (freq, occ) pairs */
    int (*entries)[2] = malloc((size_t)count * sizeof(*entries));
    if (!entries) {
        fprintf(stderr, "discover: out of memory for freq dump\n");
        fclose(fp);
        return;
    }

    int idx = 0;
    for (int h = 0; h < HASH_SIZE; h++) {
        struct llist *tmp = hp[h];
        while (tmp) {
            entries[idx][0] = tmp->freq;
            entries[idx][1] = tmp->lastocc;
            idx++;
            tmp = tmp->next;
        }
    }

    /* Sort by frequency descending */
    qsort(entries, (size_t)count, sizeof(*entries), cmp_freq_desc);

    /* Write */
    char lmer_buf[256];
    for (int i = 0; i < count; i++) {
        int occ = entries[i][1];
        for (int x = 0; x < l; x++)
            lmer_buf[x] = sequence[occ + x];
        reverse_if_necessary(lmer_buf, l);
        for (int x = 0; x < l; x++)
            fputc(num_to_char(lmer_buf[x]), fp);
        fprintf(fp, "\t%d\t%d\n", entries[i][0], entries[i][1]);
    }

    free(entries);
    fclose(fp);
    fprintf(stderr, "discover: wrote %d l-mers to %s\n", count, path);
}

/* ================================================================
 * Public API: discover_families
 * ================================================================ */

CandidateList *discover_families(const Genome *genome,
                                 const DiscoverParams *params)
{
    struct llist **headptr;
    struct llist *besttmp;
    int x;

    /* ---- Set parameters from DiscoverParams ---- */
    l              = params->l;
    L              = params->L;
    MAXOFFSET      = params->MAXOFFSET;
    MAXN           = params->MAXN;
    MAXR           = params->MAXR;
    MATCH          = params->MATCH;
    MISMATCH       = params->MISMATCH;
    GAP            = params->GAP;
    CAPPENALTY     = params->CAPPENALTY;
    MINIMPROVEMENT = params->MINIMPROVEMENT;
    WHEN_TO_STOP   = params->WHEN_TO_STOP;
    MAXENTROPY     = params->MAXENTROPY;
    GOODLENGTH     = params->GOODLENGTH;
    MINTHRESH      = params->MINTHRESH;
    TANDEMDIST     = params->TANDEMDIST;
    VERBOSE        = params->VERBOSE;

    /* ---- Setup genome data ---- */
    orig_length = (int)genome->length;
    disc_num_sequences = genome->num_sequences;

    /* Auto-compute l if not specified (use original length) */
    if (l <= 0)
        l = (int)ceil(1.0 + log((double)orig_length) / log(4.0));

    fprintf(stderr, "discover: l=%d, L=%d, MAXOFFSET=%d, MINTHRESH=%d, GOODLENGTH=%d\n",
            l, L, MAXOFFSET, MINTHRESH, GOODLENGTH);

    /* Build doubled genome: forward + PADLENGTH padding + RC */
    int doubled_len;
    sequence_owned = build_doubled_genome(genome->sequence, orig_length, &doubled_len);
    if (!sequence_owned) {
        fprintf(stderr, "discover: out of memory for doubled genome\n");
        return NULL;
    }
    sequence = sequence_owned;
    length = doubled_len;

    /* Copy boundaries to int array (RepeatScout uses int) */
    disc_boundaries = malloc((disc_num_sequences + 1) * sizeof(int));
    if (!disc_boundaries) {
        fprintf(stderr, "discover: out of memory for boundaries\n");
        free(sequence_owned); sequence_owned = NULL;
        return NULL;
    }
    for (int i = 0; i <= disc_num_sequences; i++)
        disc_boundaries[i] = (int)genome->boundaries[i];

    /* Allocate masking array (covers doubled genome) */
    removed = calloc(length, sizeof(char));
    if (!removed) {
        fprintf(stderr, "discover: out of memory for masking array\n");
        free(disc_boundaries);
        free(sequence_owned); sequence_owned = NULL;
        return NULL;
    }

    /* ---- Allocate workspace ---- */
    if (allocate_workspace() < 0) {
        fprintf(stderr, "discover: out of memory for workspace\n");
        free(removed);
        free(disc_boundaries);
        free(sequence_owned); sequence_owned = NULL;
        return NULL;
    }

    /* ---- Build hash table ---- */
    headptr = malloc(HASH_SIZE * sizeof(*headptr));
    if (!headptr) {
        fprintf(stderr, "discover: out of memory for headptr\n");
        free_workspace();
        return NULL;
    }

    if (params->freq_file) {
        fprintf(stderr, "discover: reading frequency table from %s\n", params->freq_file);
        build_headptr_from_freq(headptr, params->freq_file);
    } else {
        fprintf(stderr, "discover: counting l-mers internally...\n");
        build_headptr_internal(headptr);
    }

    /* Dump frequency table if requested (before trim, same as build_lmer_table) */
    if (params->freq_output) {
        dump_freq_table(headptr, params->freq_output);
    }

    trim_headptr(headptr);
    build_all_pos(headptr);

    /* ---- Allocate result ---- */
    CandidateList *result = malloc(sizeof(CandidateList));
    if (!result) {
        fprintf(stderr, "discover: out of memory for result\n");
        free_headptr(headptr);
        free_workspace();
        return NULL;
    }
    result->cap_families = 1024;
    result->families = malloc(result->cap_families * sizeof(CandidateFamily));
    if (!result->families) {
        fprintf(stderr, "discover: out of memory for families\n");
        free(result);
        free_headptr(headptr);
        free_workspace();
        return NULL;
    }
    result->num_families = 0;

    /* ---- Main discovery loop ---- */
    R = 0;
    prevbestfreq = 1000000000;
    prevbesthash = 0;

    while (1) {
        besttmp = find_besttmp(headptr);

        if (besttmp == NULL || besttmp->freq < MINTHRESH) {
            if (VERBOSE)
                fprintf(stderr, "discover: stopped at R=%d (no more frequent l-mers)\n", R);
            break;
        }

        /* Initialize seed in master */
        for (x = 0; x < l; x++)
            master[L + x] = sequence[besttmp->lastocc + x];

        /* Build occurrence list */
        build_pos(besttmp);

        if (N < MINTHRESH) {
            if (VERBOSE >= 2)
                fprintf(stderr, "discover: R=%d N=%d < MINTHRESH=%d, updating freq\n",
                        R, N, MINTHRESH);
            besttmp->freq = N;
            continue;
        }

        /* Initialize offset tracking */
        for (int n = 0; n < N; n++) {
            best_left_offset[n] = -1;
            save_left_offset[n] = -1;
            best_right_offset[n] = -1;
            save_right_offset[n] = -1;
        }

        /* Allocate masters[R] if needed */
        if (masters_allocated[R] == 0) {
            masters[R] = malloc((2 * L + l) * sizeof(char));
            if (!masters[R]) {
                fprintf(stderr, "discover: out of memory for masters[%d]\n", R);
                break;
            }
            masters_allocated[R] = 1;
        }

        /* Extend right and left */
        extend_right();
        extend_left();

        int cons_len = masterend[R] - masterstart[R] + 1;

        if (cons_len >= GOODLENGTH) {
            /* Good family: save and collect instances */
            if (VERBOSE) {
                fprintf(stderr, "discover: R=%d N=%d length=%d\n", R, N, cons_len);
            }

            /* Grow result array if needed */
            if (result->num_families >= result->cap_families) {
                result->cap_families *= 2;
                CandidateFamily *tmp_f = realloc(result->families,
                    result->cap_families * sizeof(CandidateFamily));
                if (!tmp_f) {
                    fprintf(stderr, "discover: out of memory growing families\n");
                    break;
                }
                result->families = tmp_f;
            }

            /* Create CandidateFamily */
            CandidateFamily *fam = &result->families[result->num_families];
            memset(fam, 0, sizeof(*fam));
            fam->id = (uid_t)result->num_families;

            /* Copy consensus */
            fam->consensus = malloc(cons_len * sizeof(char));
            if (!fam->consensus) {
                fprintf(stderr, "discover: out of memory for consensus\n");
                break;
            }
            for (x = 0; x < cons_len; x++)
                fam->consensus[x] = masters[R][masterstart[R] + x];
            fam->consensus_length = cons_len;

            /* Set default fields */
            fam->topology = TOPO_LINEAR;
            fam->component_id = 0;
            fam->estimated_copies = N;

            /* Collect instances from extension occurrences */
            collect_instances_from_extension(fam);

            result->num_families++;

            R++;
            if (R == MAXR) break;

            /* Mask this family */
            mask_headptr(headptr);
        } else {
            /* Short family: mask but don't record */
            R++;
            if (R == MAXR) break;
            mask_headptr(headptr);
            R--;
        }
    }

    fprintf(stderr, "discover: found %d families (R=%d total including short)\n",
            result->num_families, R);

    /* ---- Cleanup ---- */
    free_headptr(headptr);
    free_workspace();

    return result;
}
