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
 * Thread-safe: all mutable state lives in DiscoverContext (heap-allocated).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "discover_internal.h"
#include "kmer.h"

/* ================================================================
 * Cross-module helpers (decls in discover_internal.h).  These are
 * NOT static — discover_mask.c calls them.
 * ================================================================ */

char compl_base(char c)
{
    if (c == DNA_N) return DNA_N;
    return 3 - c;
}

/* ================================================================
 * Hash functions (symmetric w.r.t. reverse complement)
 * ================================================================ */

int hash_function(DiscoverContext *C, const char *lmer)
{
    int x, ans, ans2;

    for (x = 0; x < C->l; x++)
        if (lmer[x] == DNA_N) return -1;

    ans = 0;
    for (x = 0; x < C->l; x++)
        ans = (4 * ans + (lmer[x] % 4)) % (int)C->hash_size;
    ans2 = 0;
    for (x = 0; x < C->l; x++)
        ans2 = (4 * ans2 + ((3 - lmer[C->l - 1 - x]) % 4)) % (int)C->hash_size;
    if (ans2 > ans) ans = ans2;
    return ans;
}

int smallhash_function(const char *lmer)
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
 * L-mer matching functions (cross-module — see discover_internal.h)
 * ================================================================ */

int lmermatch(DiscoverContext *C, const char *lmer1, const char *lmer2)
{
    for (int x = 0; x < C->l; x++)
        if (lmer1[x] != lmer2[x]) return 0;
    return 1;
}

int lmermatchrc(DiscoverContext *C, const char *lmer1, const char *lmer2)
{
    for (int x = 0; x < C->l; x++)
        if (lmer1[x] + lmer2[C->l - 1 - x] != 3) return 0;
    return 1;
}

int lmermatcheither(DiscoverContext *C, const char *lmer1, const char *lmer2)
{
    for (int x = 0; x < C->l; x++)
        if (lmer1[x] != lmer2[x]) return lmermatchrc(C, lmer1, lmer2);
    return 1;
}

/* ================================================================
 * Entropy (Shannon, natural log — returns negative value)
 * ================================================================ */

static double compute_entropy(DiscoverContext *C, const char *lmer)
{
    int count[4] = {0, 0, 0, 0};
    double answer = 0.0;

    for (int x = 0; x < C->l; x++)
        count[(int)lmer[x]] += 1;

    for (int x = 0; x < 4; x++) {
        if (count[x] == 0) continue;
        double y = (double)count[x] / (double)C->l;
        answer += y * log(y);
    }
    return answer;
}

/* ================================================================
 * Periodicity filter — reject l-mers with short tandem period
 * ================================================================
 *
 * The Shannon entropy filter (compute_entropy) catches l-mers with
 * unbalanced base composition such as (A)n or (CG)n.  It does NOT
 * catch periodic l-mers with balanced composition such as (GAGA)n,
 * (ATCATC)n, or (CAG)n: their per-base distribution can be uniform
 * yet the sequence is a simple tandem repeat that should not seed
 * a discovery.  Without this guard, simple-repeat seeds dominate
 * the high-frequency end of the seed pool and produce dozens of
 * spurious candidate families.
 *
 * Returns the smallest period in [1, l/2] for which the l-mer is
 * >= PERIODIC_MATCH_PCT identical to itself shifted by that period,
 * or 0 if no such period exists.
 */
#define PERIODIC_MATCH_PCT 85

static int is_periodic_lmer(const char *lmer, int l)
{
    if (l < 4) return 0;
    for (int p = 1; p <= l / 2; p++) {
        int compares = l - p;
        if (compares < 4) continue;
        int matches = 0;
        for (int i = 0; i < compares; i++)
            if (lmer[i] == lmer[i + p]) matches++;
        if (matches * 100 >= compares * PERIODIC_MATCH_PCT) return p;
    }
    return 0;
}

/* ================================================================
 * Sequence boundary helper (FIXED — no hardcoded override)
 * ================================================================ */

static void get_boundaries(DiscoverContext *C, int n, gpos_t *bStart_out, gpos_t *bEnd_out)
{
    gpos_t bStart = PADLENGTH;
    gpos_t bEnd;

    if (C->upperBoundI[n] == -1) {
        /* Defensive: position past all boundaries */
        bStart = PADLENGTH;
        bEnd = C->length;
    } else {
        if (C->upperBoundI[n] > 0)
            bStart = C->disc_boundaries[C->upperBoundI[n] - 1] + PADLENGTH;
        bEnd = C->disc_boundaries[C->upperBoundI[n]] + PADLENGTH;
    }

    /* RepeatScout bug: bStart=0; bEnd=50000000; — NOT applied here */

    *bStart_out = bStart;
    *bEnd_out = bEnd;
}

/* ================================================================
 * Find sequence index for a padded position
 * ================================================================ */

static int find_seq_index(DiscoverContext *C, gpos_t padded_pos)
{
    gpos_t raw = padded_pos - PADLENGTH;
    if (raw < 0) return -1;

    for (int i = 0; i < C->disc_num_sequences; i++) {
        if (i == 0) {
            if (raw < C->disc_boundaries[0]) return 0;
        } else {
            if (raw >= C->disc_boundaries[i - 1] && raw < C->disc_boundaries[i])
                return i;
        }
    }
    return C->disc_num_sequences - 1;
}

/* ================================================================
 * Hash table: read from .freq file (RepeatScout format)
 * ================================================================ */

static void build_headptr_from_freq(DiscoverContext *C, struct llist **headptr, const char *freq_file)
{
    FILE *fp;
    char string[1000];
    int thisfreq, x, h;
    gpos_t thisocc;
    struct llist *tmp;

    fp = fopen(freq_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "discover: could not open freq file %s\n", freq_file);
        exit(1);
    }

    for (h = 0; h < C->hash_size; h++)
        headptr[h] = NULL;

    while (1) {
        if (fscanf(fp, "%s", string) != 1) break;
        if (fscanf(fp, "%d", &thisfreq) != 1) {
            fprintf(stderr, "discover: error reading frequency from %s\n", freq_file);
            exit(1);
        }
        if (fscanf(fp, "%" SCNd64, &thisocc) != 1) {
            fprintf(stderr, "discover: error reading occurrence from %s\n", freq_file);
            exit(1);
        }

        if ((int)strlen(string) != C->l) {
            fprintf(stderr, "discover: l-mer length mismatch: expected %d, got %d\n",
                    C->l, (int)strlen(string));
            exit(1);
        }

        for (x = 0; x < C->l; x++)
            string[x] = char_to_num(string[x]);

        h = hash_function(C, string);
        if (h < 0) continue;

        /* Check if already present (duplicate = end of useful data) */
        tmp = headptr[h];
        while (tmp != NULL) {
            if (lmermatcheither(C, C->sequence + tmp->lastocc, string) == 1)
                break;
            tmp = tmp->next;
        }
        if (tmp != NULL) break; /* already in table — done */

        /* Entropy filter */
        if (compute_entropy(C, C->sequence + thisocc) > C->MAXENTROPY) continue;
        /* Periodicity filter — reject simple-repeat seeds (e.g. GAGAGA) */
        if (is_periodic_lmer(C->sequence + thisocc, C->l)) continue;

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

static void build_headptr_internal(DiscoverContext *C, struct llist **headptr)
{
    int h;
    gpos_t x;
    struct llist *tmp;

    for (h = 0; h < C->hash_size; h++)
        headptr[h] = NULL;

    /*
     * Count l-mer frequencies with per-strand TANDEMDIST filtering,
     * matching RepeatScout's build_lmer_table behavior.
     * Only scan the forward copy (not the RC copy in doubled genome).
     */
    gpos_t count_end = C->orig_length - C->l;
    for (x = PADLENGTH; x <= count_end; x++) {
        h = hash_function(C, C->sequence + x);
        if (h < 0) continue;

        /* Check if already in table */
        tmp = headptr[h];
        while (tmp != NULL) {
            if (lmermatch(C, C->sequence + tmp->lastplusocc, C->sequence + x)) {
                /* Forward match: TANDEMDIST check against lastplusocc */
                if (x - tmp->lastplusocc >= C->TANDEMDIST)
                    tmp->freq++;
                tmp->lastplusocc = x;
                tmp->lastocc = x;
                break;
            } else if (lmermatchrc(C, C->sequence + tmp->lastplusocc, C->sequence + x)) {
                /* Reverse complement match: TANDEMDIST check against lastminusocc */
                if (x - tmp->lastminusocc >= C->TANDEMDIST)
                    tmp->freq++;
                tmp->lastminusocc = x;
                tmp->lastocc = x;
                break;
            }
            tmp = tmp->next;
        }

        if (tmp == NULL) {
            /* New entry — apply entropy + periodicity filters */
            if (compute_entropy(C, C->sequence + x) > C->MAXENTROPY) continue;
            if (is_periodic_lmer(C->sequence + x, C->l)) continue;

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
 * Hash table: parallel l-mer counting via kmer.c striped-lock hash
 *
 * Option γ adapter: replaces the O(N) serial scan in
 * build_headptr_internal with a parallel kmer_count() call, then
 * translates surviving KmerEntry records into the discover.c llist
 * structure expected by the rest of the pipeline.
 *
 * Scan range matches build_headptr_internal: only the forward copy
 * of the doubled genome ([PADLENGTH, orig_length)).
 * ================================================================ */

static void build_headptr_parallel(DiscoverContext *C,
                                   struct llist **headptr,
                                   int num_threads)
{
    int h;

    for (h = 0; h < C->hash_size; h++)
        headptr[h] = NULL;

    /* Build a temporary Genome wrapping the forward (non-doubled) copy.
     * kmer.c's kmer_count scans [PADLENGTH, length - k + 1), so length
     * must be C->orig_length (NOT C->length which covers the RC copy). */
    Genome tmp_genome;
    memset(&tmp_genome, 0, sizeof(tmp_genome));
    tmp_genome.sequence      = C->sequence_owned;
    tmp_genome.length        = C->orig_length;
    tmp_genome.raw_length    = C->orig_length - PADLENGTH;
    tmp_genome.boundaries    = C->disc_boundaries;
    tmp_genome.num_sequences = C->disc_num_sequences;
    tmp_genome.sequence_ids  = NULL;

    KmerTable *kt = kmer_count(&tmp_genome, C->l, C->TANDEMDIST, num_threads);
    if (kt == NULL) {
        fprintf(stderr, "discover: parallel kmer_count failed\n");
        exit(1);
    }

    /* Decode buffer for unpacking 2-bit kmers back to 0-3 char[] */
    char decoded[64];
    if (C->l > 64) {
        fprintf(stderr, "discover: l=%d exceeds decode buffer (64)\n", C->l);
        kmer_free(kt);
        exit(1);
    }

    for (size_t b = 0; b < kt->table_size; b++) {
        KmerEntry *e = kt->buckets[b];
        while (e != NULL) {
            /* Skip below-threshold early to avoid pointless decode/alloc */
            if (e->frequency < C->MINTHRESH) {
                e = e->next;
                continue;
            }

            /* Decode packed canonical kmer back to char[l] (0-3 encoding).
             * MSB-first packing means the highest 2 bits encode index 0;
             * unpack from index l-1 downward. Use a local copy so we
             * don't mutate the entry. */
            uint64_t k = e->kmer;
            for (int i = C->l - 1; i >= 0; i--) {
                decoded[i] = (char)(k & 3);
                k >>= 2;
            }

            /* Apply discover-side filters that kmer.c does not know about.
             * compute_entropy and is_periodic_lmer accept char[] in
             * 0-3 encoding, exactly what we just produced. */
            if (compute_entropy(C, decoded) > C->MAXENTROPY) {
                e = e->next;
                continue;
            }
            if (is_periodic_lmer(decoded, C->l)) {
                e = e->next;
                continue;
            }

            int hh = hash_function(C, decoded);
            if (hh < 0) {
                e = e->next;
                continue;
            }

            /* Defensive: skip if already present in headptr[hh].
             * kmer.c canonicalises and deduplicates, so this should
             * not happen, but a stray hash collision via the asymmetric
             * canonical forms could in principle produce one. */
            struct llist *existing = headptr[hh];
            int dup = 0;
            while (existing != NULL) {
                if (lmermatcheither(C, C->sequence + existing->lastocc, decoded)) {
                    dup = 1;
                    break;
                }
                existing = existing->next;
            }
            if (dup) {
                e = e->next;
                continue;
            }

            /* Choose lastocc: prefer plus, fall back to minus if plus is
             * the sentinel (-1000000), in which case the kmer was only
             * seen on the RC strand. */
            gpos_t lastocc;
            if (e->last_plus_occ >= PADLENGTH)
                lastocc = e->last_plus_occ;
            else
                lastocc = e->last_minus_occ;

            struct llist *node = malloc(sizeof(*node));
            if (node == NULL) {
                fprintf(stderr, "discover: out of memory in build_headptr_parallel\n");
                kmer_free(kt);
                exit(1);
            }
            node->freq         = e->frequency;
            node->lastocc      = lastocc;
            node->lastplusocc  = e->last_plus_occ;
            node->lastminusocc = e->last_minus_occ;
            node->pos          = NULL;
            node->next         = headptr[hh];
            headptr[hh] = node;

            e = e->next;
        }
    }

    kmer_free(kt);
}

/* ================================================================
 * Trim: remove l-mers with freq < MINTHRESH, reset freq to 0
 * ================================================================ */

static void trim_headptr(DiscoverContext *C, struct llist **headptr)
{
    int h;
    struct llist *tmp, *prevtmp, *nexttmp;

    for (h = 0; h < C->hash_size; h++) {
        prevtmp = NULL;
        tmp = headptr[h];
        while (tmp != NULL) {
            if (tmp->freq >= C->MINTHRESH) {
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

static void build_all_pos(DiscoverContext *C, struct llist **headptr)
{
    gpos_t x, pos1, pos2;
    int h;
    struct llist *tmp;
    struct posllist *postmp, *postmp2, *prevpostmp, *nextpostmp;
    int currBoundary = 0;

    if (C->VERBOSE) fprintf(stderr, "discover: building all positions...\n");

    /* Pass 1: scan genome, add positions to matching l-mers */
    for (x = C->l - 1; x < C->length - C->l + 1; x++) {
        /* Check sequence boundaries */
        if (C->disc_boundaries[currBoundary] > 0) {
            if (x == C->disc_boundaries[currBoundary] + PADLENGTH)
                currBoundary++;
            if (x + C->l > C->disc_boundaries[currBoundary] + PADLENGTH)
                continue;
        }

        h = hash_function(C, C->sequence + x);
        if (h < 0) continue;

        tmp = headptr[h];
        while (tmp != NULL) {
            if (lmermatcheither(C, C->sequence + tmp->lastocc, C->sequence + x)) {
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

    if (C->VERBOSE) fprintf(stderr, "discover: TANDEMDIST filtering...\n");

    /* Pass 2: remove position pairs within TANDEMDIST (same strand only) */
    for (h = 0; h < C->hash_size; h++) {
        tmp = headptr[h];
        while (tmp != NULL) {
            /* Mark within-TANDEMDIST same-strand pairs as negative */
            postmp = tmp->pos;
            while (postmp != NULL && postmp->next != NULL) {
                pos1 = postmp->this;
                postmp2 = postmp;
                while ((postmp2 = postmp2->next) != NULL) {
                    pos2 = postmp2->this;
                    if (pos1 - pos2 >= C->TANDEMDIST)
                        break;
                    if (lmermatch(C, C->sequence + pos1, C->sequence + pos2)) {
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
                /* Remove (TANDEMDIST-pruned) — use permanent sentinel
                 * so this position never becomes available to any
                 * future seed regardless of MAX_FAMILY_CLAIMS. */
                nextpostmp = postmp->next;
                tmp->freq--;
                C->removed[-postmp->this] = CLAIM_PERMANENT;
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

    if (C->VERBOSE) fprintf(stderr, "discover: position building complete\n");
}

/* ================================================================
 * Seed selection: find highest-frequency l-mer
 * (Exact port of RepeatScout find_besttmp, lines 748-788)
 * ================================================================ */

static struct llist *find_besttmp(DiscoverContext *C, struct llist **headptr)
{
    int h;
    struct llist *tmp, *besttmp;
    int bestfreq;

    /* Try to match prevbestfreq first (locality optimization) */
    for (h = C->prevbesthash; h < C->hash_size; h++) {
        tmp = headptr[h];
        while (tmp != NULL) {
            if (tmp->freq == C->prevbestfreq) {
                C->prevbesthash = h;
                return tmp;
            }
            tmp = tmp->next;
        }
    }

    /* Global search for max freq */
    besttmp = NULL;
    bestfreq = 0;
    for (h = 0; h < C->hash_size; h++) {
        tmp = headptr[h];
        while (tmp != NULL) {
            if (tmp->freq > bestfreq) {
                besttmp = tmp;
                bestfreq = tmp->freq;
                C->prevbesthash = h;
            }
            tmp = tmp->next;
        }
    }
    C->prevbestfreq = bestfreq;
    return besttmp;
}

/* ================================================================
 * Build occurrence list for a seed
 * (Exact port of RepeatScout build_pos, lines 996-1032)
 * ================================================================ */

static void build_pos(DiscoverContext *C, struct llist *besttmp)
{
    struct posllist *postmp;
    int x;

    postmp = besttmp->pos;
    C->N = 0;
    while (postmp != NULL) {
        /* Coverage-aware: allow positions claimed by < MAX_FAMILY_CLAIMS
         * families.  CLAIM_PERMANENT sentinel keeps TANDEMDIST-pruned
         * positions excluded. */
        if (C->removed[postmp->this] < MAX_FAMILY_CLAIMS) {
            C->pos[C->N] = postmp->this;

            /* Strand determination: lmermatch = forward, else reverse */
            if (lmermatch(C, C->sequence + besttmp->lastocc, C->sequence + postmp->this))
                C->rev[C->N] = 0;
            else
                C->rev[C->N] = 1;

            /* Find upper boundary index */
            C->upperBoundI[C->N] = -1;
            x = 0;
            while (C->disc_boundaries[x] != 0) {
                if (C->disc_boundaries[x] + PADLENGTH > C->pos[C->N]) {
                    C->upperBoundI[C->N] = x;
                    break;
                }
                x++;
            }

            C->N++;
            if (C->N == C->MAXN) break;
        }
        postmp = postmp->next;
    }
}

/* ================================================================
 * DP scoring: compute_score_right
 * (Exact port of RepeatScout lines 1877-1969, with boundary fix)
 * ================================================================ */

static int compute_score_right(DiscoverContext *C, int y, int n, int offset, char a)
{
    int oldoffset, tempscore, ans, ismatch, x;
    gpos_t bStart, bEnd;

    get_boundaries(C, n, &bStart, &bEnd);

    /* Boundary check */
    if (C->rev[n]) {
        if (C->pos[n] - (offset + y - C->L - C->l) - 1 < bStart)
            return 0;
    } else {
        if (C->pos[n] + offset + y - C->L >= bEnd)
            return 0;
    }

    ans = -1000000000;

    /* Case A: gap in sequence (oldoffset = offset+1) */
    if (offset < C->MAXOFFSET) {
        oldoffset = offset + 1;
        tempscore = C->score[(y - 1) % 2][n][oldoffset + C->MAXOFFSET] + C->GAP;
        if (tempscore > ans) ans = tempscore;
    }

    /* Case B: diagonal match/mismatch (oldoffset = offset) */
    oldoffset = offset;
    tempscore = C->score[(y - 1) % 2][n][oldoffset + C->MAXOFFSET];
    if (C->rev[n]) {
        if (a == compl_base(C->sequence[C->pos[n] - (offset + y - C->L - C->l) - 1]))
            tempscore += C->MATCH;
        else
            tempscore += C->MISMATCH;
    } else {
        if (a == C->sequence[C->pos[n] + offset + y - C->L])
            tempscore += C->MATCH;
        else
            tempscore += C->MISMATCH;
    }
    if (tempscore > ans) ans = tempscore;

    /* Case C: multiple gaps in consensus (oldoffset < offset) */
    for (oldoffset = -C->MAXOFFSET; oldoffset < offset; oldoffset++) {
        ismatch = 0;
        for (x = oldoffset; x <= offset; x++) {
            if (C->rev[n]) {
                if (a == compl_base(C->sequence[C->pos[n] - (x + y - C->L - C->l) - 1]))
                    ismatch = 1;
            } else {
                if (a == C->sequence[C->pos[n] + x + y - C->L])
                    ismatch = 1;
            }
        }
        tempscore = C->score[(y - 1) % 2][n][oldoffset + C->MAXOFFSET];
        tempscore += (offset - oldoffset) * C->GAP;
        if (ismatch) tempscore += C->MATCH;
        else tempscore += C->MISMATCH;
        if (tempscore > ans) ans = tempscore;
    }

    return ans;
}

/* ================================================================
 * DP scoring: compute_score_left
 * (Exact port of RepeatScout lines 1976-2069, with boundary fix)
 * ================================================================ */

static int compute_score_left(DiscoverContext *C, int w, int n, int offset, char a)
{
    int oldoffset, tempscore, ans, ismatch, x;
    gpos_t bStart, bEnd;

    get_boundaries(C, n, &bStart, &bEnd);

    /* Boundary check (note: directions are swapped vs right) */
    if (C->rev[n]) {
        if (C->pos[n] - (offset + w - C->L - C->l) - 1 >= bEnd)
            return 0;
    } else {
        if (C->pos[n] + offset + w - C->L < bStart)
            return 0;
    }

    ans = -1000000000;

    /* Case A: gap in sequence (oldoffset = offset-1) */
    if (offset > -C->MAXOFFSET) {
        oldoffset = offset - 1;
        tempscore = C->score[(w + 1) % 2][n][oldoffset + C->MAXOFFSET] + C->GAP;
        if (tempscore > ans) ans = tempscore;
    }

    /* Case B: diagonal (oldoffset = offset) */
    oldoffset = offset;
    tempscore = C->score[(w + 1) % 2][n][oldoffset + C->MAXOFFSET];
    if (C->rev[n]) {
        if (a == compl_base(C->sequence[C->pos[n] - (offset + w - C->L - C->l) - 1]))
            tempscore += C->MATCH;
        else
            tempscore += C->MISMATCH;
    } else {
        if (a == C->sequence[C->pos[n] + offset + w - C->L])
            tempscore += C->MATCH;
        else
            tempscore += C->MISMATCH;
    }
    if (tempscore > ans) ans = tempscore;

    /* Case C: multiple gaps in consensus (oldoffset > offset) */
    for (oldoffset = offset + 1; oldoffset <= C->MAXOFFSET; oldoffset++) {
        ismatch = 0;
        for (x = offset; x <= oldoffset; x++) {
            if (C->rev[n]) {
                if (a == compl_base(C->sequence[C->pos[n] - (x + w - C->L - C->l) - 1]))
                    ismatch = 1;
            } else {
                if (a == C->sequence[C->pos[n] + x + w - C->L])
                    ismatch = 1;
            }
        }
        tempscore = C->score[(w + 1) % 2][n][oldoffset + C->MAXOFFSET];
        tempscore += (oldoffset - offset) * C->GAP;
        if (ismatch) tempscore += C->MATCH;
        else tempscore += C->MISMATCH;
        if (tempscore > ans) ans = tempscore;
    }

    return ans;
}

/* ================================================================
 * Checkpoint: compute_totalbestscore_right
 * (Exact port of RepeatScout lines 1722-1776)
 * ================================================================ */

static void compute_totalbestscore_right(DiscoverContext *C, int y)
{
    int n, bestscore, offset;

    C->nrepeatocc = 0;
    C->nactiverepeatocc = 0;
    C->totalbestscore = 0;

    for (n = 0; n < C->N; n++) {
        bestscore = C->bestbestscore[n] + C->CAPPENALTY;
        if (bestscore < 0) bestscore = 0;

        int bestscore_offset = -100000;
        for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
            if (C->score[y % 2][n][offset + C->MAXOFFSET] > bestscore) {
                bestscore = C->score[y % 2][n][offset + C->MAXOFFSET];
                bestscore_offset = offset;
            }
        }

        if (bestscore > 0) C->nrepeatocc++;
        if (bestscore > C->bestbestscore[n] + C->CAPPENALTY) C->nactiverepeatocc++;
        if (bestscore > C->bestbestscore[n]) {
            if (bestscore_offset >= -C->MAXOFFSET)
                C->best_right_offset[n] = bestscore_offset + y;
            C->bestbestscore[n] = bestscore;
        }
        C->totalbestscore += bestscore;
    }

    /* MINIMPROVEMENT checkpoint */
    if ((C->totalbestscore >= C->besttotalbestscore + (y - C->besty) * C->MINIMPROVEMENT)
        && (C->totalbestscore > C->besttotalbestscore)) {
        C->besty = y;
        C->besttotalbestscore = C->totalbestscore;
        C->bestnrepeatocc = C->nrepeatocc;
        C->bestnactiverepeatocc = C->nactiverepeatocc;
        for (n = 0; n < C->N; n++) {
            C->save_right_offset[n] = C->best_right_offset[n];
            C->savebestscore[n] = C->bestbestscore[n];
            for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++)
                C->score_of_besty[n][offset + C->MAXOFFSET] = C->score[y % 2][n][offset + C->MAXOFFSET];
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

static void compute_totalbestscore_left(DiscoverContext *C, int w)
{
    int n, bestscore, offset;

    C->nrepeatocc = 0;
    C->nactiverepeatocc = 0;
    C->totalbestscore = 0;

    for (n = 0; n < C->N; n++) {
        int bestscore_offset = -100000;
        bestscore = C->bestbestscore[n] + C->CAPPENALTY;
        if (bestscore < 0) bestscore = 0;

        for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
            if (C->score[w % 2][n][offset + C->MAXOFFSET] > bestscore) {
                bestscore = C->score[w % 2][n][offset + C->MAXOFFSET];
                bestscore_offset = offset;
            }
        }

        if (bestscore > 0) C->nrepeatocc++;
        if (bestscore > C->bestbestscore[n] + C->CAPPENALTY) C->nactiverepeatocc++;
        if (bestscore > C->bestbestscore[n]) {
            if (bestscore_offset >= -C->MAXOFFSET)
                C->best_left_offset[n] = bestscore_offset + w;
            C->bestbestscore[n] = bestscore;
        }
        C->totalbestscore += bestscore;
    }

    /* MINIMPROVEMENT checkpoint — only save offsets, NOT scores/score_of_besty */
    if ((C->totalbestscore >= C->besttotalbestscore + (C->bestw - w) * C->MINIMPROVEMENT)
        && (C->totalbestscore > C->besttotalbestscore)) {
        C->bestw = w;
        C->besttotalbestscore = C->totalbestscore;
        C->bestnrepeatocc = C->nrepeatocc;
        C->bestnactiverepeatocc = C->nactiverepeatocc;
        for (n = 0; n < C->N; n++)
            C->save_left_offset[n] = C->best_left_offset[n];
    }
}

/* ================================================================
 * extend_right
 * (Exact port of RepeatScout lines 1034-1117)
 * ================================================================ */

static void extend_right(DiscoverContext *C)
{
    int y, n, bestscore, tempscore, offset;
    char a, besta;
    int newtotalbestscore, newtotalbestscore_a;
    besta = 0;

    /* Initialize score at seed's last base */
    y = C->L + C->l - 1;
    for (n = 0; n < C->N; n++) {
        for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
            C->score[y % 2][n][offset + C->MAXOFFSET] = C->l * C->MATCH;
            if (offset < 0)
                C->score[y % 2][n][offset + C->MAXOFFSET] += -offset * C->GAP;
            if (offset > 0)
                C->score[y % 2][n][offset + C->MAXOFFSET] += offset * C->GAP;
        }
        C->bestbestscore[n] = C->l * C->MATCH;
    }

    /* Initialize checkpoint */
    C->besty = C->L + C->l - 1;
    for (n = 0; n < C->N; n++) {
        for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++)
            C->score_of_besty[n][offset + C->MAXOFFSET] = C->score[y % 2][n][offset + C->MAXOFFSET];
    }
    C->besttotalbestscore = 0;
    compute_totalbestscore_right(C, y);

    /* Extend right, one base at a time */
    for (y = C->L + C->l; y < 2 * C->L + C->l; y++) {
        newtotalbestscore = 0;
        for (a = 0; a < 4; a++) {
            newtotalbestscore_a = 0;
            for (n = 0; n < C->N; n++) {
                bestscore = C->bestbestscore[n] + C->CAPPENALTY;
                if (bestscore < 0) bestscore = 0;
                for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
                    tempscore = compute_score_right(C, y, n, offset, a);
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

        C->master[y] = besta;
        for (n = 0; n < C->N; n++) {
            for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++)
                C->score[y % 2][n][offset + C->MAXOFFSET] = compute_score_right(C, y, n, offset, besta);
        }
        compute_totalbestscore_right(C, y);

        /* Adaptive stop (思路 4): give long consensi a longer quiet
         * window proportional to extension already achieved.  Without
         * this, long elements with internal "信息低谷" (low-info
         * 100-300 bp stretches) trigger early stop even though the
         * extension has already produced a meaningful long consensus. */
        int extended_so_far = C->besty - (C->L + C->l - 1);
        int adaptive_when_to_stop = C->WHEN_TO_STOP;
        if (extended_so_far / 10 > adaptive_when_to_stop)
            adaptive_when_to_stop = extended_so_far / 10;
        if (y - C->besty >= adaptive_when_to_stop) break;
    }

    if (y == 2 * C->L + C->l && C->VERBOSE)
        fprintf(stderr, "Warning: extended right all the way to %d\n", y);

    y = C->besty;
    C->totalbestscore = C->besttotalbestscore;
    C->nrepeatocc = C->bestnrepeatocc;
    C->nactiverepeatocc = C->bestnactiverepeatocc;
}

/* ================================================================
 * extend_left
 * (Exact port of RepeatScout lines 1119-1205)
 * ================================================================ */

static void extend_left(DiscoverContext *C)
{
    int w, n, bestscore, tempscore, offset;
    char a, besta;
    int newtotalbestscore, newtotalbestscore_a;
    besta = 0;

    /* Initialize from right extension checkpoint */
    w = C->L;
    for (n = 0; n < C->N; n++) {
        bestscore = C->savebestscore[n] + C->CAPPENALTY;
        for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
            if (C->score_of_besty[n][offset + C->MAXOFFSET] > bestscore)
                bestscore = C->score_of_besty[n][offset + C->MAXOFFSET];
        }
        for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
            C->score[w % 2][n][offset + C->MAXOFFSET] = bestscore;
            if (offset < 0)
                C->score[w % 2][n][offset + C->MAXOFFSET] += -offset * C->GAP;
            if (offset > 0)
                C->score[w % 2][n][offset + C->MAXOFFSET] += offset * C->GAP;
        }
        C->bestbestscore[n] = bestscore;
    }

    /* Initialize checkpoint */
    C->bestw = C->L;
    compute_totalbestscore_left(C, w);

    /* Extend left */
    for (w = C->L - 1; w >= 0; w--) {
        newtotalbestscore = 0;
        for (a = 0; a < 4; a++) {
            newtotalbestscore_a = 0;
            for (n = 0; n < C->N; n++) {
                bestscore = C->bestbestscore[n] + C->CAPPENALTY;
                if (bestscore < 0) bestscore = 0;
                for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++) {
                    tempscore = compute_score_left(C, w, n, offset, a);
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

        C->master[w] = besta;
        for (n = 0; n < C->N; n++) {
            for (offset = -C->MAXOFFSET; offset <= C->MAXOFFSET; offset++)
                C->score[w % 2][n][offset + C->MAXOFFSET] = compute_score_left(C, w, n, offset, besta);
        }
        compute_totalbestscore_left(C, w);

        /* Adaptive stop (思路 4): mirror of extend_right. */
        int extended_left_so_far = C->L - C->bestw;
        int adaptive_when_to_stop_l = C->WHEN_TO_STOP;
        if (extended_left_so_far / 10 > adaptive_when_to_stop_l)
            adaptive_when_to_stop_l = extended_left_so_far / 10;
        if (w - C->bestw <= -adaptive_when_to_stop_l) break;
    }

    if (w == -1 && C->VERBOSE)
        fprintf(stderr, "Warning: extended left all the way to %d\n", w);

    w = C->bestw;
    C->totalbestscore = C->besttotalbestscore;
    C->nrepeatocc = C->bestnrepeatocc;
    C->nactiverepeatocc = C->bestnactiverepeatocc;

    C->masterstart[C->R] = C->bestw;
    C->masterend[C->R] = C->besty;
    for (w = C->bestw; w <= C->besty; w++)
        C->masters[C->R][w] = C->master[w];
}

/* ================================================================
 * Masking subsystem moved to discover_mask.c (M3#8 split).
 * mask_headptr() is the single public entry; see discover_internal.h.
 * ================================================================ */


/* ================================================================
 * Workspace allocation and deallocation
 * ================================================================ */

static int allocate_workspace(DiscoverContext *C)
{
    int x, n;

    C->master = malloc((2 * C->L + C->l + 1) * sizeof(char));
    if (!C->master) return -1;
    C->master[2 * C->L + C->l] = '\0';

    C->masters_allocated = calloc(C->MAXR, sizeof(int));
    C->masterstart = malloc(C->MAXR * sizeof(int));
    C->masterend = malloc(C->MAXR * sizeof(int));
    if (!C->masters_allocated || !C->masterstart || !C->masterend) return -1;

    C->masters = malloc(C->MAXR * sizeof(*C->masters));
    if (!C->masters) return -1;

    C->rev = malloc(C->MAXN * sizeof(char));
    C->upperBoundI = malloc(C->MAXN * sizeof(int));
    C->pos = malloc(C->MAXN * sizeof(*C->pos));
    C->bestbestscore = malloc(C->MAXN * sizeof(int));
    C->savebestscore = malloc(C->MAXN * sizeof(int));
    C->best_left_offset = malloc(C->MAXN * sizeof(int));
    C->save_left_offset = malloc(C->MAXN * sizeof(int));
    C->best_right_offset = malloc(C->MAXN * sizeof(int));
    C->save_right_offset = malloc(C->MAXN * sizeof(int));
    if (!C->rev || !C->upperBoundI || !C->pos || !C->bestbestscore || !C->savebestscore ||
        !C->best_left_offset || !C->save_left_offset ||
        !C->best_right_offset || !C->save_right_offset) return -1;

    /* score[2][MAXN][2*MAXOFFSET+1] */
    C->score = malloc(2 * sizeof(*C->score));
    if (!C->score) return -1;
    for (x = 0; x < 2; x++) {
        C->score[x] = malloc(C->MAXN * sizeof(*C->score[x]));
        if (!C->score[x]) return -1;
        for (n = 0; n < C->MAXN; n++) {
            C->score[x][n] = malloc((2 * C->MAXOFFSET + 1) * sizeof(*C->score[x][n]));
            if (!C->score[x][n]) return -1;
        }
    }

    /* score_of_besty[MAXN][2*MAXOFFSET+1] */
    C->score_of_besty = malloc(C->MAXN * sizeof(*C->score_of_besty));
    if (!C->score_of_besty) return -1;
    for (n = 0; n < C->MAXN; n++) {
        C->score_of_besty[n] = malloc((2 * C->MAXOFFSET + 1) * sizeof(*C->score_of_besty[n]));
        if (!C->score_of_besty[n]) return -1;
    }

    /* maskscore[2][2*MAXOFFSET+1] */
    C->maskscore = malloc(2 * sizeof(*C->maskscore));
    if (!C->maskscore) return -1;
    for (x = 0; x < 2; x++) {
        C->maskscore[x] = malloc((2 * C->MAXOFFSET + 1) * sizeof(*C->maskscore[x]));
        if (!C->maskscore[x]) return -1;
    }

    return 0;
}

static void free_workspace(DiscoverContext *C)
{
    int x, n;

    free(C->master); C->master = NULL;
    if (C->masters) {
        for (int r = 0; r < C->MAXR; r++) {
            if (C->masters_allocated && C->masters_allocated[r])
                free(C->masters[r]);
        }
        free(C->masters); C->masters = NULL;
    }
    free(C->masters_allocated); C->masters_allocated = NULL;
    free(C->masterstart); C->masterstart = NULL;
    free(C->masterend); C->masterend = NULL;

    free(C->rev); C->rev = NULL;
    free(C->upperBoundI); C->upperBoundI = NULL;
    free(C->pos); C->pos = NULL;
    free(C->bestbestscore); C->bestbestscore = NULL;
    free(C->savebestscore); C->savebestscore = NULL;
    free(C->best_left_offset); C->best_left_offset = NULL;
    free(C->save_left_offset); C->save_left_offset = NULL;
    free(C->best_right_offset); C->best_right_offset = NULL;
    free(C->save_right_offset); C->save_right_offset = NULL;

    if (C->score) {
        for (x = 0; x < 2; x++) {
            if (C->score[x]) {
                for (n = 0; n < C->MAXN; n++)
                    free(C->score[x][n]);
                free(C->score[x]);
            }
        }
        free(C->score); C->score = NULL;
    }

    if (C->score_of_besty) {
        for (n = 0; n < C->MAXN; n++)
            free(C->score_of_besty[n]);
        free(C->score_of_besty); C->score_of_besty = NULL;
    }

    if (C->maskscore) {
        for (x = 0; x < 2; x++)
            free(C->maskscore[x]);
        free(C->maskscore); C->maskscore = NULL;
    }

    free(C->removed); C->removed = NULL;
    free(C->disc_boundaries); C->disc_boundaries = NULL;
    free(C->sequence_owned); C->sequence_owned = NULL;
    C->sequence = NULL;
}

static void free_headptr(struct llist **headptr, gpos_t hash_size)
{
    int h;
    struct llist *tmp, *nexttmp;
    struct posllist *postmp, *nextpostmp;

    if (!headptr) return;

    for (h = 0; h < hash_size; h++) {
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

static void collect_instances_from_extension(DiscoverContext *C, CandidateFamily *fam)
{
    int cons_len = C->masterend[C->R] - C->masterstart[C->R] + 1;
    int n, i;

    fam->instances = malloc(C->N * sizeof(Instance));
    if (!fam->instances) {
        fprintf(stderr, "discover: out of memory for instances\n");
        fam->num_instances = 0;
        fam->cap_instances = 0;
        return;
    }
    fam->cap_instances = C->N;
    fam->num_instances = 0;

    for (n = 0; n < C->N; n++) {
        gpos_t genome_start, genome_end;

        if (C->rev[n] == 0) {
            /* Forward strand */
            genome_start = C->pos[n] + C->bestw - C->L;
            genome_end = C->pos[n] + C->besty - C->L + 1;
        } else {
            /* Reverse strand: extension reads genome in opposite direction */
            genome_start = C->pos[n] - (C->besty - C->L - C->l) - 1;
            genome_end = C->pos[n] - (C->bestw - C->L - C->l);
        }

        /* Bounds check — only keep instances from forward copy */
        if (genome_start < 0 || genome_end > C->orig_length)
            continue;

        glen_t aligned_length = genome_end - genome_start;
        if (aligned_length <= 0 || aligned_length > 100000)
            continue;

        /* Count mismatches by direct comparison */
        int edits = 0;
        int compare_len = (int)aligned_length < cons_len ? (int)aligned_length : cons_len;
        for (i = 0; i < compare_len; i++) {
            char seq_base;
            if (C->rev[n] == 0) {
                seq_base = C->sequence[genome_start + i];
            } else {
                seq_base = compl_base(C->sequence[genome_end - 1 - i]);
            }
            char cons_base = C->masters[C->R][C->masterstart[C->R] + i];
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
        inst->score = C->bestbestscore[n];
        inst->strand = C->rev[n] ? -1 : 1;
        inst->seq_index = find_seq_index(C, genome_start);
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

static char *build_doubled_genome(const char *orig_seq, glen_t orig_len, glen_t *new_len)
{
    glen_t doubled_len = 2 * orig_len + PADLENGTH;
    char *doubled = malloc((size_t)doubled_len * sizeof(char));
    if (!doubled) return NULL;

    /* Copy forward strand */
    memcpy(doubled, orig_seq, (size_t)orig_len);

    /* PADLENGTH spacer (DNA_N = 99) */
    memset(doubled + orig_len, 99, PADLENGTH);

    /* Reverse complement */
    for (glen_t x = 0; x < orig_len; x++) {
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
    params->GOODLENGTH     = 30;       /* plan: 50->30, MDL decides */
    params->MINTHRESH      = 2;        /* plan: 3->2, let MDL decide */
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

struct freq_entry {
    int freq;
    gpos_t occ;
};

static int cmp_freq_desc(const void *a, const void *b)
{
    const struct freq_entry *fa = (const struct freq_entry *)a;
    const struct freq_entry *fb = (const struct freq_entry *)b;
    if (fa->freq > fb->freq) return -1;
    if (fa->freq < fb->freq) return  1;
    return 0;
}

static void dump_freq_table(DiscoverContext *C, struct llist **hp, const char *path)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "discover: could not open %s for writing\n", path);
        return;
    }

    /* Count entries */
    int count = 0;
    for (gpos_t h = 0; h < C->hash_size; h++) {
        struct llist *tmp = hp[h];
        while (tmp) { count++; tmp = tmp->next; }
    }

    /* Collect (freq, occ) pairs */
    struct freq_entry *entries = malloc((size_t)count * sizeof(*entries));
    if (!entries) {
        fprintf(stderr, "discover: out of memory for freq dump\n");
        fclose(fp);
        return;
    }

    int idx = 0;
    for (gpos_t h = 0; h < C->hash_size; h++) {
        struct llist *tmp = hp[h];
        while (tmp) {
            entries[idx].freq = tmp->freq;
            entries[idx].occ = tmp->lastocc;
            idx++;
            tmp = tmp->next;
        }
    }

    /* Sort by frequency descending */
    qsort(entries, (size_t)count, sizeof(*entries), cmp_freq_desc);

    /* Write */
    char lmer_buf[256];
    for (int i = 0; i < count; i++) {
        gpos_t occ = entries[i].occ;
        for (int x = 0; x < C->l; x++)
            lmer_buf[x] = C->sequence[occ + x];
        reverse_if_necessary(lmer_buf, C->l);
        for (int x = 0; x < C->l; x++)
            fputc(num_to_char(lmer_buf[x]), fp);
        fprintf(fp, "\t%d\t%" PRId64 "\n", entries[i].freq, entries[i].occ);
    }

    free(entries);
    fclose(fp);
    fprintf(stderr, "discover: wrote %d l-mers to %s\n", count, path);
}

/* ================================================================
 * Public API: discover_families
 * ================================================================ */

CandidateList *discover_families(const Genome *genome,
                                 const DiscoverParams *params,
                                 int num_threads)
{
    struct llist **headptr;
    struct llist *besttmp;
    int x;

    /* ---- Allocate context ---- */
    DiscoverContext *C = calloc(1, sizeof(DiscoverContext));
    if (!C) return NULL;

    /* ---- Set parameters from DiscoverParams ---- */
    C->l              = params->l;
    C->L              = params->L;
    C->MAXOFFSET      = params->MAXOFFSET;
    C->MAXN           = params->MAXN;
    C->MAXR           = params->MAXR;
    C->MATCH          = params->MATCH;
    C->MISMATCH       = params->MISMATCH;
    C->GAP            = params->GAP;
    C->CAPPENALTY     = params->CAPPENALTY;
    C->MINIMPROVEMENT = params->MINIMPROVEMENT;
    C->WHEN_TO_STOP   = params->WHEN_TO_STOP;
    C->MAXENTROPY     = params->MAXENTROPY;
    C->GOODLENGTH     = params->GOODLENGTH;
    C->MINTHRESH      = params->MINTHRESH;
    C->TANDEMDIST     = params->TANDEMDIST;
    C->VERBOSE        = params->VERBOSE;

    /* ---- Setup genome data ---- */
    C->orig_length = genome->length;
    C->disc_num_sequences = genome->num_sequences;

    /* Auto-compute l if not specified (use original length) */
    if (C->l <= 0)
        C->l = (int)ceil(1.0 + log((double)C->orig_length) / log(4.0));

    /* Dynamic hash size: at least HASH_SIZE_MIN, but scale with genome/l.
     * 4 * genome_len / l gives ~4 entries per bucket on average for a
     * fully loaded table.  Round up to next odd number (not necessarily
     * prime, but good enough for this chained hash). */
    {
        gpos_t dynamic_size = 4 * C->orig_length / C->l;
        if (dynamic_size % 2 == 0) dynamic_size++;   /* keep odd */
        C->hash_size = (dynamic_size > HASH_SIZE_MIN) ? dynamic_size : HASH_SIZE_MIN;
    }
    fprintf(stderr, "discover: hash_size=%" PRId64 " (genome=%" PRId64 " l=%d)\n",
            C->hash_size, C->orig_length, C->l);

    fprintf(stderr, "discover: l=%d, L=%d, MAXOFFSET=%d, MINTHRESH=%d, GOODLENGTH=%d\n",
            C->l, C->L, C->MAXOFFSET, C->MINTHRESH, C->GOODLENGTH);

    /* Build doubled genome: forward + PADLENGTH padding + RC */
    glen_t doubled_len;
    C->sequence_owned = build_doubled_genome(genome->sequence, C->orig_length, &doubled_len);
    if (!C->sequence_owned) {
        fprintf(stderr, "discover: out of memory for doubled genome\n");
        free(C);
        return NULL;
    }
    C->sequence = C->sequence_owned;
    C->length = doubled_len;

    /* Copy boundaries to gpos_t array */
    C->disc_boundaries = malloc((C->disc_num_sequences + 1) * sizeof(gpos_t));
    if (!C->disc_boundaries) {
        fprintf(stderr, "discover: out of memory for boundaries\n");
        free(C->sequence_owned);
        free(C);
        return NULL;
    }
    for (int i = 0; i <= C->disc_num_sequences; i++)
        C->disc_boundaries[i] = genome->boundaries[i];

    /* Allocate masking array (covers doubled genome) */
    C->removed = calloc((size_t)C->length, sizeof(char));
    if (!C->removed) {
        fprintf(stderr, "discover: out of memory for masking array\n");
        free(C->disc_boundaries);
        free(C->sequence_owned);
        free(C);
        return NULL;
    }

    /* ---- Allocate workspace ---- */
    if (allocate_workspace(C) < 0) {
        fprintf(stderr, "discover: out of memory for workspace\n");
        free(C->removed);
        free(C->disc_boundaries);
        free(C->sequence_owned);
        free(C);
        return NULL;
    }

    /* ---- Build hash table ---- */
    headptr = calloc((size_t)C->hash_size, sizeof(*headptr));
    if (!headptr) {
        fprintf(stderr, "discover: out of memory for headptr\n");
        free_workspace(C);
        free(C);
        return NULL;
    }

    if (params->freq_file) {
        fprintf(stderr, "discover: reading frequency table from %s\n", params->freq_file);
        build_headptr_from_freq(C, headptr, params->freq_file);
    } else if (num_threads > 1) {
        fprintf(stderr, "discover: counting l-mers in parallel (%d threads)...\n",
                num_threads);
        build_headptr_parallel(C, headptr, num_threads);
    } else {
        fprintf(stderr, "discover: counting l-mers internally...\n");
        build_headptr_internal(C, headptr);
    }

    /* Dump frequency table if requested (before trim, same as build_lmer_table) */
    if (params->freq_output) {
        dump_freq_table(C, headptr, params->freq_output);
    }

    trim_headptr(C, headptr);
    build_all_pos(C, headptr);

    /* ---- Allocate result ---- */
    CandidateList *result = malloc(sizeof(CandidateList));
    if (!result) {
        fprintf(stderr, "discover: out of memory for result\n");
        free_headptr(headptr, C->hash_size);
        free_workspace(C);
        free(C);
        return NULL;
    }
    result->cap_families = 1024;
    result->families = malloc(result->cap_families * sizeof(CandidateFamily));
    if (!result->families) {
        fprintf(stderr, "discover: out of memory for families\n");
        free(result);
        free_headptr(headptr, C->hash_size);
        free_workspace(C);
        free(C);
        return NULL;
    }
    result->num_families = 0;

    /* ---- Main discovery loop ---- */
    C->R = 0;
    C->prevbestfreq = 1000000000;
    C->prevbesthash = 0;

    while (1) {
        besttmp = find_besttmp(C, headptr);

        if (besttmp == NULL || besttmp->freq < C->MINTHRESH) {
            if (C->VERBOSE)
                fprintf(stderr, "discover: stopped at R=%d (no more frequent l-mers)\n", C->R);
            break;
        }

        /* Initialize seed in master */
        for (x = 0; x < C->l; x++)
            C->master[C->L + x] = C->sequence[besttmp->lastocc + x];

        /* Build occurrence list */
        build_pos(C, besttmp);

        if (C->N < C->MINTHRESH) {
            if (C->VERBOSE >= 2)
                fprintf(stderr, "discover: R=%d N=%d < MINTHRESH=%d, updating freq\n",
                        C->R, C->N, C->MINTHRESH);
            besttmp->freq = C->N;
            continue;
        }

        /* Initialize offset tracking */
        for (int n = 0; n < C->N; n++) {
            C->best_left_offset[n] = -1;
            C->save_left_offset[n] = -1;
            C->best_right_offset[n] = -1;
            C->save_right_offset[n] = -1;
        }

        /* Allocate masters[R] if needed */
        if (C->masters_allocated[C->R] == 0) {
            C->masters[C->R] = malloc((2 * C->L + C->l) * sizeof(char));
            if (!C->masters[C->R]) {
                fprintf(stderr, "discover: out of memory for masters[%d]\n", C->R);
                break;
            }
            C->masters_allocated[C->R] = 1;
        }

        /* Extend right and left */
        extend_right(C);
        extend_left(C);

        int cons_len = C->masterend[C->R] - C->masterstart[C->R] + 1;

        if (cons_len >= C->GOODLENGTH) {
            /* Good family: save and collect instances */
            if (C->VERBOSE) {
                fprintf(stderr, "discover: R=%d N=%d length=%d\n", C->R, C->N, cons_len);
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
            fam->id = (mdl_uid_t)result->num_families;

            /* Copy consensus */
            fam->consensus = malloc(cons_len * sizeof(char));
            if (!fam->consensus) {
                fprintf(stderr, "discover: out of memory for consensus\n");
                break;
            }
            for (x = 0; x < cons_len; x++)
                fam->consensus[x] = C->masters[C->R][C->masterstart[C->R] + x];
            fam->consensus_length = cons_len;

            /* Set default fields */
            fam->topology = TOPO_LINEAR;
            fam->component_id = 0;
            fam->estimated_copies = C->N;

            /* Collect instances from extension occurrences */
            collect_instances_from_extension(C, fam);

            result->num_families++;

            C->R++;
            if (C->R == C->MAXR) break;

            /* Mask this family */
            mask_headptr(C, headptr);
        } else {
            /* Short family: mask but don't record */
            C->R++;
            if (C->R == C->MAXR) break;
            mask_headptr(C, headptr);
            C->R--;
        }
    }

    fprintf(stderr, "discover: found %d families (R=%d total including short)\n",
            result->num_families, C->R);

    /* ---- Cleanup ---- */
    free_headptr(headptr, C->hash_size);
    free_workspace(C);
    free(C);

    return result;
}
