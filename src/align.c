#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <pthread.h>

#include "align.h"

/* ================================================================
 * Internal constants
 * ================================================================ */

#define MAX_SEED_HITS    50000   /* cap per-family seed hits */
#define MAX_CONS_KMERS   10000   /* cap consensus k-mer set entries */
#define BAND_WIDTH       (2 * ALIGN_MAXOFFSET + 1)
#define ALIGN_MAX_EXTENSION 10000  /* max bases to extend per direction */

/* ================================================================
 * Consensus k-mer hash set (small, for one family at a time)
 * ================================================================ */

typedef struct ckmer_entry {
    uint64_t canon;            /* canonical k-mer */
    uint64_t fwd;              /* forward (non-canonical) packed k-mer */
    int      cons_pos;         /* position in consensus */
    struct ckmer_entry *next;
} CKmerEntry;

typedef struct {
    CKmerEntry **buckets;
    size_t       size;         /* prime table size */
    int          count;
} CKmerSet;

static CKmerSet *ckmer_set_create(int expected)
{
    CKmerSet *cs = calloc(1, sizeof(CKmerSet));
    if (!cs) return NULL;
    /* Use ~2x expected entries for low collision */
    size_t sz = (size_t)(expected * 2 + 1);
    if (sz < 127) sz = 127;
    cs->size = sz;
    cs->buckets = calloc(sz, sizeof(CKmerEntry *));
    if (!cs->buckets) { free(cs); return NULL; }
    return cs;
}

static void ckmer_set_free(CKmerSet *cs)
{
    if (!cs) return;
    for (size_t i = 0; i < cs->size; i++) {
        CKmerEntry *e = cs->buckets[i];
        while (e) {
            CKmerEntry *next = e->next;
            free(e);
            e = next;
        }
    }
    free(cs->buckets);
    free(cs);
}

static void ckmer_set_insert(CKmerSet *cs, uint64_t canon, uint64_t fwd,
                              int cons_pos)
{
    size_t h = (size_t)(canon % cs->size);
    /* Check for duplicate canonical k-mer — keep first occurrence */
    for (CKmerEntry *e = cs->buckets[h]; e; e = e->next) {
        if (e->canon == canon) return;
    }
    CKmerEntry *e = malloc(sizeof(CKmerEntry));
    if (!e) return;
    e->canon = canon;
    e->fwd = fwd;
    e->cons_pos = cons_pos;
    e->next = cs->buckets[h];
    cs->buckets[h] = e;
    cs->count++;
}

static CKmerEntry *ckmer_set_lookup(const CKmerSet *cs, uint64_t canon)
{
    size_t h = (size_t)(canon % cs->size);
    for (CKmerEntry *e = cs->buckets[h]; e; e = e->next) {
        if (e->canon == canon) return e;
    }
    return NULL;
}

/* ================================================================
 * Step 1: Multi-k-mer seeding via position index lookup
 * ================================================================ */

/*
 * For each unique consensus k-mer, look up its positions from the
 * KmerTable position index.  O(sum of freqs) instead of O(genome_length).
 */
static int seed_genome_scan(const CandidateFamily *fam, const KmerTable *kt,
                            int k, SeedHit *hits, int max_hits)
{
    int cons_len = fam->consensus_length;
    if (cons_len < k) return 0;

    int n_cons_kmers = cons_len - k + 1;
    if (n_cons_kmers > MAX_CONS_KMERS) n_cons_kmers = MAX_CONS_KMERS;

    /* Use CKmerSet to deduplicate canonical k-mers in the consensus
     * (same k-mer appearing at multiple positions → only process once) */
    CKmerSet *processed = ckmer_set_create(n_cons_kmers);
    if (!processed) return 0;

    int n_hits = 0;

    for (int i = 0; i < n_cons_kmers && n_hits < max_hits; i++) {
        uint64_t fwd = kmer_pack(fam->consensus + i, k);
        if (fwd == UINT64_MAX) continue;
        uint64_t canon = kmer_canonical(fwd, k);

        /* Skip if we've already processed this canonical k-mer */
        if (ckmer_set_lookup(processed, canon)) continue;
        ckmer_set_insert(processed, canon, fwd, i);

        KmerEntry *entry = kmer_lookup(kt, canon);
        if (!entry || entry->num_positions == 0) continue;

        int cons_is_rc = (fwd != canon);

        for (int pi = 0; pi < entry->num_positions && n_hits < max_hits; pi++) {
            gpos_t sp = entry->positions[pi];
            gpos_t pos = (sp < 0) ? -sp : sp;
            int genome_is_rc = (sp < 0);

            int8_t strand = (cons_is_rc == genome_is_rc) ? +1 : -1;

            hits[n_hits].genome_pos = pos;
            hits[n_hits].cons_pos = i;
            hits[n_hits].strand = strand;
            n_hits++;
        }
    }

    ckmer_set_free(processed);
    return n_hits;
}

/* ================================================================
 * Step 2: Seed clustering (cluster_seed_hits)
 * ================================================================ */

/* Sort seed hits by (strand, genome_pos) */
static int cmp_seed_hits(const void *a, const void *b)
{
    const SeedHit *ha = (const SeedHit *)a;
    const SeedHit *hb = (const SeedHit *)b;
    if (ha->strand != hb->strand) return ha->strand - hb->strand;
    if (ha->genome_pos < hb->genome_pos) return -1;
    if (ha->genome_pos > hb->genome_pos) return 1;
    return 0;
}

/*
 * Cluster nearby same-strand hits into candidate instances.
 * For each cluster, pick the hit closest to consensus center as anchor.
 * Returns number of clusters. Fills anchors[] with selected anchor hits.
 */
static int cluster_seed_hits(SeedHit *hits, int n_hits, int cons_len,
                             SeedHit *anchors, int max_anchors)
{
    if (n_hits == 0) return 0;

    qsort(hits, (size_t)n_hits, sizeof(SeedHit), cmp_seed_hits);

    glen_t merge_dist = (glen_t)(cons_len * 1.5);
    int cons_mid = cons_len / 2;
    int n_anchors = 0;

    int cluster_start = 0;
    for (int i = 1; i <= n_hits && n_anchors < max_anchors; i++) {
        int end_cluster = 0;
        if (i == n_hits) {
            end_cluster = 1;
        } else if (hits[i].strand != hits[i - 1].strand) {
            end_cluster = 1;
        } else if (hits[i].genome_pos - hits[i - 1].genome_pos > merge_dist) {
            end_cluster = 1;
        }

        if (end_cluster) {
            /* Pick anchor closest to consensus center */
            int best = cluster_start;
            int best_dist = abs(hits[cluster_start].cons_pos - cons_mid);
            for (int j = cluster_start + 1; j < i; j++) {
                int d = abs(hits[j].cons_pos - cons_mid);
                if (d < best_dist) {
                    best_dist = d;
                    best = j;
                }
            }
            anchors[n_anchors++] = hits[best];
            cluster_start = i;
        }
    }

    return n_anchors;
}

/* ================================================================
 * Step 3: Banded alignment (align_banded)
 * ================================================================ */

/*
 * Access a genome base, with strand handling.
 * For forward strand: genome[pos]
 * For reverse strand: complement of genome[pos]
 */
static inline char genome_base(const Genome *genome, gpos_t pos, int8_t strand)
{
    if (pos < 0 || pos >= genome->length) return DNA_N;
    char b = genome->sequence[pos];
    if (strand < 0) return dna_complement(b);
    return b;
}

/*
 * Compute the genome position for consensus position i, given an anchor.
 * For forward strand: anchor_genome + (i - anchor_cons)
 * For reverse strand: anchor_genome - (i - anchor_cons)
 *   (reverse strand reads genome backwards)
 */
static inline gpos_t genome_pos_for_cons(gpos_t anchor_genome, int anchor_cons,
                                         int cons_i, int8_t strand)
{
    if (strand > 0)
        return anchor_genome + (cons_i - anchor_cons);
    else
        return anchor_genome - (cons_i - anchor_cons);
}

/*
 * Banded DP alignment extending from a seed anchor.
 * Extends right from anchor, then left, using RepeatScout's scoring.
 *
 * Returns an AlignedInstance with boundaries and score.
 * Returns .score = INT_MIN on failure.
 */
static AlignedInstance align_banded(const CandidateFamily *fam,
                                   const Genome *genome,
                                   const SeedHit *anchor)
{
    AlignedInstance result;
    memset(&result, 0, sizeof(result));
    result.score = INT_MIN;
    result.strand = anchor->strand;

    int cons_len = fam->consensus_length;
    int anchor_cons = anchor->cons_pos;
    gpos_t anchor_genome = anchor->genome_pos;

    /* DP arrays: two rows, each of width BAND_WIDTH.
     * offset index: o ∈ [0, BAND_WIDTH-1], actual offset = o - MAXOFFSET */
    int dp[2][BAND_WIDTH];
    int prev_row, curr_row;

    /* --- Extend RIGHT from anchor --- */
    int best_right_score = 0;
    int best_right_i = anchor_cons;
    int right_end_cons = anchor_cons; /* exclusive */

    /* Initialize DP at anchor position */
    prev_row = 0;
    for (int o = 0; o < BAND_WIDTH; o++)
        dp[prev_row][o] = ALIGN_CAPPENALTY;
    dp[prev_row][ALIGN_MAXOFFSET] = 0; /* offset 0 = exact alignment */

    /* Score the anchor position itself */
    {
        gpos_t gp = genome_pos_for_cons(anchor_genome, anchor_cons,
                                        anchor_cons, anchor->strand);
        char gb = genome_base(genome, gp, anchor->strand);
        char cb = fam->consensus[anchor_cons];
        if (gb != DNA_N && cb != DNA_N) {
            int s = (gb == cb) ? ALIGN_MATCH : ALIGN_MISMATCH;
            dp[prev_row][ALIGN_MAXOFFSET] = s;
            best_right_score = s;
        }
        right_end_cons = anchor_cons + 1;
    }

    int stall = 0;
    for (int i = anchor_cons + 1; i < cons_len; i++) {
        curr_row = 1 - prev_row;

        for (int o = 0; o < BAND_WIDTH; o++)
            dp[curr_row][o] = ALIGN_CAPPENALTY;

        for (int o = 0; o < BAND_WIDTH; o++) {
            int actual_offset = o - ALIGN_MAXOFFSET;
            gpos_t gp = genome_pos_for_cons(anchor_genome, anchor_cons,
                                            i, anchor->strand);
            gp += actual_offset;

            char gb = genome_base(genome, gp, anchor->strand);
            char cb = fam->consensus[i];

            if (gb == DNA_N || cb == DNA_N) continue;

            int diag = (gb == cb) ? ALIGN_MATCH : ALIGN_MISMATCH;

            /* Match/mismatch from same offset */
            int val = dp[prev_row][o] + diag;

            /* Gap in genome: offset shifts by +1 (previous o+1) */
            if (o + 1 < BAND_WIDTH && dp[prev_row][o + 1] + ALIGN_GAP > val)
                val = dp[prev_row][o + 1] + ALIGN_GAP;

            /* Gap in consensus: offset shifts by -1 (previous o-1) */
            if (o > 0 && dp[prev_row][o - 1] + ALIGN_GAP > val)
                val = dp[prev_row][o - 1] + ALIGN_GAP;

            /* Cap penalty floor */
            if (val < best_right_score + ALIGN_CAPPENALTY)
                val = best_right_score + ALIGN_CAPPENALTY;

            dp[curr_row][o] = val;

            if (val > best_right_score) {
                best_right_score = val;
                best_right_i = i;
                stall = 0;
            }
        }

        stall++;
        if (stall >= ALIGN_WHEN_TO_STOP) break;

        prev_row = curr_row;
    }
    right_end_cons = best_right_i + 1;

    /* --- Extend LEFT from anchor --- */
    int best_left_score = 0;
    int best_left_i = anchor_cons;
    int left_start_cons = anchor_cons;

    /* Re-initialize DP */
    prev_row = 0;
    for (int o = 0; o < BAND_WIDTH; o++)
        dp[prev_row][o] = ALIGN_CAPPENALTY;
    dp[prev_row][ALIGN_MAXOFFSET] = 0;

    stall = 0;
    for (int i = anchor_cons - 1; i >= 0; i--) {
        curr_row = 1 - prev_row;

        for (int o = 0; o < BAND_WIDTH; o++)
            dp[curr_row][o] = ALIGN_CAPPENALTY;

        for (int o = 0; o < BAND_WIDTH; o++) {
            int actual_offset = o - ALIGN_MAXOFFSET;
            gpos_t gp = genome_pos_for_cons(anchor_genome, anchor_cons,
                                            i, anchor->strand);
            gp += actual_offset;

            char gb = genome_base(genome, gp, anchor->strand);
            char cb = fam->consensus[i];

            if (gb == DNA_N || cb == DNA_N) continue;

            int diag = (gb == cb) ? ALIGN_MATCH : ALIGN_MISMATCH;

            int val = dp[prev_row][o] + diag;

            if (o + 1 < BAND_WIDTH && dp[prev_row][o + 1] + ALIGN_GAP > val)
                val = dp[prev_row][o + 1] + ALIGN_GAP;

            if (o > 0 && dp[prev_row][o - 1] + ALIGN_GAP > val)
                val = dp[prev_row][o - 1] + ALIGN_GAP;

            if (val < best_left_score + ALIGN_CAPPENALTY)
                val = best_left_score + ALIGN_CAPPENALTY;

            dp[curr_row][o] = val;

            if (val > best_left_score) {
                best_left_score = val;
                best_left_i = i;
                stall = 0;
            }
        }

        stall++;
        if (stall >= ALIGN_WHEN_TO_STOP) break;

        prev_row = curr_row;
    }
    left_start_cons = best_left_i;

    /* Compute aligned region */
    int aligned_cons_len = right_end_cons - left_start_cons;
    if (aligned_cons_len < 10) return result; /* too short */

    int total_score = best_right_score + best_left_score;
    if (total_score <= 0) return result;

    /* Compute genome coordinates */
    gpos_t gstart = genome_pos_for_cons(anchor_genome, anchor_cons,
                                        left_start_cons, anchor->strand);
    gpos_t gend = genome_pos_for_cons(anchor_genome, anchor_cons,
                                      right_end_cons - 1, anchor->strand);
    if (gstart > gend) {
        gpos_t tmp = gstart;
        gstart = gend;
        gend = tmp;
    }
    gend++; /* exclusive */

    /* Count actual edits in the aligned region */
    int edits = 0;
    int compared = 0;
    for (int i = left_start_cons; i < right_end_cons; i++) {
        gpos_t gp = genome_pos_for_cons(anchor_genome, anchor_cons,
                                        i, anchor->strand);
        char gb = genome_base(genome, gp, anchor->strand);
        char cb = fam->consensus[i];
        if (gb == DNA_N || cb == DNA_N) continue;
        compared++;
        if (gb != cb) edits++;
    }

    float div = (compared > 0) ? (float)edits / (float)compared : 1.0f;
    if (div > ALIGN_MAX_DIVERGENCE) return result;

    /* Fill result */
    result.genome_start = gstart;
    result.genome_end = gend;
    result.cons_start = left_start_cons;
    result.cons_end = right_end_cons;
    result.num_edits = edits;
    result.divergence = div;
    result.score = total_score;
    result.strand = anchor->strand;

    return result;
}

/* ================================================================
 * Step 4: align_collect_instances — orchestrator
 * ================================================================ */

int align_collect_instances(CandidateFamily *fam, const Genome *genome,
                            const KmerTable *kt, int k)
{
    /* Free old instances */
    free(fam->instances);
    fam->instances = NULL;
    fam->num_instances = 0;
    fam->cap_instances = 0;

    if (fam->consensus_length < k) return 0;

    /* Step 1: Multi-k-mer genome scan */
    SeedHit *hits = malloc((size_t)MAX_SEED_HITS * sizeof(SeedHit));
    if (!hits) return 0;

    int n_hits = seed_genome_scan(fam, kt, k, hits, MAX_SEED_HITS);
    if (n_hits == 0) {
        free(hits);
        return 0;
    }

    /* Step 2: Cluster hits into candidate instances */
    int max_anchors = DEFAULT_MAXN;
    SeedHit *anchors = malloc((size_t)max_anchors * sizeof(SeedHit));
    if (!anchors) {
        free(hits);
        return 0;
    }

    int n_anchors = cluster_seed_hits(hits, n_hits, fam->consensus_length,
                                      anchors, max_anchors);
    free(hits);

    if (n_anchors == 0) {
        free(anchors);
        return 0;
    }

    /* Pre-allocate instance array */
    int cap = (n_anchors < 64) ? 64 : n_anchors;
    if (cap > DEFAULT_MAXN) cap = DEFAULT_MAXN;
    fam->instances = malloc((size_t)cap * sizeof(Instance));
    if (!fam->instances) {
        free(anchors);
        fam->cap_instances = 0;
        return 0;
    }
    fam->cap_instances = cap;

    /* Step 3: Banded alignment for each anchor */
    for (int i = 0; i < n_anchors && fam->num_instances < DEFAULT_MAXN; i++) {
        AlignedInstance ai = align_banded(fam, genome, &anchors[i]);
        if (ai.score == INT_MIN) continue;

        /* Boundary check */
        glen_t alen = ai.genome_end - ai.genome_start;
        if (ai.genome_start < PADLENGTH) continue;
        if (ai.genome_end > genome->length) continue;
        if (!genome_check_boundary(genome, ai.genome_start, alen))
            continue;

        /* Check for duplicate: skip if overlapping an existing instance */
        int is_dup = 0;
        for (int j = 0; j < fam->num_instances; j++) {
            Instance *existing = &fam->instances[j];
            gpos_t ex_end = existing->position + existing->aligned_length;
            /* Overlap if they share > 50% of the shorter one */
            gpos_t ov_start = (ai.genome_start > existing->position)
                              ? ai.genome_start : existing->position;
            gpos_t ov_end = (ai.genome_end < ex_end)
                            ? ai.genome_end : ex_end;
            if (ov_end > ov_start) {
                glen_t overlap = ov_end - ov_start;
                glen_t shorter = (alen < existing->aligned_length)
                                 ? alen : existing->aligned_length;
                if (overlap > shorter / 2) {
                    is_dup = 1;
                    break;
                }
            }
        }
        if (is_dup) continue;

        /* Grow array if needed */
        if (fam->num_instances >= fam->cap_instances) {
            int new_cap = fam->cap_instances * 2;
            if (new_cap > DEFAULT_MAXN) new_cap = DEFAULT_MAXN;
            Instance *tmp = realloc(fam->instances,
                                    (size_t)new_cap * sizeof(Instance));
            if (!tmp) break;
            fam->instances = tmp;
            fam->cap_instances = new_cap;
        }

        Instance *inst = &fam->instances[fam->num_instances];
        inst->position = ai.genome_start;
        inst->aligned_length = ai.genome_end - ai.genome_start;
        inst->cons_start = ai.cons_start;
        inst->cons_end = ai.cons_end;
        inst->num_edits = ai.num_edits;
        inst->divergence = ai.divergence;
        inst->score = ai.score;
        inst->strand = ai.strand;
        inst->seq_index = genome_get_seq_index(genome, ai.genome_start);
        fam->num_instances++;
    }

    free(anchors);
    return fam->num_instances;
}

/* ================================================================
 * Step 5: Consensus rebuild via majority voting
 * ================================================================ */

int align_rebuild_consensus(CandidateFamily *fam, const Genome *genome)
{
    if (fam->num_instances < 2) return 0;

    int cons_len = fam->consensus_length;
    int changed = 0;

    for (int p = 0; p < cons_len; p++) {
        int count[4] = {0, 0, 0, 0};

        for (int j = 0; j < fam->num_instances; j++) {
            Instance *inst = &fam->instances[j];

            /* Only count if this instance covers consensus position p */
            if (p < inst->cons_start || p >= inst->cons_end)
                continue;

            /* Weight by alignment score (higher quality → more influence) */
            int weight = inst->score;
            if (weight < 1) weight = 1;

            /* Map consensus position to genome position.
             * For forward strand: genome_pos = position + (p - cons_start)
             * For reverse strand: genome_pos = position + (cons_end - 1 - p) */
            gpos_t gp;
            if (inst->strand > 0) {
                gp = inst->position + (p - inst->cons_start);
            } else {
                gp = inst->position + (inst->cons_end - 1 - p);
            }

            if (gp < 0 || gp >= genome->length) continue;

            char base = genome->sequence[gp];
            if (inst->strand < 0)
                base = dna_complement(base);

            if (base >= 0 && base <= 3)
                count[(int)base] += weight;
        }

        /* Find majority base */
        int total = count[0] + count[1] + count[2] + count[3];
        if (total == 0) continue;

        int best_base = 0;
        for (int b = 1; b < 4; b++) {
            if (count[b] > count[best_base])
                best_base = b;
        }

        if (fam->consensus[p] != (char)best_base) {
            fam->consensus[p] = (char)best_base;
            changed++;
        }
    }

    return changed;
}

/* ================================================================
 * Step 5b: Consensus extension via flanking context
 * ================================================================ */

/*
 * Extend consensus in one direction using majority voting over instance
 * flanking genome bases.  direction: +1 = right, -1 = left.
 * Returns number of bases extended.
 */
static int extend_direction(CandidateFamily *fam, const Genome *genome,
                            int direction, int max_ext)
{
    int extended = 0;
    int best_score = 0;
    int score = 0;
    int stall = 0;
    int best_extended = 0;

    if (max_ext > ALIGN_MAX_EXTENSION) max_ext = ALIGN_MAX_EXTENSION;

    char *ext_buf = malloc((size_t)max_ext);
    if (!ext_buf) return 0;

    for (int offset = 0; offset < max_ext; offset++) {
        /* Consensus position being probed (beyond current bounds) */
        int p;
        if (direction > 0)
            p = fam->consensus_length + offset;
        else
            p = -(offset + 1);

        int count[4] = {0, 0, 0, 0};

        for (int j = 0; j < fam->num_instances; j++) {
            Instance *inst = &fam->instances[j];

            /* Only use instances that reach near the edge being extended */
            if (direction > 0) {
                if (inst->cons_end < fam->consensus_length - EXTENSION_SLACK) continue;
            } else {
                if (inst->cons_start > EXTENSION_SLACK) continue;
            }

            /* Compute genome position using existing linear mapping */
            gpos_t gp;
            if (inst->strand > 0) {
                gp = inst->position + (p - inst->cons_start);
            } else {
                gp = inst->position + (inst->cons_end - 1 - p);
            }

            if (gp < 0 || gp >= genome->length) continue;

            char base = genome->sequence[gp];
            if (inst->strand < 0)
                base = dna_complement(base);

            if (base >= 0 && base <= 3)
                count[(int)base]++;
        }

        int total = count[0] + count[1] + count[2] + count[3];
        if (total < 2) break;  /* insufficient coverage */

        /* Dynamic minimum support: require more instances for longer extensions.
         * Prevents segdup over-extension: 2 identical copies score +1 at every
         * position (100% agreement), so WHEN_TO_STOP never triggers. This check
         * ensures that extensions > 1000bp need ≥ 3 supporting instances. */
        if (extended > 1000 && total < 3) break;
        if (extended > 5000 && total < 5) break;

        /* Find majority base */
        int best_base = 0;
        for (int b = 1; b < 4; b++) {
            if (count[b] > count[best_base])
                best_base = b;
        }

        /* Hard stop if no clear majority (< 50%) */
        if (count[best_base] * 2 < total) break;

        ext_buf[extended] = (char)best_base;
        extended++;

        /* Score tracking: +1 if strong support (>=60%), -1 otherwise */
        score += (count[best_base] * 10 >= total * 6) ? 1 : -1;

        if (score > best_score) {
            best_score = score;
            best_extended = extended;
            stall = 0;
        } else {
            stall++;
        }

        if (stall >= ALIGN_WHEN_TO_STOP) break;

        /* Cap penalty: hard stop if score drops too far below best */
        if (score < best_score + ALIGN_CAPPENALTY) break;
    }

    /* Trim to best scoring position */
    extended = best_extended;

    if (extended == 0) {
        free(ext_buf);
        return 0;
    }

    /* Reallocate consensus */
    int new_len = fam->consensus_length + extended;
    char *new_cons = realloc(fam->consensus, (size_t)new_len);
    if (!new_cons) {
        free(ext_buf);
        return 0;
    }

    if (direction > 0) {
        /* Append to right */
        memcpy(new_cons + fam->consensus_length, ext_buf, (size_t)extended);
    } else {
        /* Prepend to left: shift existing consensus right, insert reversed */
        memmove(new_cons + extended, new_cons, (size_t)fam->consensus_length);
        for (int i = 0; i < extended; i++)
            new_cons[i] = ext_buf[extended - 1 - i];

        /* Shift all instance coordinates to match new consensus origin */
        for (int j = 0; j < fam->num_instances; j++) {
            fam->instances[j].cons_start += extended;
            fam->instances[j].cons_end += extended;
        }
    }

    fam->consensus = new_cons;
    fam->consensus_length = new_len;

    free(ext_buf);
    return extended;
}

int align_extend_consensus(CandidateFamily *fam, const Genome *genome)
{
    if (fam->num_instances < 2) return 0;

    /* Extension cap based on instance count.
     * More instances = higher confidence = allow longer extension per iteration.
     * Few instances (2-3) get conservative cap to prevent segdup runaway.
     * Many instances (10+) get full extension — the WHEN_TO_STOP and CAPPENALTY
     * mechanisms naturally stop at the TE boundary.
     * This replaces the old consensus-length-based cap that limited growth of
     * initially short graph unitigs (median 222bp → stuck at short lengths). */
    int cap;
    if (fam->num_instances <= 3)
        cap = 500;
    else if (fam->num_instances <= 9)
        cap = 3000;
    else
        cap = ALIGN_MAX_EXTENSION;

    int ext_right = extend_direction(fam, genome, +1, cap);
    int ext_left  = extend_direction(fam, genome, -1, cap);

    return ext_right + ext_left;
}

/* ================================================================
 * Step 6: Iteration loop
 * ================================================================ */

int align_refine_family(CandidateFamily *fam, const Genome *genome,
                        const KmerTable *kt, int k, int max_iter)
{
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        align_collect_instances(fam, genome, kt, k);

        if (fam->num_instances < 2)
            break;

        int extended = align_extend_consensus(fam, genome);
        int changed = align_rebuild_consensus(fam, genome);

        /* Converged: no extension and consensus stable */
        if (extended == 0 && changed == 0)
            break;
        if (extended == 0 && changed < fam->consensus_length / 100 + 1)
            break;
    }

    return fam->num_instances;
}

/* ================================================================
 * Step 7: Batch refinement (parallel)
 * ================================================================ */

/* Shared state for worker threads */
typedef struct {
    CandidateList  *cl;
    const Genome   *genome;
    const KmerTable *kt;
    int             k;
    int             verbose;
    int             next_family;    /* atomic: next family index to process */
    int             done_count;     /* atomic: families completed so far */
    pthread_mutex_t progress_lock;  /* protects stderr progress output */
} RefineWorkerState;

static void *refine_worker(void *arg)
{
    RefineWorkerState *state = (RefineWorkerState *)arg;

    while (1) {
        int idx = __atomic_fetch_add(&state->next_family, 1, __ATOMIC_SEQ_CST);
        if (idx >= state->cl->num_families)
            break;

        CandidateFamily *fam = &state->cl->families[idx];
        align_refine_family(fam, state->genome, state->kt, state->k,
                            ALIGN_MAX_ITERATIONS);

        int done = __atomic_add_fetch(&state->done_count, 1, __ATOMIC_SEQ_CST);
        if (state->verbose && done % 100 == 0) {
            pthread_mutex_lock(&state->progress_lock);
            fprintf(stderr, "  Refined %d / %d families\r",
                    done, state->cl->num_families);
            pthread_mutex_unlock(&state->progress_lock);
        }
    }

    return NULL;
}

void align_refine_all(CandidateList *cl, const Genome *genome,
                      const KmerTable *kt, int k, int num_threads,
                      int verbose)
{
    if (!cl || cl->num_families == 0) return;

    if (num_threads < 1) num_threads = 1;

    /* Cap threads at family count */
    if (num_threads > cl->num_families)
        num_threads = cl->num_families;

    if (num_threads == 1) {
        /* Sequential path (no thread overhead) */
        int refined = 0, lost = 0;
        for (int i = 0; i < cl->num_families; i++) {
            CandidateFamily *fam = &cl->families[i];
            align_refine_family(fam, genome, kt, k, ALIGN_MAX_ITERATIONS);

            if (fam->num_instances >= 2) refined++;
            else lost++;

            if (verbose && (i + 1) % 100 == 0)
                fprintf(stderr, "  Refined %d / %d families\r",
                        i + 1, cl->num_families);
        }
        if (verbose)
            fprintf(stderr, "  Alignment refinement: %d families refined, "
                    "%d lost all instances\n", refined, lost);
        return;
    }

    /* Parallel path */
    if (verbose)
        fprintf(stderr, "  Using %d threads\n", num_threads);

    RefineWorkerState state;
    state.cl = cl;
    state.genome = genome;
    state.kt = kt;
    state.k = k;
    state.verbose = verbose;
    state.next_family = 0;
    state.done_count = 0;
    pthread_mutex_init(&state.progress_lock, NULL);

    pthread_t *threads = malloc((size_t)num_threads * sizeof(pthread_t));
    if (!threads) {
        fprintf(stderr, "WARNING: thread allocation failed, falling back to single-threaded\n");
        align_refine_all(cl, genome, kt, k, 1, verbose);
        return;
    }

    /* Launch workers */
    for (int t = 0; t < num_threads; t++) {
        if (pthread_create(&threads[t], NULL, refine_worker, &state) != 0) {
            fprintf(stderr, "WARNING: pthread_create failed for thread %d, "
                    "continuing with %d threads\n", t, t);
            num_threads = t;
            break;
        }
    }

    /* Wait for all workers */
    for (int t = 0; t < num_threads; t++)
        pthread_join(threads[t], NULL);

    free(threads);
    pthread_mutex_destroy(&state.progress_lock);

    /* Count results */
    int refined = 0, lost = 0;
    for (int i = 0; i < cl->num_families; i++) {
        if (cl->families[i].num_instances >= 2) refined++;
        else lost++;
    }

    if (verbose)
        fprintf(stderr, "  Alignment refinement: %d families refined, "
                "%d lost all instances\n", refined, lost);
}
