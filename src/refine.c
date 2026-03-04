#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "refine.h"
#include "align.h"
#include "mdl.h"

/* ================================================================
 * K-mer profile for fast screening
 * ================================================================ */

/* Profile size: 4^8 = 65536 bits = 1024 uint64_t words */
#define PROFILE_BITS  65536
#define PROFILE_WORDS (PROFILE_BITS / 64)

/*
 * Compute k-mer profile (bitset) for a consensus sequence.
 * Each k-mer present sets its corresponding bit.
 */
static void compute_kmer_profile(const char *consensus, int len, int k,
                                 uint64_t *profile)
{
    memset(profile, 0, PROFILE_WORDS * sizeof(uint64_t));
    if (len < k) return;

    for (int i = 0; i <= len - k; i++) {
        uint32_t hash = 0;
        int valid = 1;
        for (int j = 0; j < k; j++) {
            char b = consensus[i + j];
            if (b == DNA_N) { valid = 0; break; }
            hash = (hash << 2) | (uint32_t)b;
        }
        if (!valid) continue;
        profile[hash / 64] |= (uint64_t)1 << (hash % 64);
    }
}

/*
 * Compute Jaccard similarity between two k-mer profiles.
 */
static float kmer_jaccard(const uint64_t *a, const uint64_t *b)
{
    int intersection = 0, union_count = 0;
    for (int i = 0; i < PROFILE_WORDS; i++) {
        intersection += __builtin_popcountll(a[i] & b[i]);
        union_count  += __builtin_popcountll(a[i] | b[i]);
    }
    if (union_count == 0) return 0.0f;
    return (float)intersection / (float)union_count;
}

/* ================================================================
 * Reverse complement of numeric consensus
 * ================================================================ */

static void revcomp_consensus(const char *src, char *dst, int len)
{
    for (int i = 0; i < len; i++)
        dst[len - 1 - i] = dna_complement(src[i]);
}

/* ================================================================
 * Semi-global alignment (free end-gaps on the shorter sequence)
 * ================================================================ */

/*
 * Semi-global DP alignment of query (shorter) against target (longer).
 * Free leading and trailing gaps on the query side.
 * Outputs identity, coverage (of shorter), and aligned bp count.
 * Returns 1 if alignment meets the given thresholds, 0 otherwise.
 */
static int semiglobal_align(const char *query, int qlen,
                            const char *target, int tlen,
                            float min_identity, float min_coverage,
                            int min_aligned,
                            float *identity_out, float *coverage_out,
                            int *aligned_out)
{
    /* Guard against huge DP matrices */
    int64_t dp_cells = (int64_t)(qlen + 1) * (int64_t)(tlen + 1);
    if (dp_cells > REFINE_MAX_DP_CELLS) {
        if (identity_out) *identity_out = 0;
        if (coverage_out) *coverage_out = 0;
        if (aligned_out) *aligned_out = 0;
        return 0;
    }

    /* Allocate DP matrix: (qlen+1) x (tlen+1) */
    int rows = qlen + 1;
    int cols = tlen + 1;
    int *dp = malloc((size_t)rows * (size_t)cols * sizeof(int));
    if (!dp) return 0;

    /* Initialization:
     * dp[0][j] = 0 for all j (free leading gaps on query)
     * dp[i][0] = i * GAP (penalize leading gaps on target) */
    for (int j = 0; j < cols; j++)
        dp[j] = 0;  /* row 0 */
    for (int i = 1; i < rows; i++)
        dp[i * cols] = i * REFINE_GAP;

    /* Fill */
    for (int i = 1; i < rows; i++) {
        for (int j = 1; j < cols; j++) {
            int match_score = (query[i-1] == target[j-1] &&
                               query[i-1] != DNA_N)
                              ? REFINE_MATCH : REFINE_MISMATCH;
            int diag = dp[(i-1) * cols + (j-1)] + match_score;
            int up   = dp[(i-1) * cols + j] + REFINE_GAP;
            int left = dp[i * cols + (j-1)] + REFINE_GAP;
            int best = diag;
            if (up > best) best = up;
            if (left > best) best = left;
            dp[i * cols + j] = best;
        }
    }

    /* Find best score in last row (free trailing gaps on query) */
    int best_score = dp[qlen * cols];
    int best_j = 0;
    for (int j = 1; j < cols; j++) {
        if (dp[qlen * cols + j] > best_score) {
            best_score = dp[qlen * cols + j];
            best_j = j;
        }
    }

    /* Traceback to compute identity statistics */
    int matches = 0, mismatches = 0, gaps = 0;
    int i = qlen, j = best_j;

    while (i > 0 && j > 0) {
        int current = dp[i * cols + j];
        int match_score = (query[i-1] == target[j-1] &&
                           query[i-1] != DNA_N)
                          ? REFINE_MATCH : REFINE_MISMATCH;

        if (current == dp[(i-1) * cols + (j-1)] + match_score) {
            if (match_score == REFINE_MATCH)
                matches++;
            else
                mismatches++;
            i--; j--;
        } else if (current == dp[(i-1) * cols + j] + REFINE_GAP) {
            gaps++;
            i--;
        } else {
            gaps++;
            j--;
        }
    }
    /* Remaining query bases are gaps (leading gaps on target side) */
    gaps += i;

    free(dp);

    int aligned_len = matches + mismatches + gaps;
    if (aligned_len < min_aligned) {
        if (identity_out) *identity_out = 0;
        if (coverage_out) *coverage_out = 0;
        if (aligned_out) *aligned_out = 0;
        return 0;
    }

    float identity = (aligned_len > 0)
                     ? (float)matches / (float)aligned_len : 0.0f;
    float coverage = (qlen > 0)
                     ? (float)(matches + mismatches) / (float)qlen : 0.0f;

    if (identity_out) *identity_out = identity;
    if (coverage_out) *coverage_out = coverage;
    if (aligned_out) *aligned_out = aligned_len;

    return (identity >= min_identity && coverage >= min_coverage &&
            aligned_len >= min_aligned);
}

/* ================================================================
 * Instance overlap detection
 * ================================================================ */

/* Comparator for sorting instances by position */
static int cmp_instance_pos(const void *x, const void *y)
{
    gpos_t pa = ((const Instance *)x)->position;
    gpos_t pb = ((const Instance *)y)->position;
    if (pa < pb) return -1;
    if (pa > pb) return  1;
    return 0;
}

/*
 * Check if instances of two families significantly overlap in the genome.
 * Returns 1 if overlap_frac >= REFINE_OVERLAP_RELAX, 0 otherwise.
 *
 * For large families (>500 instances), uses sorted merge + binary search
 * for O((n+m)log(n+m)) instead of O(n*m).
 */
static int check_instance_overlap(const CandidateFamily *a,
                                  const CandidateFamily *b)
{
    if (a->num_instances == 0 || b->num_instances == 0)
        return 0;

    int overlapping = 0;

    if (a->num_instances > 500 || b->num_instances > 500) {
        /* Sorted approach: sort b instances by position, binary search for each a */
        Instance *b_sorted = malloc((size_t)b->num_instances * sizeof(Instance));
        if (!b_sorted) goto fallback;
        memcpy(b_sorted, b->instances, (size_t)b->num_instances * sizeof(Instance));
        qsort(b_sorted, (size_t)b->num_instances, sizeof(Instance), cmp_instance_pos);

        for (int i = 0; i < a->num_instances; i++) {
            gpos_t a_start = a->instances[i].position;
            gpos_t a_end   = a_start + a->instances[i].aligned_length;

            /* Binary search for first b instance that could overlap */
            int lo = 0, hi = b->num_instances - 1;
            while (lo < hi) {
                int mid = lo + (hi - lo) / 2;
                gpos_t b_end = b_sorted[mid].position + b_sorted[mid].aligned_length;
                if (b_end <= a_start)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            /* Scan forward from lo */
            for (int j = lo; j < b->num_instances; j++) {
                gpos_t b_start = b_sorted[j].position;
                if (b_start >= a_end) break; /* no more can overlap */
                gpos_t b_end2 = b_start + b_sorted[j].aligned_length;

                gpos_t ovl_start = (a_start > b_start) ? a_start : b_start;
                gpos_t ovl_end   = (a_end < b_end2) ? a_end : b_end2;
                if (ovl_end > ovl_start) {
                    glen_t ovl_len = ovl_end - ovl_start;
                    glen_t shorter = a->instances[i].aligned_length;
                    if (b_sorted[j].aligned_length < shorter)
                        shorter = b_sorted[j].aligned_length;
                    if (shorter > 0 && (float)ovl_len / (float)shorter >= 0.50f) {
                        overlapping++;
                        break;
                    }
                }
            }
        }
        free(b_sorted);
    } else {
fallback:
        /* Small families: linear scan (original O(n*m) approach) */
        for (int i = 0; i < a->num_instances; i++) {
            gpos_t a_start = a->instances[i].position;
            gpos_t a_end   = a_start + a->instances[i].aligned_length;

            for (int j = 0; j < b->num_instances; j++) {
                gpos_t b_start = b->instances[j].position;
                gpos_t b_end   = b_start + b->instances[j].aligned_length;

                gpos_t ovl_start = (a_start > b_start) ? a_start : b_start;
                gpos_t ovl_end   = (a_end < b_end) ? a_end : b_end;
                if (ovl_end > ovl_start) {
                    glen_t ovl_len = ovl_end - ovl_start;
                    glen_t shorter = a->instances[i].aligned_length;
                    if (b->instances[j].aligned_length < shorter)
                        shorter = b->instances[j].aligned_length;
                    if (shorter > 0 && (float)ovl_len / (float)shorter >= 0.50f) {
                        overlapping++;
                        break;
                    }
                }
            }
        }
    }

    int min_inst = a->num_instances;
    if (b->num_instances < min_inst) min_inst = b->num_instances;

    return ((float)overlapping / (float)min_inst >= REFINE_OVERLAP_RELAX);
}

/* ================================================================
 * Union-Find
 * ================================================================ */

static int uf_find(int *parent, int x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]]; /* path compression */
        x = parent[x];
    }
    return x;
}

static void uf_unite(int *parent, int *rank, int a, int b)
{
    a = uf_find(parent, a);
    b = uf_find(parent, b);
    if (a == b) return;
    if (rank[a] < rank[b]) { int t = a; a = b; b = t; }
    parent[b] = a;
    if (rank[a] == rank[b]) rank[a]++;
}

/* ================================================================
 * Parallel merge worker
 * ================================================================ */

typedef struct { int i, j; } MergePair;

typedef struct {
    const CandidateList *cl;
    uint64_t           **profiles;
    int                  n;
    int                 *next_row;     /* shared atomic row counter */
    int                  max_cons_len;
    /* Per-thread output */
    MergePair           *pairs;
    int                  num_pairs;
    int                  cap_pairs;
} MergeWorker;

static void *merge_worker_fn(void *arg)
{
    MergeWorker *w = (MergeWorker *)arg;

    /* Thread-local revcomp buffer */
    char *rc_buf = malloc((size_t)(w->max_cons_len + 1));
    if (!rc_buf) return NULL;

    while (1) {
        int i = __atomic_fetch_add(w->next_row, 1, __ATOMIC_SEQ_CST);
        if (i >= w->n) break;

        for (int j = i + 1; j < w->n; j++) {
            float jaccard = kmer_jaccard(w->profiles[i], w->profiles[j]);
            if (jaccard < REFINE_MIN_JACCARD)
                continue;

            const CandidateFamily *fi = &w->cl->families[i];
            const CandidateFamily *fj = &w->cl->families[j];
            const char *query, *target;
            int qlen, tlen;

            if (fi->consensus_length <= fj->consensus_length) {
                query = fi->consensus; qlen = fi->consensus_length;
                target = fj->consensus; tlen = fj->consensus_length;
            } else {
                query = fj->consensus; qlen = fj->consensus_length;
                target = fi->consensus; tlen = fi->consensus_length;
            }

            if (qlen < REFINE_MIN_ALIGNED)
                continue;

            int should_merge = 0;

            /* Check if DP would exceed cell limit for long consensus pairs */
            int64_t dp_cells = (int64_t)(qlen + 1) * (int64_t)(tlen + 1);
            if (dp_cells > REFINE_MAX_DP_CELLS) {
                /* Fallback: Jaccard + instance overlap for long families */
                if (jaccard >= 0.80f && check_instance_overlap(fi, fj))
                    should_merge = 1;
            } else {
                float fwd_identity = 0, fwd_coverage = 0;
                int fwd_aligned = 0;
                int fwd_pass = semiglobal_align(query, qlen, target, tlen,
                                                REFINE_MIN_IDENTITY,
                                                REFINE_MIN_COVERAGE,
                                                REFINE_MIN_ALIGNED,
                                                &fwd_identity, &fwd_coverage,
                                                &fwd_aligned);

                float rc_identity = 0, rc_coverage = 0;
                int rc_aligned = 0;
                revcomp_consensus(target, rc_buf, tlen);
                int rc_pass = semiglobal_align(query, qlen, rc_buf, tlen,
                                               REFINE_MIN_IDENTITY,
                                               REFINE_MIN_COVERAGE,
                                               REFINE_MIN_ALIGNED,
                                               &rc_identity, &rc_coverage,
                                               &rc_aligned);

                float best_identity = fwd_identity;
                float best_coverage = fwd_coverage;
                if (rc_identity > fwd_identity) {
                    best_identity = rc_identity;
                    best_coverage = rc_coverage;
                }

                if (fwd_pass || rc_pass) {
                    should_merge = 1;
                } else if (best_identity >= REFINE_RELAXED_IDENTITY &&
                           best_coverage >= REFINE_RELAXED_COVERAGE) {
                    if (check_instance_overlap(fi, fj))
                        should_merge = 1;
                }
            }

            if (should_merge) {
                if (w->num_pairs >= w->cap_pairs) {
                    int nc = w->cap_pairs == 0 ? 256 : w->cap_pairs * 2;
                    MergePair *tmp = realloc(w->pairs,
                                             (size_t)nc * sizeof(MergePair));
                    if (!tmp) continue;
                    w->pairs = tmp;
                    w->cap_pairs = nc;
                }
                w->pairs[w->num_pairs].i = i;
                w->pairs[w->num_pairs].j = j;
                w->num_pairs++;
            }
        }
    }

    free(rc_buf);
    return NULL;
}

/* ================================================================
 * Main merge orchestrator
 * ================================================================ */

int refine_merge_families(CandidateList *cl, const Genome *genome,
                          const KmerTable *kt, int k, int verbose,
                          int num_threads)
{
    int n = cl->num_families;
    if (n < 2) return 0;

    /* --- Step 1: Compute k-mer profiles for screening --- */
    uint64_t **profiles = malloc((size_t)n * sizeof(uint64_t *));
    if (!profiles) return 0;
    for (int i = 0; i < n; i++) {
        profiles[i] = malloc(PROFILE_WORDS * sizeof(uint64_t));
        if (!profiles[i]) {
            for (int j = 0; j < i; j++) free(profiles[j]);
            free(profiles);
            return 0;
        }
        compute_kmer_profile(cl->families[i].consensus,
                             cl->families[i].consensus_length,
                             REFINE_SCREEN_K, profiles[i]);
    }

    /* --- Step 2: Union-Find setup --- */
    int *parent = malloc((size_t)n * sizeof(int));
    int *uf_rank = malloc((size_t)n * sizeof(int));
    if (!parent || !uf_rank) {
        for (int i = 0; i < n; i++) free(profiles[i]);
        free(profiles); free(parent); free(uf_rank);
        return 0;
    }
    for (int i = 0; i < n; i++) { parent[i] = i; uf_rank[i] = 0; }

    int n_merges = 0;

    int max_cons_len = 0;
    for (int i = 0; i < n; i++)
        if (cl->families[i].consensus_length > max_cons_len)
            max_cons_len = cl->families[i].consensus_length;

    /* --- Step 3: Screen pairs and align --- */
    if (num_threads <= 1) {
        /* Sequential path */
        char *rc_buf = malloc((size_t)(max_cons_len + 1));
        if (!rc_buf) {
            for (int i = 0; i < n; i++) free(profiles[i]);
            free(profiles); free(parent); free(uf_rank);
            return 0;
        }

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (uf_find(parent, i) == uf_find(parent, j))
                    continue;

                float jaccard = kmer_jaccard(profiles[i], profiles[j]);
                if (jaccard < REFINE_MIN_JACCARD)
                    continue;

                CandidateFamily *fi = &cl->families[i];
                CandidateFamily *fj = &cl->families[j];
                const char *query, *target;
                int qlen, tlen;

                if (fi->consensus_length <= fj->consensus_length) {
                    query = fi->consensus; qlen = fi->consensus_length;
                    target = fj->consensus; tlen = fj->consensus_length;
                } else {
                    query = fj->consensus; qlen = fj->consensus_length;
                    target = fi->consensus; tlen = fi->consensus_length;
                }

                if (qlen < REFINE_MIN_ALIGNED) continue;

                /* Check if DP would exceed cell limit for long consensus pairs */
                int64_t dp_cells = (int64_t)(qlen + 1) * (int64_t)(tlen + 1);
                if (dp_cells > REFINE_MAX_DP_CELLS) {
                    /* Fallback: Jaccard + instance overlap for long families */
                    if (jaccard >= 0.80f && check_instance_overlap(fi, fj)) {
                        uf_unite(parent, uf_rank, i, j);
                        n_merges++;
                        if (verbose >= 2)
                            fprintf(stderr, "  Merge (long-fallback): F%d (len=%d) +"
                                    " F%d (len=%d) jaccard=%.2f\n",
                                    fi->id, fi->consensus_length,
                                    fj->id, fj->consensus_length, jaccard);
                    }
                    continue;
                }

                float fwd_identity = 0, fwd_coverage = 0;
                int fwd_aligned = 0;
                int fwd_pass = semiglobal_align(query, qlen, target, tlen,
                    REFINE_MIN_IDENTITY, REFINE_MIN_COVERAGE,
                    REFINE_MIN_ALIGNED, &fwd_identity, &fwd_coverage,
                    &fwd_aligned);

                float rc_identity = 0, rc_coverage = 0;
                int rc_aligned = 0;
                revcomp_consensus(target, rc_buf, tlen);
                int rc_pass = semiglobal_align(query, qlen, rc_buf, tlen,
                    REFINE_MIN_IDENTITY, REFINE_MIN_COVERAGE,
                    REFINE_MIN_ALIGNED, &rc_identity, &rc_coverage,
                    &rc_aligned);

                float best_identity = fwd_identity;
                float best_coverage = fwd_coverage;
                if (rc_identity > fwd_identity) {
                    best_identity = rc_identity;
                    best_coverage = rc_coverage;
                }

                if (fwd_pass || rc_pass) {
                    uf_unite(parent, uf_rank, i, j);
                    n_merges++;
                    if (verbose >= 2)
                        fprintf(stderr, "  Merge: F%d (len=%d, n=%d) +"
                                " F%d (len=%d, n=%d) identity=%.1f%%"
                                " coverage=%.1f%%\n",
                                fi->id, fi->consensus_length,
                                fi->num_instances, fj->id,
                                fj->consensus_length, fj->num_instances,
                                best_identity * 100, best_coverage * 100);
                    continue;
                }

                if (best_identity >= REFINE_RELAXED_IDENTITY &&
                    best_coverage >= REFINE_RELAXED_COVERAGE) {
                    if (check_instance_overlap(fi, fj)) {
                        uf_unite(parent, uf_rank, i, j);
                        n_merges++;
                        if (verbose >= 2)
                            fprintf(stderr, "  Merge (relaxed): F%d + F%d"
                                    " identity=%.1f%% coverage=%.1f%%\n",
                                    fi->id, fj->id,
                                    best_identity * 100,
                                    best_coverage * 100);
                    }
                }
            }
        }
        free(rc_buf);
    } else {
        /* Parallel path: distribute rows across threads */
        int shared_next_row = 0;
        MergeWorker *workers = calloc((size_t)num_threads,
                                       sizeof(MergeWorker));
        pthread_t *threads = malloc((size_t)num_threads * sizeof(pthread_t));
        if (!workers || !threads) {
            free(workers); free(threads);
            /* Fallback to sequential */
            for (int i = 0; i < n; i++) free(profiles[i]);
            free(profiles); free(parent); free(uf_rank);
            return refine_merge_families(cl, genome, kt, k, verbose, 1);
        }

        for (int t = 0; t < num_threads; t++) {
            workers[t].cl = cl;
            workers[t].profiles = profiles;
            workers[t].n = n;
            workers[t].next_row = &shared_next_row;
            workers[t].max_cons_len = max_cons_len;
        }

        int launched = 0;
        for (int t = 0; t < num_threads; t++) {
            if (pthread_create(&threads[t], NULL, merge_worker_fn,
                               &workers[t]) != 0)
                break;
            launched++;
        }
        for (int t = launched; t < num_threads; t++)
            merge_worker_fn(&workers[t]);
        for (int t = 0; t < launched; t++)
            pthread_join(threads[t], NULL);

        /* Apply collected merge pairs to union-find (sequential) */
        for (int t = 0; t < num_threads; t++) {
            for (int p = 0; p < workers[t].num_pairs; p++) {
                int i = workers[t].pairs[p].i;
                int j = workers[t].pairs[p].j;
                if (uf_find(parent, i) != uf_find(parent, j)) {
                    uf_unite(parent, uf_rank, i, j);
                    n_merges++;
                }
            }
            free(workers[t].pairs);
        }
        free(workers);
        free(threads);
    }

    for (int i = 0; i < n; i++) free(profiles[i]);
    free(profiles);

    if (n_merges == 0) {
        free(parent); free(uf_rank);
        return 0;
    }

    /* --- Step 4: Merge groups — transfer instances to representative --- */
    /* For each group, representative = family with most instances */
    int *representative = malloc((size_t)n * sizeof(int));
    if (!representative) { free(parent); free(uf_rank); return n_merges; }

    /* Initialize: each family is its own representative */
    for (int i = 0; i < n; i++)
        representative[i] = i;

    /* Find representative for each group (most instances) */
    for (int i = 0; i < n; i++) {
        int root = uf_find(parent, i);
        if (cl->families[i].num_instances >
            cl->families[representative[root]].num_instances) {
            representative[root] = i;
        }
    }

    /* Transfer instances from non-representative families */
    for (int i = 0; i < n; i++) {
        int root = uf_find(parent, i);
        int rep = representative[root];
        if (i == rep) continue;

        CandidateFamily *dst = &cl->families[rep];
        CandidateFamily *src = &cl->families[i];

        if (src->num_instances == 0) continue;

        /* Grow destination instance array */
        int need = dst->num_instances + src->num_instances;
        if (need > dst->cap_instances) {
            int new_cap = need + 64;
            Instance *tmp = realloc(dst->instances,
                                    (size_t)new_cap * sizeof(Instance));
            if (!tmp) continue; /* skip on alloc failure */
            dst->instances = tmp;
            dst->cap_instances = new_cap;
        }

        memcpy(dst->instances + dst->num_instances,
               src->instances,
               (size_t)src->num_instances * sizeof(Instance));
        dst->num_instances += src->num_instances;

        /* Mark source as absorbed */
        free(src->instances);
        src->instances = NULL;
        src->num_instances = 0;
        src->cap_instances = 0;
    }

    /* --- Step 5: Re-refine merged families --- */
    for (int i = 0; i < n; i++) {
        int root = uf_find(parent, i);
        int rep = representative[root];
        if (i != rep) continue;
        if (i == root && root == rep) {
            /* Only re-refine if this family actually absorbed others */
            int absorbed_any = 0;
            for (int j = 0; j < n; j++) {
                if (j != i && uf_find(parent, j) == root) {
                    absorbed_any = 1;
                    break;
                }
            }
            if (absorbed_any) {
                if (verbose >= 2)
                    fprintf(stderr, "  Re-refining merged F%d (n=%d)\n",
                            cl->families[i].id, cl->families[i].num_instances);
                align_refine_family(&cl->families[i], genome, kt, k,
                                    ALIGN_MAX_ITERATIONS);
            }
        }
    }

    free(representative);
    free(parent);
    free(uf_rank);

    /* --- Step 6: Compact — remove absorbed families --- */
    int write_idx = 0;
    for (int i = 0; i < cl->num_families; i++) {
        if (cl->families[i].num_instances > 0) {
            if (write_idx != i)
                cl->families[write_idx] = cl->families[i];
            cl->families[write_idx].id = (uid_t)write_idx;
            write_idx++;
        } else {
            /* Free consensus of absorbed family */
            free(cl->families[i].consensus);
        }
    }
    cl->num_families = write_idx;

    return n_merges;
}

/* ================================================================
 * Phase 6.3: Family splitting (bimodal divergence detection)
 * ================================================================ */

/*
 * Otsu's method: find the divergence threshold that maximizes
 * inter-class variance, indicating bimodal distribution.
 * Returns the optimal threshold and the bimodality score
 * (inter-class variance / total variance).
 */
static float otsu_threshold(const float *divs, int n, float *bimodality_out,
                            int *valley_pass_out)
{
    *valley_pass_out = 1; /* default: pass */
    if (n < 2) { *bimodality_out = 0; return 0; }

    /* Find range */
    float dmin = divs[0], dmax = divs[0];
    for (int i = 1; i < n; i++) {
        if (divs[i] < dmin) dmin = divs[i];
        if (divs[i] > dmax) dmax = divs[i];
    }
    if (dmax - dmin < 0.01f) { *bimodality_out = 0; return dmin; }

    /* Build histogram */
    int hist[REFINE_DIV_BINS];
    memset(hist, 0, sizeof(hist));
    float bin_width = (dmax - dmin) / (float)REFINE_DIV_BINS;
    for (int i = 0; i < n; i++) {
        int bin = (int)((divs[i] - dmin) / bin_width);
        if (bin >= REFINE_DIV_BINS) bin = REFINE_DIV_BINS - 1;
        if (bin < 0) bin = 0;
        hist[bin]++;
    }

    /* Total mean */
    double total_mean = 0;
    for (int i = 0; i < n; i++) total_mean += divs[i];
    total_mean /= n;

    /* Total variance */
    double total_var = 0;
    for (int i = 0; i < n; i++) {
        double d = divs[i] - total_mean;
        total_var += d * d;
    }
    total_var /= n;
    if (total_var < 1e-10) { *bimodality_out = 0; return (float)total_mean; }

    /* Sweep thresholds to maximize inter-class variance */
    double best_icv = 0;
    int best_bin = 0;
    int w0 = 0;
    double sum0 = 0;
    double sum_total = 0;
    for (int i = 0; i < REFINE_DIV_BINS; i++)
        sum_total += hist[i] * (dmin + (i + 0.5f) * bin_width);

    for (int t = 0; t < REFINE_DIV_BINS - 1; t++) {
        w0 += hist[t];
        int w1 = n - w0;
        if (w0 == 0 || w1 == 0) continue;

        sum0 += hist[t] * (dmin + (t + 0.5f) * bin_width);
        double mu0 = sum0 / w0;
        double mu1 = (sum_total - sum0) / w1;
        double icv = (double)w0 * (double)w1 * (mu0 - mu1) * (mu0 - mu1)
                     / ((double)n * (double)n);

        if (icv > best_icv) {
            best_icv = icv;
            best_bin = t;
        }
    }

    float bimodality = (float)(best_icv / total_var);
    *bimodality_out = bimodality;

    /* Valley depth check for borderline bimodality (0.20-0.40 range).
     * For clear bimodality (>= 0.40), skip this check to protect
     * legitimate splits like Alu S/J subfamilies. */
    if (bimodality >= REFINE_BIMODALITY_THRESH && bimodality < 0.40f) {
        /* Find lower mode height (max bin count in left partition) */
        int lower_mode_height = 0;
        for (int i = 0; i <= best_bin; i++)
            if (hist[i] > lower_mode_height)
                lower_mode_height = hist[i];

        /* Valley height = histogram value at the threshold bin */
        int valley_height = hist[best_bin];

        /* Require valley to be deep enough: valley < lower_mode / 2 */
        if (lower_mode_height > 0 &&
            valley_height >= lower_mode_height / 2) {
            *valley_pass_out = 0; /* valley not deep enough */
        }
    }

    return dmin + (best_bin + 1) * bin_width;
}

/*
 * Build consensus from a subset of instances via majority voting.
 * Allocates and returns a new consensus; caller must free.
 * Returns consensus length, or 0 on failure.
 */
static int build_subset_consensus(const CandidateFamily *fam,
                                  const Instance *insts, int n_insts,
                                  const Genome *genome,
                                  char **consensus_out)
{
    if (n_insts < 2 || fam->consensus_length <= 0) return 0;

    int cons_len = fam->consensus_length;
    char *cons = malloc((size_t)cons_len);
    if (!cons) return 0;
    memcpy(cons, fam->consensus, (size_t)cons_len);

    /* Weighted majority voting over the subset instances */
    for (int p = 0; p < cons_len; p++) {
        int count[4] = {0, 0, 0, 0};

        for (int j = 0; j < n_insts; j++) {
            const Instance *inst = &insts[j];
            if (p < inst->cons_start || p >= inst->cons_end)
                continue;

            /* Weight by alignment score */
            int weight = inst->score;
            if (weight < 1) weight = 1;

            gpos_t gp;
            if (inst->strand > 0)
                gp = inst->position + (p - inst->cons_start);
            else
                gp = inst->position + (inst->cons_end - 1 - p);

            if (gp < 0 || gp >= genome->length) continue;

            char base = genome->sequence[gp];
            if (inst->strand < 0)
                base = dna_complement(base);

            if (base >= 0 && base <= 3)
                count[(int)base] += weight;
        }

        int total = count[0] + count[1] + count[2] + count[3];
        if (total == 0) continue;

        int best_base = 0;
        for (int b = 1; b < 4; b++)
            if (count[b] > count[best_base]) best_base = b;

        cons[p] = (char)best_base;
    }

    *consensus_out = cons;
    return cons_len;
}

int refine_split_families(CandidateList *cl, const Genome *genome,
                          const KmerTable *kt, int k,
                          glen_t genome_len, int verbose,
                          int num_families)
{
    int n_splits = 0;
    int orig_n = cl->num_families;

    for (int fi = 0; fi < orig_n; fi++) {
        CandidateFamily *fam = &cl->families[fi];

        if (fam->num_instances < REFINE_MIN_SPLIT_INSTANCES)
            continue;

        /* Collect divergence values */
        float *divs = malloc((size_t)fam->num_instances * sizeof(float));
        if (!divs) continue;
        for (int j = 0; j < fam->num_instances; j++)
            divs[j] = fam->instances[j].divergence;

        /* Otsu's method to find split threshold */
        float bimodality;
        int valley_pass;
        float threshold = otsu_threshold(divs, fam->num_instances, &bimodality,
                                         &valley_pass);

        if (bimodality < REFINE_BIMODALITY_THRESH) {
            free(divs);
            continue;
        }

        /* Valley depth check: reject borderline splits without clear valley */
        if (!valley_pass) {
            if (verbose >= 2)
                fprintf(stderr, "  Skip split F%d: bimodality=%.2f but "
                        "valley not deep enough\n", fam->id, bimodality);
            free(divs);
            continue;
        }

        /* Split instances into two groups */
        int n_lo = 0, n_hi = 0;
        for (int j = 0; j < fam->num_instances; j++) {
            if (divs[j] <= threshold) n_lo++;
            else n_hi++;
        }

        if (n_lo < REFINE_MIN_CLUSTER_SIZE ||
            n_hi < REFINE_MIN_CLUSTER_SIZE) {
            free(divs);
            continue;
        }

        /* Check mean divergence gap */
        float mean_lo = 0, mean_hi = 0;
        for (int j = 0; j < fam->num_instances; j++) {
            if (divs[j] <= threshold) mean_lo += divs[j];
            else mean_hi += divs[j];
        }
        mean_lo /= n_lo;
        mean_hi /= n_hi;

        if (mean_hi - mean_lo < REFINE_MIN_DIV_GAP) {
            free(divs);
            continue;
        }

        /* Allocate instance arrays for each group */
        Instance *lo_insts = malloc((size_t)n_lo * sizeof(Instance));
        Instance *hi_insts = malloc((size_t)n_hi * sizeof(Instance));
        if (!lo_insts || !hi_insts) {
            free(lo_insts); free(hi_insts); free(divs);
            continue;
        }

        int li = 0, hi = 0;
        for (int j = 0; j < fam->num_instances; j++) {
            if (divs[j] <= threshold)
                lo_insts[li++] = fam->instances[j];
            else
                hi_insts[hi++] = fam->instances[j];
        }
        free(divs);

        /* Build consensus for each group */
        char *lo_cons = NULL, *hi_cons = NULL;
        int lo_clen = build_subset_consensus(fam, lo_insts, n_lo, genome,
                                             &lo_cons);
        int hi_clen = build_subset_consensus(fam, hi_insts, n_hi, genome,
                                             &hi_cons);

        if (lo_clen <= 0 || hi_clen <= 0) {
            free(lo_insts); free(hi_insts);
            free(lo_cons); free(hi_cons);
            continue;
        }

        /* Score original family */
        mdl_score_family(fam, genome_len, num_families);
        double orig_score = fam->mdl_score;

        /* Create temporary sub-families and score them */
        CandidateFamily sub_lo, sub_hi;
        memset(&sub_lo, 0, sizeof(sub_lo));
        memset(&sub_hi, 0, sizeof(sub_hi));

        sub_lo.consensus = lo_cons;
        sub_lo.consensus_length = lo_clen;
        sub_lo.instances = lo_insts;
        sub_lo.num_instances = n_lo;
        sub_lo.cap_instances = n_lo;

        sub_hi.consensus = hi_cons;
        sub_hi.consensus_length = hi_clen;
        sub_hi.instances = hi_insts;
        sub_hi.num_instances = n_hi;
        sub_hi.cap_instances = n_hi;

        /* Use num_families+1 for scoring splits (one more family after split) */
        mdl_score_family(&sub_lo, genome_len, num_families + 1);
        mdl_score_family(&sub_hi, genome_len, num_families + 1);
        double split_score = sub_lo.mdl_score + sub_hi.mdl_score;

        if (split_score <= orig_score) {
            /* Split doesn't improve MDL, reject */
            free(lo_cons); free(hi_cons);
            free(lo_insts); free(hi_insts);
            continue;
        }

        /* Accept split: replace original with sub_lo, add sub_hi as new */
        if (verbose)
            fprintf(stderr, "  Split F%d: %d instances (div %.1f%%) + "
                    "%d instances (div %.1f%%), MDL %.1f -> %.1f\n",
                    fam->id, n_lo, mean_lo * 100, n_hi, mean_hi * 100,
                    orig_score, split_score);

        /* Replace original family with low-divergence group */
        free(fam->consensus);
        free(fam->instances);
        fam->consensus = lo_cons;
        fam->consensus_length = lo_clen;
        fam->instances = lo_insts;
        fam->num_instances = n_lo;
        fam->cap_instances = n_lo;
        fam->mdl_score = sub_lo.mdl_score;
        fam->model_cost = sub_lo.model_cost;

        /* Add high-divergence group as new family */
        if (cl->num_families >= cl->cap_families) {
            int new_cap = cl->cap_families * 2;
            CandidateFamily *tmp = realloc(cl->families,
                                           (size_t)new_cap * sizeof(CandidateFamily));
            if (!tmp) { free(hi_cons); free(hi_insts); continue; }
            cl->families = tmp;
            cl->cap_families = new_cap;
            /* Re-acquire pointer after realloc */
            fam = &cl->families[fi];
        }

        CandidateFamily *new_fam = &cl->families[cl->num_families];
        memset(new_fam, 0, sizeof(CandidateFamily));
        new_fam->id = (uid_t)cl->num_families;
        new_fam->consensus = hi_cons;
        new_fam->consensus_length = hi_clen;
        new_fam->instances = hi_insts;
        new_fam->num_instances = n_hi;
        new_fam->cap_instances = n_hi;
        new_fam->component_id = fam->component_id;
        new_fam->topology = fam->topology;
        new_fam->estimated_copies = (freq_t)n_hi;
        new_fam->mdl_score = sub_hi.mdl_score;
        new_fam->model_cost = sub_hi.model_cost;
        cl->num_families++;

        /* Re-refine both split families */
        align_refine_family(&cl->families[fi], genome, kt, k,
                            ALIGN_MAX_ITERATIONS);
        align_refine_family(&cl->families[cl->num_families - 1], genome, kt, k,
                            ALIGN_MAX_ITERATIONS);

        n_splits++;
    }

    return n_splits;
}

/* ================================================================
 * Phase 6.4: Pruning marginal families after MDL selection
 * ================================================================ */

int refine_prune_families(CandidateList *cl, glen_t genome_len, int verbose,
                          int num_families)
{
    /* Count accepted families */
    int n_accepted = 0;
    for (int i = 0; i < cl->num_families; i++)
        if (cl->families[i].mdl_score > 0) n_accepted++;

    if (n_accepted <= 1) return 0;

    /* Build (index, score) pairs for accepted families, sorted ascending */
    typedef struct { int idx; double score; } IdxScore;
    IdxScore *order = malloc((size_t)n_accepted * sizeof(IdxScore));
    if (!order) return 0;

    int oi = 0;
    for (int i = 0; i < cl->num_families; i++) {
        if (cl->families[i].mdl_score > 0) {
            order[oi].idx = i;
            order[oi].score = cl->families[i].mdl_score;
            oi++;
        }
    }

    /* Sort by score ascending (weakest first) */
    for (int i = 0; i < n_accepted - 1; i++) {
        for (int j = i + 1; j < n_accepted; j++) {
            if (order[j].score < order[i].score) {
                IdxScore tmp = order[i];
                order[i] = order[j];
                order[j] = tmp;
            }
        }
    }

    /* Build coverage count array */
    uint8_t *cov = calloc((size_t)genome_len, 1);
    if (!cov) { free(order); return 0; }

    for (int i = 0; i < cl->num_families; i++) {
        if (cl->families[i].mdl_score <= 0) continue;
        CandidateFamily *f = &cl->families[i];
        for (int j = 0; j < f->num_instances; j++) {
            gpos_t start = f->instances[j].position;
            int alen = (int)f->instances[j].aligned_length;
            for (int p = 0; p < alen; p++) {
                gpos_t gp = start + p;
                if (gp >= 0 && gp < genome_len && cov[gp] < 255)
                    cov[gp]++;
            }
        }
    }

    int n_pruned = 0;

    /* Try pruning each accepted family, weakest first */
    for (int k = 0; k < n_accepted; k++) {
        int fi = order[k].idx;
        CandidateFamily *f = &cl->families[fi];
        if (f->mdl_score <= 0) continue; /* already pruned in this pass */

        /* Compute exclusive savings: only count positions where cov == 1 */
        double exclusive_savings = 0;
        int exclusive_instances = 0;

        for (int j = 0; j < f->num_instances; j++) {
            Instance *inst = &f->instances[j];
            gpos_t start = inst->position;
            int alen = (int)inst->aligned_length;

            int excl_bases = 0;
            for (int p = 0; p < alen; p++) {
                gpos_t gp = start + p;
                if (gp >= 0 && gp < genome_len && cov[gp] == 1)
                    excl_bases++;
            }

            /* Skip instances with <25% exclusive coverage */
            if (excl_bases < alen / 4) continue;

            int edits = (int)(inst->divergence * excl_bases + 0.5f);
            double lit = 2.0 * excl_bases;
            double enc = mdl_instance_cost_full(excl_bases, edits,
                                                f->consensus_length,
                                                num_families);
            exclusive_savings += (lit - enc);
            exclusive_instances++;
        }

        double exclusive_score = exclusive_savings - f->model_cost;

        if (exclusive_instances == 0) {
            /* No exclusive coverage at all — purely redundant */
            if (verbose)
                fprintf(stderr, "  Pruned F%d: excl_score=%.1f "
                        "(model=%.1f, excl_inst=%d)\n",
                        f->id, exclusive_score, f->model_cost,
                        exclusive_instances);

            /* Decrement coverage counts */
            for (int j = 0; j < f->num_instances; j++) {
                gpos_t start = f->instances[j].position;
                int alen = (int)f->instances[j].aligned_length;
                for (int p = 0; p < alen; p++) {
                    gpos_t gp = start + p;
                    if (gp >= 0 && gp < genome_len && cov[gp] > 0)
                        cov[gp]--;
                }
            }

            f->mdl_score = 0;
            n_pruned++;
        }
    }

    free(cov);
    free(order);
    return n_pruned;
}

/* ================================================================
 * Phase 6b: Fragment assembly (assemble adjacent TE fragments)
 *
 * Non-overlapping fragments from the same TE share zero k-mers
 * (they're from different parts of the consensus). Use spatial
 * co-occurrence instead: if instances of two families consistently
 * appear near each other in the genome with consistent orientation,
 * they're likely fragments of a single TE.
 * ================================================================ */

/* Instance position entry for sweep-line */
typedef struct {
    gpos_t start;
    gpos_t end;
    int    family_idx;
    int    instance_idx;
    int8_t strand;
} InstanceEntry;

/* Co-occurrence pair accumulator */
typedef struct {
    int fam_a;
    int fam_b;
    int count;           /* number of co-occurrences */
    int same_dir_count;  /* number with consistent relative orientation */
    int a_inside_b;      /* instances of A spatially inside B */
    int b_inside_a;      /* instances of B spatially inside A */
} CooccPair;

/* Gap observation from one co-occurring pair */
typedef struct {
    int gap;             /* signed: positive = gap, negative = overlap */
    int8_t a_strand;     /* strand of family A instance */
} GapObs;

static int cmp_instance_entry(const void *x, const void *y)
{
    gpos_t sa = ((const InstanceEntry *)x)->start;
    gpos_t sb = ((const InstanceEntry *)y)->start;
    if (sa < sb) return -1;
    if (sa > sb) return  1;
    return 0;
}

static int cmp_int(const void *a, const void *b)
{
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    return (ia > ib) - (ia < ib);
}


int refine_assemble_fragments(CandidateList *cl, const Genome *genome,
                              const KmerTable *kt, int k,
                              glen_t genome_len, int verbose,
                              int num_threads)
{
    (void)num_threads; /* future: parallelize sweep */
    (void)genome;      /* used by fill_gap_majority if enabled */
    (void)kt;
    (void)k;
    int n = cl->num_families;
    if (n < 2) return 0;

    /* ---- Step 1: Build sorted instance entry array ---- */
    int total_instances = 0;
    for (int i = 0; i < n; i++)
        total_instances += cl->families[i].num_instances;

    if (total_instances == 0) return 0;

    /* Subsample: cap each family at 10K instances for efficiency */
    int max_per_family = 10000;
    int entries_count = 0;
    for (int i = 0; i < n; i++) {
        int ni = cl->families[i].num_instances;
        entries_count += (ni > max_per_family) ? max_per_family : ni;
    }

    InstanceEntry *entries = malloc((size_t)entries_count * sizeof(InstanceEntry));
    if (!entries) return 0;

    int ei = 0;
    for (int i = 0; i < n; i++) {
        CandidateFamily *f = &cl->families[i];
        int ni = f->num_instances;
        int step = 1;
        if (ni > max_per_family)
            step = ni / max_per_family;

        for (int j = 0; j < ni; j += step) {
            if (ei >= entries_count) break;
            Instance *inst = &f->instances[j];
            entries[ei].start = inst->position;
            entries[ei].end = inst->position + inst->aligned_length;
            entries[ei].family_idx = i;
            entries[ei].instance_idx = j;
            entries[ei].strand = inst->strand;
            ei++;
        }
    }
    entries_count = ei;

    qsort(entries, (size_t)entries_count, sizeof(InstanceEntry),
          cmp_instance_entry);

    /* ---- Step 2: Compute proximity distance D ---- */
    /* D = median_consensus_length * 2, clamped to [500, 5000] */
    int64_t cons_len_sum = 0;
    for (int i = 0; i < n; i++)
        cons_len_sum += cl->families[i].consensus_length;
    int median_cons = (int)(cons_len_sum / n);
    int D = median_cons * 2;
    if (D < 500) D = 500;
    if (D > 5000) D = 5000;

    /* ---- Step 3: Sweep-line co-occurrence counting ---- */
    /* Use a hash map of pair -> count. For simplicity, use a flat array
     * for pairs since n is moderate (typically < 100K families post-merge). */
    /* For efficiency, limit to pairs with n_families <= 50K */
    if (n > 50000) {
        free(entries);
        if (verbose)
            fprintf(stderr, "  Fragment assembly: too many families (%d), skipping\n", n);
        return 0;
    }

    /* Pair tracking: use dynamic array of CooccPair, hashed by (min, max) */
    int cap_pairs = 1024;
    int num_pairs = 0;
    CooccPair *pairs = calloc((size_t)cap_pairs, sizeof(CooccPair));
    if (!pairs) { free(entries); return 0; }

    /* Hash table for pair lookup: map (fam_a, fam_b) -> index in pairs[] */
    /* Simple open-addressing hash */
    int ht_size = 16 * 1024;  /* must be power of 2 */
    int *ht = malloc((size_t)ht_size * sizeof(int));
    if (!ht) { free(entries); free(pairs); return 0; }
    memset(ht, -1, (size_t)ht_size * sizeof(int));

    #define PAIR_HASH(a, b) (((unsigned)(a) * 2654435761u ^ (unsigned)(b) * 40503u) & (unsigned)(ht_size - 1))

    for (int i = 0; i < entries_count; i++) {
        int fa = entries[i].family_idx;

        /* Look forward within distance D */
        for (int j = i + 1; j < entries_count; j++) {
            if (entries[j].start - entries[i].start > D) break;
            int fb = entries[j].family_idx;
            if (fa == fb) continue;

            /* Canonicalize pair order */
            int pa = fa < fb ? fa : fb;
            int pb = fa < fb ? fb : fa;

            /* Check containment for nesting detection */
            int a_in_b = (entries[i].start >= entries[j].start &&
                          entries[i].end <= entries[j].end);
            int b_in_a = (entries[j].start >= entries[i].start &&
                          entries[j].end <= entries[i].end);

            /* Check orientation consistency */
            int same_dir = (entries[i].strand == entries[j].strand);

            /* Find or create pair in hash table */
            unsigned h = PAIR_HASH(pa, pb);
            int found = -1;
            for (int probe = 0; probe < ht_size; probe++) {
                unsigned idx = (h + (unsigned)probe) & (unsigned)(ht_size - 1);
                if (ht[idx] == -1) {
                    /* Empty slot — create new pair */
                    if (num_pairs >= cap_pairs) {
                        int nc = cap_pairs * 2;
                        CooccPair *tmp = realloc(pairs, (size_t)nc * sizeof(CooccPair));
                        if (!tmp) break;
                        pairs = tmp;
                        cap_pairs = nc;
                    }
                    ht[idx] = num_pairs;
                    memset(&pairs[num_pairs], 0, sizeof(CooccPair));
                    pairs[num_pairs].fam_a = pa;
                    pairs[num_pairs].fam_b = pb;
                    found = num_pairs;
                    num_pairs++;
                    break;
                }
                if (pairs[ht[idx]].fam_a == pa && pairs[ht[idx]].fam_b == pb) {
                    found = ht[idx];
                    break;
                }
            }

            if (found >= 0) {
                pairs[found].count++;
                if (same_dir) pairs[found].same_dir_count++;
                if (fa < fb) {
                    if (a_in_b) pairs[found].a_inside_b++;
                    if (b_in_a) pairs[found].b_inside_a++;
                } else {
                    if (a_in_b) pairs[found].b_inside_a++;
                    if (b_in_a) pairs[found].a_inside_b++;
                }
            }
        }
    }

    #undef PAIR_HASH
    free(ht);

    /* ---- Step 4: Classify pairs and filter ---- */
    /* Identify ADJACENT pairs for assembly (not NESTED) */
    typedef struct {
        int fam_a;
        int fam_b;
        int count;
    } AssemblyCandidate;

    int cap_candidates = 256;
    int num_candidates = 0;
    AssemblyCandidate *candidates_arr = malloc((size_t)cap_candidates *
                                                sizeof(AssemblyCandidate));
    if (!candidates_arr) { free(entries); free(pairs); return 0; }

    for (int i = 0; i < num_pairs; i++) {
        CooccPair *p = &pairs[i];

        /* Require ≥3 co-occurrences */
        if (p->count < 3) continue;

        /* Nesting guard: skip if ≥50% of one is inside the other */
        CandidateFamily *fa = &cl->families[p->fam_a];
        CandidateFamily *fb = &cl->families[p->fam_b];
        int min_inst = fa->num_instances;
        if (fb->num_instances < min_inst) min_inst = fb->num_instances;

        if (min_inst > 0) {
            float frac_a_in_b = (float)p->a_inside_b / (float)min_inst;
            float frac_b_in_a = (float)p->b_inside_a / (float)min_inst;
            if (frac_a_in_b >= 0.50f || frac_b_in_a >= 0.50f) {
                if (verbose >= 2)
                    fprintf(stderr, "  Nesting detected: F%d/F%d "
                            "(a_in_b=%.0f%%, b_in_a=%.0f%%)\n",
                            fa->id, fb->id,
                            frac_a_in_b * 100, frac_b_in_a * 100);
                continue;
            }
        }

        /* Size ratio guard: skip if one is much shorter */
        float size_ratio = (fa->consensus_length < fb->consensus_length)
            ? (float)fa->consensus_length / (float)fb->consensus_length
            : (float)fb->consensus_length / (float)fa->consensus_length;
        if (size_ratio < 0.10f) continue;

        /* Orientation consistency: ≥80% same direction */
        if (p->count > 0 &&
            (float)p->same_dir_count / (float)p->count < 0.80f)
            continue;

        /* Passed all filters — assembly candidate */
        if (num_candidates >= cap_candidates) {
            int nc = cap_candidates * 2;
            AssemblyCandidate *tmp = realloc(candidates_arr,
                                              (size_t)nc * sizeof(AssemblyCandidate));
            if (!tmp) continue;
            candidates_arr = tmp;
            cap_candidates = nc;
        }
        candidates_arr[num_candidates].fam_a = p->fam_a;
        candidates_arr[num_candidates].fam_b = p->fam_b;
        candidates_arr[num_candidates].count = p->count;
        num_candidates++;
    }

    free(pairs);

    if (num_candidates == 0) {
        free(entries);
        free(candidates_arr);
        return 0;
    }

    if (verbose)
        fprintf(stderr, "  Fragment assembly: %d candidate pairs from %d co-occurring\n",
                num_candidates, num_pairs);

    /* ---- Step 5: Union-Find for transitive assembly ---- */
    int *uf_parent = malloc((size_t)n * sizeof(int));
    int *uf_rnk = malloc((size_t)n * sizeof(int));
    if (!uf_parent || !uf_rnk) {
        free(entries); free(candidates_arr);
        free(uf_parent); free(uf_rnk);
        return 0;
    }
    for (int i = 0; i < n; i++) { uf_parent[i] = i; uf_rnk[i] = 0; }

    /* Track chain length per group to cap at 10 */
    int *chain_len = calloc((size_t)n, sizeof(int));
    if (!chain_len) {
        free(entries); free(candidates_arr);
        free(uf_parent); free(uf_rnk);
        return 0;
    }
    for (int i = 0; i < n; i++) chain_len[i] = 1;

    int n_assemblies = 0;

    /* ---- Step 6: For each candidate pair, compute gap and assemble ---- */
    for (int ci = 0; ci < num_candidates; ci++) {
        int fa_idx = candidates_arr[ci].fam_a;
        int fb_idx = candidates_arr[ci].fam_b;

        /* Check union-find: don't merge if already in same group */
        int root_a = uf_find(uf_parent, fa_idx);
        int root_b = uf_find(uf_parent, fb_idx);
        if (root_a == root_b) continue;

        /* Cap chain length */
        if (chain_len[root_a] + chain_len[root_b] > 10) continue;

        CandidateFamily *fa = &cl->families[fa_idx];
        CandidateFamily *fb = &cl->families[fb_idx];

        /* Skip dead families (absorbed by prior assembly in this pass) */
        if (fa->consensus == NULL || fb->consensus == NULL) continue;
        if (fa->num_instances < 2 || fb->num_instances < 2) continue;

        /* Compute gap distribution from co-occurring instances.
         * Find instances of A and B that are near each other. */
        int cap_gaps = 256;
        int n_gaps = 0;
        int *gap_arr = malloc((size_t)cap_gaps * sizeof(int));
        if (!gap_arr) continue;

        /* Linear scan of sorted entries to find nearby A/B pairs.
         * Check both orderings: A before B, and B before A. */
        for (int i = 0; i < entries_count; i++) {
            int fi_fam = entries[i].family_idx;
            if (fi_fam != fa_idx && fi_fam != fb_idx) continue;

            int target_fam = (fi_fam == fa_idx) ? fb_idx : fa_idx;

            for (int j = i + 1; j < entries_count; j++) {
                if (entries[j].start - entries[i].start > D) break;
                if (entries[j].family_idx != target_fam) continue;

                /* Gap: end of first to start of second */
                int gap = (int)(entries[j].start - entries[i].end);

                /* Only count reasonable gaps */
                if (gap < -500 || gap > D) continue;

                if (n_gaps >= cap_gaps) {
                    int nc = cap_gaps * 2;
                    int *tmp = realloc(gap_arr, (size_t)nc * sizeof(int));
                    if (!tmp) break;
                    gap_arr = tmp;
                    cap_gaps = nc;
                }
                gap_arr[n_gaps++] = gap;
            }
        }

        if (n_gaps < 3) { free(gap_arr); continue; }

        /* Compute median gap */
        qsort(gap_arr, (size_t)n_gaps, sizeof(int), cmp_int);
        int median_gap = gap_arr[n_gaps / 2];

        /* Compute MAD for gap consistency check */
        int *abs_devs = malloc((size_t)n_gaps * sizeof(int));
        if (!abs_devs) { free(gap_arr); continue; }
        for (int i = 0; i < n_gaps; i++) {
            int dev = gap_arr[i] - median_gap;
            abs_devs[i] = (dev < 0) ? -dev : dev;
        }
        qsort(abs_devs, (size_t)n_gaps, sizeof(int), cmp_int);
        int mad = abs_devs[n_gaps / 2];
        free(abs_devs);
        free(gap_arr);

        /* Max allowable gap: median + 2*MAD, clamped to 500 */
        int max_gap = median_gap + 2 * mad;
        if (max_gap > 500) max_gap = 500;

        /* Reject if gap too large or too negative (TSD overlap limit) */
        if (median_gap > max_gap || median_gap < -20) continue;

        /* Reject if gap distribution is too noisy (MAD > 100) */
        if (mad > 100) continue;

        /* ---- Build assembled consensus ---- */
        /* Determine order: which family goes first?
         * Use the order that was observed more frequently. */
        int a_first_count = 0, b_first_count = 0;
        for (int i = 0; i < entries_count; i++) {
            if (entries[i].family_idx != fa_idx) continue;
            for (int j = i + 1; j < entries_count; j++) {
                if (entries[j].start - entries[i].start > D) break;
                if (entries[j].family_idx == fb_idx) a_first_count++;
            }
        }
        for (int i = 0; i < entries_count; i++) {
            if (entries[i].family_idx != fb_idx) continue;
            for (int j = i + 1; j < entries_count; j++) {
                if (entries[j].start - entries[i].start > D) break;
                if (entries[j].family_idx == fa_idx) b_first_count++;
            }
        }

        /* Swap so that "first" is the one that tends to come first */
        CandidateFamily *first_fam = fa;
        CandidateFamily *second_fam = fb;
        if (b_first_count > a_first_count) {
            first_fam = fb;
            second_fam = fa;
        }

        int gap_fill_len;  /* actual gap/overlap in assembled consensus */
        int new_cons_len;
        char *new_cons;

        if (median_gap <= 0) {
            /* Overlap or abutting */
            int overlap = -median_gap;
            if (overlap > first_fam->consensus_length ||
                overlap > second_fam->consensus_length)
                continue;

            gap_fill_len = 0;  /* no fill for overlap */
            new_cons_len = first_fam->consensus_length +
                           second_fam->consensus_length - overlap;
            new_cons = malloc((size_t)new_cons_len);
            if (!new_cons) continue;

            memcpy(new_cons, first_fam->consensus,
                   (size_t)first_fam->consensus_length);
            memcpy(new_cons + first_fam->consensus_length,
                   second_fam->consensus + overlap,
                   (size_t)(second_fam->consensus_length - overlap));
        } else {
            /* Positive gap: N-padding */
            gap_fill_len = median_gap;
            new_cons_len = first_fam->consensus_length + gap_fill_len +
                           second_fam->consensus_length;
            new_cons = malloc((size_t)new_cons_len);
            if (!new_cons) continue;

            memcpy(new_cons, first_fam->consensus,
                   (size_t)first_fam->consensus_length);
            memset(new_cons + first_fam->consensus_length,
                   DNA_N, (size_t)gap_fill_len);
            memcpy(new_cons + first_fam->consensus_length + gap_fill_len,
                   second_fam->consensus,
                   (size_t)second_fam->consensus_length);
        }

        /* ---- MDL validation: assembled must beat sum of parts ---- */
        int R_est = cl->num_families;
        if (R_est < 2) R_est = 2;
        mdl_score_family(fa, genome_len, R_est);
        mdl_score_family(fb, genome_len, R_est);
        double sum_individual = fa->mdl_score + fb->mdl_score;

        /* Create temporary assembled family with merged instances.
         * Adjust B's instance coordinates for offset in assembled consensus. */
        int b_offset = first_fam->consensus_length + gap_fill_len;
        if (median_gap <= 0) b_offset = first_fam->consensus_length;

        int merged_n = first_fam->num_instances + second_fam->num_instances;
        Instance *merged_inst = malloc((size_t)merged_n * sizeof(Instance));
        if (!merged_inst) { free(new_cons); continue; }

        /* Copy first family's instances (coordinates unchanged) */
        memcpy(merged_inst, first_fam->instances,
               (size_t)first_fam->num_instances * sizeof(Instance));

        /* Copy second family's instances with adjusted cons_start/cons_end */
        for (int j = 0; j < second_fam->num_instances; j++) {
            Instance inst = second_fam->instances[j];
            inst.cons_start += b_offset;
            inst.cons_end   += b_offset;
            if (inst.cons_end > new_cons_len)
                inst.cons_end = new_cons_len;
            merged_inst[first_fam->num_instances + j] = inst;
        }

        CandidateFamily assembled;
        memset(&assembled, 0, sizeof(assembled));
        assembled.consensus = new_cons;
        assembled.consensus_length = new_cons_len;
        assembled.instances = merged_inst;
        assembled.num_instances = merged_n;
        assembled.cap_instances = merged_n;

        /* Score assembled family (no align_refine — use merged instances) */
        mdl_score_family(&assembled, genome_len, R_est - 1);

        if (assembled.mdl_score <= sum_individual) {
            free(assembled.consensus);
            free(assembled.instances);
            continue;
        }

        /* Accept assembly */
        if (verbose)
            fprintf(stderr, "  Assembled F%d (len=%d, n=%d) + F%d (len=%d, n=%d) "
                    "-> len=%d n=%d gap=%d MDL: %.1f+%.1f=%.1f -> %.1f\n",
                    first_fam->id, first_fam->consensus_length,
                    first_fam->num_instances,
                    second_fam->id, second_fam->consensus_length,
                    second_fam->num_instances,
                    assembled.consensus_length, assembled.num_instances,
                    median_gap, fa->mdl_score, fb->mdl_score,
                    sum_individual, assembled.mdl_score);

        /* Replace fa with assembled, mark fb as absorbed */
        free(fa->consensus);
        free(fa->instances);
        fa->consensus = assembled.consensus;
        fa->consensus_length = assembled.consensus_length;
        fa->instances = assembled.instances;
        fa->num_instances = assembled.num_instances;
        fa->cap_instances = assembled.cap_instances;
        fa->mdl_score = assembled.mdl_score;
        fa->model_cost = assembled.model_cost;

        /* Mark fb as dead */
        free(fb->consensus);
        free(fb->instances);
        fb->consensus = NULL;
        fb->instances = NULL;
        fb->num_instances = 0;
        fb->cap_instances = 0;
        fb->consensus_length = 0;

        /* Update union-find */
        int new_chain = chain_len[root_a] + chain_len[root_b];
        uf_unite(uf_parent, uf_rnk, fa_idx, fb_idx);
        int new_root = uf_find(uf_parent, fa_idx);
        chain_len[new_root] = new_chain;

        n_assemblies++;
    }

    free(entries);
    free(candidates_arr);
    free(chain_len);

    /* Compact: remove dead families */
    if (n_assemblies > 0) {
        int write_idx = 0;
        for (int i = 0; i < cl->num_families; i++) {
            if (cl->families[i].num_instances > 0 &&
                cl->families[i].consensus != NULL) {
                if (write_idx != i)
                    cl->families[write_idx] = cl->families[i];
                cl->families[write_idx].id = (uid_t)write_idx;
                write_idx++;
            }
        }
        cl->num_families = write_idx;
    }

    free(uf_parent);
    free(uf_rnk);

    return n_assemblies;
}
