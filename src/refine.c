#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "refine.h"
#include "align.h"
#include "mdl.h"

/* Runtime-configurable max DP cells (default 10M ~40MB) */
int64_t g_refine_max_dp_cells = 10 * 1000 * 1000;
const char *g_refine_split_audit_path = NULL;

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
    if (dp_cells > g_refine_max_dp_cells) {
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
 * Compute the fraction of SHORTER family's instances fully contained
 * (≥80% of the shorter instance overlaps a longer instance) within
 * LONGER family's instances.
 *
 * Used as the nested-element merge gate (P1 follow-up).  When two
 * consensi pass 80-80-80 at 100% identity / 100% coverage and the
 * shorter is fully embedded in the longer, the question is biological:
 *   - If most of the shorter family's instances are nested in the
 *     longer's → the shorter is just a fragment view of the longer →
 *     merge.
 *   - If many of the shorter's instances are SOLO (no overlapping
 *     longer instance) → the shorter is a distinct family that
 *     happens to be sequence-similar to a contained region of the
 *     longer (e.g. SINE inside LINE; the SINE also exists solo) →
 *     DO NOT merge.
 *
 * Returns: containment fraction in [0, 1].  Computed in O(n_short ×
 * n_long) which is fine since merge candidates are already filtered.
 */
static float nested_containment_fraction(const CandidateFamily *shorter,
                                         const CandidateFamily *longer)
{
    if (shorter->num_instances == 0) return 0.0f;
    int contained = 0;
    int solo      = 0;  /* shorter instances NOT overlapping any longer inst */
    for (int i = 0; i < shorter->num_instances; i++) {
        gpos_t s_start = shorter->instances[i].position;
        gpos_t s_end   = s_start + shorter->instances[i].aligned_length;
        glen_t s_len   = shorter->instances[i].aligned_length;
        if (s_len <= 0) continue;
        int matched = 0;
        int any_overlap = 0;
        for (int j = 0; j < longer->num_instances; j++) {
            /* B+Q6 / ENG-N11: cross-chromosome guard.  Without this,
             * positions on different sequences (which share the same
             * gpos_t coordinate origin only by accident) are compared
             * as if they lived on a single concatenated genome. */
            if (shorter->instances[i].seq_index != longer->instances[j].seq_index)
                continue;
            gpos_t l_start = longer->instances[j].position;
            gpos_t l_end   = l_start + longer->instances[j].aligned_length;
            gpos_t ov_s = (s_start > l_start) ? s_start : l_start;
            gpos_t ov_e = (s_end   < l_end)   ? s_end   : l_end;
            if (ov_e <= ov_s) continue;
            any_overlap = 1;
            glen_t ov_len = ov_e - ov_s;
            if ((float)ov_len / (float)s_len >= 0.80f) {
                matched = 1;
                break;
            }
        }
        if (matched) contained++;
        else if (!any_overlap) solo++;
    }
    /* B (Tier 1.5a): "solo evidence" tightening.  When the shorter family
     * has FEWER than 2 instances outside any longer instance, the two
     * families look indistinguishable in genomic occurrence — likely a
     * structural prefix/suffix relationship (longer = X + shorter, with
     * shorter never living on its own).  In that case, do NOT let the
     * containment-fraction gate block the merge: report full containment
     * so the caller's `ctf < 0.50f` veto does not fire.  We only veto
     * (return low containment) when there is meaningful solo evidence
     * (≥ 2 truly solo shorter instances) AND most shorter instances are
     * not nested.  This catches the SINE-in-LINE case (SINE has many
     * solo instances) while letting strict-prefix/suffix structural
     * matches merge cleanly. */
    if (solo < 2) return 1.0f;
    return (float)contained / (float)shorter->num_instances;
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
                /* ENG-N11 cross-chromosome guard: positions sorted purely
                 * by gpos_t mix multiple sequences; skip any pair on
                 * different chromosomes before computing overlap. */
                if (a->instances[i].seq_index != b_sorted[j].seq_index) continue;
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
                /* ENG-N11 cross-chromosome guard: skip pairs on
                 * different sequences before computing overlap. */
                if (a->instances[i].seq_index != b->instances[j].seq_index) continue;
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
 * MDL gate for merge candidate pairs
 * ================================================================ */

/*
 * Estimate the MDL score of merging two candidate families anchored on
 * the longer family's consensus.  When two families with consensus
 * identity I are merged, instances from each family are re-encoded
 * against the merged (centroid) consensus; both halves of the merged
 * family pay an estimated half_extra = (1-I)/2 additional divergence
 * (since the centroid sits roughly equidistant between the two original
 * consensi).
 *
 * Returns the estimated MDL score (savings - model_cost) of the merged
 * family.  A value <= 0 indicates the merger would produce a family
 * that fails MDL on its own, in which case the union-find merger should
 * be vetoed: otherwise both originals are absorbed into a non-viable
 * merged family that mdl_select_library will then reject, losing both.
 *
 * This is an approximation — the true post-merger consensus is computed
 * by align_refine_family after the merge — but it is fast (no DP) and
 * conservative enough to catch obvious union-find runaway chains.
 */
static double estimate_merge_score(const CandidateFamily *a,
                                   const CandidateFamily *b,
                                   float identity, int num_families)
{
    if (a->num_instances <= 0 && b->num_instances <= 0) return 0.0;

    int big_len = (a->consensus_length >= b->consensus_length)
                  ? a->consensus_length : b->consensus_length;
    if (big_len <= 0) return 0.0;

    double model_cost = mdl_model_cost(big_len);
    double savings = 0.0;

    /* Both halves pay roughly half of the consensus-pair distance (1-I)
     * as additional divergence relative to the centroid consensus. */
    double half_extra = ((double)(1.0f - identity)) / 2.0;
    if (half_extra < 0.0) half_extra = 0.0;

    const CandidateFamily *halves[2] = {a, b};
    for (int h = 0; h < 2; h++) {
        const CandidateFamily *fam = halves[h];
        for (int i = 0; i < fam->num_instances; i++) {
            const Instance *inst = &fam->instances[i];
            int a_i = (int)inst->aligned_length;
            if (a_i <= 0) continue;

            double new_div = (double)inst->divergence + half_extra;
            if (new_div < 0.0) new_div = 0.0;
            if (new_div > 0.5) new_div = 0.5;
            int m_i = (int)(new_div * (double)a_i + 0.5);
            if (m_i < 0) m_i = 0;
            if (m_i > a_i) m_i = a_i;

            double lit = 2.0 * (double)a_i;
            double enc = mdl_instance_cost_full(a_i, m_i, big_len, num_families);
            savings += (lit - enc);
        }
    }

    return savings - model_cost;
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
    int                  num_families;  /* R estimate for MDL gate */
    int                  verbose;
    int                  mdl_vetoes;    /* per-thread veto count */
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

            /* Length-ratio guard (Stage B fix 2).
             * 80-80-80 mutual coverage permits merging vastly-different-
             * length families (e.g. 500 vs 700 bp).  The merged consensus
             * splits the difference; both sides' instances get trimmed at
             * the consensus boundary, costing ~2.36 Mb on chr4 (REFINE_
             * TRACE_REPORT bottleneck #2).  Reject when shorter / longer
             * < 0.7. */
            if ((float)qlen / (float)tlen < 0.7f)
                continue;

            int   should_merge = 0;
            float gate_identity = 0.0f;  /* best identity seen for MDL gate */

            /* Check if DP would exceed cell limit for long consensus pairs */
            int64_t dp_cells = (int64_t)(qlen + 1) * (int64_t)(tlen + 1);
            if (dp_cells > g_refine_max_dp_cells) {
                /* Fallback: Jaccard + instance overlap for long families.
                 * Use jaccard as a rough identity proxy for the MDL gate. */
                if (jaccard >= 0.80f && check_instance_overlap(fi, fj)) {
                    should_merge = 1;
                    gate_identity = jaccard;  /* approx identity */
                }
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
                    gate_identity = best_identity;
                } else if (best_identity >= REFINE_RELAXED_IDENTITY &&
                           best_coverage >= REFINE_RELAXED_COVERAGE) {
                    if (check_instance_overlap(fi, fj)) {
                        should_merge = 1;
                        gate_identity = best_identity;
                    }
                }
            }

            if (should_merge) {
                /* Nested-element gate (P1 follow-up): if the shorter
                 * family's consensus is fully contained in the longer's
                 * (100%-coverage match) but many of the shorter's
                 * INSTANCES are not nested inside longer's instances,
                 * the two are biologically distinct families that
                 * happen to be sequence-similar (e.g. SINE inside LINE
                 * + SINE solo elsewhere).  Without this check, strict
                 * 80-80-80 destroys the shorter family.
                 *
                 * Trigger only when the shorter is highly contained in
                 * the longer (coverage near 100%) — that's the nested-
                 * candidate signal.  Then require ≥80% of the shorter's
                 * instances be spatially contained; otherwise refuse. */
                const CandidateFamily *shorter_f =
                    (fi->consensus_length <= fj->consensus_length) ? fi : fj;
                const CandidateFamily *longer_f =
                    (fi->consensus_length <= fj->consensus_length) ? fj : fi;
                if (longer_f->consensus_length >=
                        shorter_f->consensus_length * 3 &&
                    shorter_f->num_instances >= 3) {
                    float ctf = nested_containment_fraction(shorter_f, longer_f);
                    if (ctf < 0.50f) {
                        w->mdl_vetoes++;
                        if (w->verbose >= 2)
                            fprintf(stderr, "  Merge nested-veto: F%u + F%u "
                                    "(shorter=%d, longer=%d, "
                                    "containment=%.2f)\n",
                                    fi->id, fj->id,
                                    shorter_f->consensus_length,
                                    longer_f->consensus_length,
                                    ctf);
                        continue;
                    }
                }

                /* MDL gate (symmetric with split / fragment-assembly):
                 * veto if the merged family wouldn't pass MDL on its own.
                 * Without this gate, union-find can chain marginally-similar
                 * pairs into a single super-family with negative MDL, which
                 * mdl_select_library then rejects — losing all originals. */
                if (gate_identity <= 0.0f) gate_identity = REFINE_RELAXED_IDENTITY;
                double merged_mdl = estimate_merge_score(fi, fj, gate_identity,
                                                         w->num_families);
                if (merged_mdl <= 0.0) {
                    w->mdl_vetoes++;
                    if (w->verbose >= 2)
                        fprintf(stderr, "  Merge MDL veto: F%u + F%u "
                                "(identity=%.2f, merged_mdl=%.1f)\n",
                                fi->id, fj->id, gate_identity, merged_mdl);
                    continue;
                }

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

                /* Length-ratio guard (Stage B fix 2; see merge_worker_fn). */
                if ((float)qlen / (float)tlen < 0.7f)
                    continue;

                /* Check if DP would exceed cell limit for long consensus pairs */
                int64_t dp_cells = (int64_t)(qlen + 1) * (int64_t)(tlen + 1);
                if (dp_cells > g_refine_max_dp_cells) {
                    /* Fallback: Jaccard + instance overlap for long families */
                    if (jaccard >= 0.80f && check_instance_overlap(fi, fj)) {
                        /* MDL gate (see merge_worker_fn for rationale) */
                        double mm = estimate_merge_score(fi, fj, jaccard, n);
                        if (mm <= 0.0) {
                            if (verbose >= 2)
                                fprintf(stderr, "  Merge MDL veto (long): "
                                        "F%d + F%d (jaccard=%.2f, "
                                        "merged_mdl=%.1f)\n",
                                        fi->id, fj->id, jaccard, mm);
                        } else {
                            uf_unite(parent, uf_rank, i, j);
                            n_merges++;
                            if (verbose >= 2)
                                fprintf(stderr, "  Merge (long-fallback): F%d (len=%d) +"
                                        " F%d (len=%d) jaccard=%.2f\n",
                                        fi->id, fi->consensus_length,
                                        fj->id, fj->consensus_length, jaccard);
                        }
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
                    /* Nested-element gate (P1 follow-up): see
                     * merge_worker_fn for rationale. */
                    const CandidateFamily *shorter_f =
                        (fi->consensus_length <= fj->consensus_length) ? fi : fj;
                    const CandidateFamily *longer_f =
                        (fi->consensus_length <= fj->consensus_length) ? fj : fi;
                    if (longer_f->consensus_length >=
                            shorter_f->consensus_length * 3 &&
                        shorter_f->num_instances >= 3) {
                        float ctf = nested_containment_fraction(shorter_f, longer_f);
                        if (ctf < 0.50f) {
                            if (verbose >= 2)
                                fprintf(stderr, "  Merge nested-veto: F%d + F%d "
                                        "(shorter=%d, longer=%d, "
                                        "containment=%.2f)\n",
                                        fi->id, fj->id,
                                        shorter_f->consensus_length,
                                        longer_f->consensus_length, ctf);
                            continue;
                        }
                    }

                    /* MDL gate */
                    double mm = estimate_merge_score(fi, fj, best_identity, n);
                    if (mm <= 0.0) {
                        if (verbose >= 2)
                            fprintf(stderr, "  Merge MDL veto: F%d + F%d "
                                    "(identity=%.2f, merged_mdl=%.1f)\n",
                                    fi->id, fj->id, best_identity, mm);
                        continue;
                    }
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
                        /* MDL gate */
                        double mm = estimate_merge_score(fi, fj, best_identity, n);
                        if (mm <= 0.0) {
                            if (verbose >= 2)
                                fprintf(stderr, "  Merge MDL veto (relaxed): "
                                        "F%d + F%d (identity=%.2f, "
                                        "merged_mdl=%.1f)\n",
                                        fi->id, fj->id, best_identity, mm);
                        } else {
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
            workers[t].num_families = n;       /* R estimate for MDL gate */
            workers[t].verbose = verbose;
            workers[t].mdl_vetoes = 0;
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
        int total_vetoes = 0;
        for (int t = 0; t < num_threads; t++) {
            for (int p = 0; p < workers[t].num_pairs; p++) {
                int i = workers[t].pairs[p].i;
                int j = workers[t].pairs[p].j;
                if (uf_find(parent, i) != uf_find(parent, j)) {
                    uf_unite(parent, uf_rank, i, j);
                    n_merges++;
                }
            }
            total_vetoes += workers[t].mdl_vetoes;
            free(workers[t].pairs);
        }
        if (verbose && total_vetoes > 0)
            fprintf(stderr, "Merge MDL gate: vetoed %d of %d candidate pairs\n",
                    total_vetoes, total_vetoes + n_merges);
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

    /* Representative = the family in the component with the most instances,
     * ties broken by LOWEST family index. The previous version initialized
     * representative[root] = root and used a strict > (no replacement on tie),
     * so an instance-count tie kept the union-find root as representative — and
     * that root is set by union-by-rank in the order merge pairs are applied,
     * which is nondeterministic under parallel pair collection. Indexing the
     * tie-break on the root therefore made the chosen representative (hence the
     * merged consensus) vary run-to-run. Iterating i in increasing order with a
     * -1 sentinel makes the lowest-index max-instance family win, deterministically. */
    for (int i = 0; i < n; i++)
        representative[i] = -1;

    for (int i = 0; i < n; i++) {
        int root = uf_find(parent, i);
        int cur = representative[root];
        if (cur < 0 ||
            cl->families[i].num_instances > cl->families[cur].num_instances) {
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
    /* Re-refine each component's representative iff it absorbed at least one
     * other family. The decision is keyed only on `rep` (deterministic: the
     * max-instance family of the component) and component membership — NOT on
     * the union-find tree root. The root is set by union-by-rank in the order
     * merge pairs are applied, which is nondeterministic under parallel pair
     * collection; a former `i == root` guard therefore made *which* merged
     * families got re-refined vary run-to-run, changing their rebuilt consensus
     * (the downstream md5 nondeterminism). */
    for (int i = 0; i < n; i++) {
        int root = uf_find(parent, i);
        int rep = representative[root];
        if (i != rep) continue;

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

    free(representative);
    free(parent);
    free(uf_rank);

    /* --- Step 6: Compact — remove absorbed families --- */
    int write_idx = 0;
    for (int i = 0; i < cl->num_families; i++) {
        if (cl->families[i].num_instances > 0) {
            if (write_idx != i)
                cl->families[write_idx] = cl->families[i];
            cl->families[write_idx].id = (mdl_uid_t)write_idx;
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

/* ----------------------------------------------------------------
 * Parallel split: result struct for Phase 1 (analysis)
 * ---------------------------------------------------------------- */

typedef struct {
    int       family_idx;    /* which family was split */
    int       valid;         /* 1 if split accepted */
    int       audited;
    int       family_id;
    int       consensus_length;
    int       num_instances;
    const char *decision;
    float     threshold;
    float     bimodality;
    int       valley_pass;
    float     div_gap;
    double    min_acceptable;
    Instance *lo_insts;
    int       n_lo;
    char     *lo_cons;
    int       lo_clen;
    double    lo_mdl_score;
    double    lo_model_cost;
    Instance *hi_insts;
    int       n_hi;
    char     *hi_cons;
    int       hi_clen;
    double    hi_mdl_score;
    double    hi_model_cost;
    float     mean_lo;       /* for verbose output */
    float     mean_hi;
    double    orig_score;
    double    split_score;
} SplitResult;

static FILE *split_audit_open(void)
{
    if (!g_refine_split_audit_path) return NULL;

    FILE *fp = fopen(g_refine_split_audit_path, "w");
    if (!fp) {
        fprintf(stderr, "WARNING: could not open split audit '%s'\n",
                g_refine_split_audit_path);
        return NULL;
    }

    fprintf(fp, "family_id\tconsensus_length\tnum_instances\tdecision\t"
                "threshold\tbimodality\tvalley_pass\tn_lo\tn_hi\t"
                "mean_lo\tmean_hi\tdiv_gap\torig_score\tsplit_score\t"
                "min_acceptable\tlo_score\thi_score\tlo_len\thi_len\n");
    return fp;
}

static void split_audit_write(FILE *fp, int family_id, int consensus_length,
                              int num_instances, const char *decision,
                              float threshold, float bimodality,
                              int valley_pass, int n_lo, int n_hi,
                              float mean_lo, float mean_hi,
                              double orig_score, double split_score,
                              double min_acceptable,
                              double lo_score, double hi_score,
                              int lo_len, int hi_len)
{
    if (!fp) return;
    float div_gap = (n_lo > 0 && n_hi > 0) ? (mean_hi - mean_lo) : 0.0f;
    fprintf(fp, "%d\t%d\t%d\t%s\t%.6f\t%.6f\t%d\t%d\t%d\t"
                "%.6f\t%.6f\t%.6f\t%.3f\t%.3f\t%.3f\t"
                "%.3f\t%.3f\t%d\t%d\n",
            family_id, consensus_length, num_instances,
            decision ? decision : "unknown",
            threshold, bimodality, valley_pass, n_lo, n_hi,
            mean_lo, mean_hi, div_gap, orig_score, split_score,
            min_acceptable, lo_score, hi_score, lo_len, hi_len);
}

static void split_audit_write_result(FILE *fp, const SplitResult *res)
{
    if (!fp || !res || !res->audited) return;
    split_audit_write(fp, res->family_id, res->consensus_length,
                      res->num_instances, res->decision,
                      res->threshold, res->bimodality, res->valley_pass,
                      res->n_lo, res->n_hi, res->mean_lo, res->mean_hi,
                      res->orig_score, res->split_score,
                      res->min_acceptable, res->lo_mdl_score,
                      res->hi_mdl_score, res->lo_clen, res->hi_clen);
}

/* ----------------------------------------------------------------
 * Phase 1 worker: parallel split analysis (read-only on cl)
 * ---------------------------------------------------------------- */

typedef struct {
    const CandidateFamily *families; /* snapshot pointer (read-only) */
    const Genome          *genome;
    glen_t                 genome_len;
    int                    num_families;
    int                    verbose;
    int                    orig_n;
    int                   *next_idx;   /* shared atomic counter */
    SplitResult           *results;    /* pre-allocated array [orig_n] */
} SplitAnalysisArgs;

static void *split_analysis_worker(void *arg)
{
    SplitAnalysisArgs *w = (SplitAnalysisArgs *)arg;

    while (1) {
        int fi = __atomic_fetch_add(w->next_idx, 1, __ATOMIC_RELAXED);
        if (fi >= w->orig_n) break;

        SplitResult *res = &w->results[fi];
        res->family_idx = fi;
        res->valid = 0;

        const CandidateFamily *fam = &w->families[fi];
        res->audited = 1;
        res->family_id = fam->id;
        res->consensus_length = fam->consensus_length;
        res->num_instances = fam->num_instances;
        res->decision = "not_evaluated";
        res->valley_pass = 1;

        if (fam->num_instances < REFINE_MIN_SPLIT_INSTANCES) {
            res->decision = "skip_too_few_instances";
            if (w->verbose >= 2)
                fprintf(stderr, "  [split] F%d: skipped: n_instances=%d < "
                        "REFINE_MIN_SPLIT_INSTANCES=%d\n",
                        fam->id, fam->num_instances, REFINE_MIN_SPLIT_INSTANCES);
            continue;
        }

        /* Collect divergence values */
        float *divs = malloc((size_t)fam->num_instances * sizeof(float));
        if (!divs) {
            res->decision = "error_alloc_divs";
            continue;
        }
        for (int j = 0; j < fam->num_instances; j++)
            divs[j] = fam->instances[j].divergence;

        /* Otsu's method to find split threshold */
        float bimodality;
        int valley_pass;
        float threshold = otsu_threshold(divs, fam->num_instances, &bimodality,
                                         &valley_pass);
        res->threshold = threshold;
        res->bimodality = bimodality;
        res->valley_pass = valley_pass;

        if (bimodality < REFINE_BIMODALITY_THRESH) {
            res->decision = "reject_bimodality";
            if (w->verbose >= 2)
                fprintf(stderr, "  [split] F%d: rejected: bimodality=%.3f < "
                        "threshold=%.3f (n=%d)\n",
                        fam->id, bimodality, REFINE_BIMODALITY_THRESH,
                        fam->num_instances);
            free(divs);
            continue;
        }

        /* Valley depth check */
        if (!valley_pass) {
            res->decision = "reject_valley";
            if (w->verbose >= 2)
                fprintf(stderr, "  [split] F%d: rejected: valley not deep enough "
                        "(bimodality=%.3f, n=%d)\n",
                        fam->id, bimodality, fam->num_instances);
            free(divs);
            continue;
        }

        /* Split instances into two groups */
        int n_lo = 0, n_hi = 0;
        for (int j = 0; j < fam->num_instances; j++) {
            if (divs[j] <= threshold) n_lo++;
            else n_hi++;
        }
        res->n_lo = n_lo;
        res->n_hi = n_hi;

        if (n_lo < REFINE_MIN_CLUSTER_SIZE ||
            n_hi < REFINE_MIN_CLUSTER_SIZE) {
            res->decision = "reject_cluster_size";
            if (w->verbose >= 2)
                fprintf(stderr, "  [split] F%d: rejected: n_lo=%d or n_hi=%d < "
                        "REFINE_MIN_CLUSTER_SIZE=%d (threshold=%.3f, n=%d)\n",
                        fam->id, n_lo, n_hi, REFINE_MIN_CLUSTER_SIZE,
                        threshold, fam->num_instances);
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
        res->mean_lo = mean_lo;
        res->mean_hi = mean_hi;
        res->div_gap = mean_hi - mean_lo;

        if (mean_hi - mean_lo < REFINE_MIN_DIV_GAP) {
            res->decision = "reject_div_gap";
            if (w->verbose >= 2)
                fprintf(stderr, "  [split] F%d: rejected: div_gap=%.3f < "
                        "REFINE_MIN_DIV_GAP=%.3f (mean_lo=%.3f mean_hi=%.3f, n=%d)\n",
                        fam->id, mean_hi - mean_lo, REFINE_MIN_DIV_GAP,
                        mean_lo, mean_hi, fam->num_instances);
            free(divs);
            continue;
        }

        /* Allocate instance arrays for each group */
        Instance *lo_insts = malloc((size_t)n_lo * sizeof(Instance));
        Instance *hi_insts = malloc((size_t)n_hi * sizeof(Instance));
        if (!lo_insts || !hi_insts) {
            res->decision = "error_alloc_split_instances";
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
        int lo_clen = build_subset_consensus(fam, lo_insts, n_lo, w->genome,
                                             &lo_cons);
        int hi_clen = build_subset_consensus(fam, hi_insts, n_hi, w->genome,
                                             &hi_cons);

        if (lo_clen <= 0 || hi_clen <= 0) {
            res->decision = "error_subset_consensus";
            free(lo_insts); free(hi_insts);
            free(lo_cons); free(hi_cons);
            continue;
        }

        /* Score original family (into a local copy to avoid racing) */
        CandidateFamily orig_copy;
        memcpy(&orig_copy, fam, sizeof(CandidateFamily));
        mdl_score_family(&orig_copy, w->genome_len, w->num_families);
        double orig_score = orig_copy.mdl_score;
        res->orig_score = orig_score;

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

        mdl_score_family(&sub_lo, w->genome_len, w->num_families + 1);
        mdl_score_family(&sub_hi, w->genome_len, w->num_families + 1);
        double split_score = sub_lo.mdl_score + sub_hi.mdl_score;
        res->lo_mdl_score = sub_lo.mdl_score;
        res->hi_mdl_score = sub_hi.mdl_score;
        res->lo_clen = lo_clen;
        res->hi_clen = hi_clen;
        res->split_score = split_score;

        /* Relaxed gate (chr4 90×80 experiment 2026-05-02): for positive
         * original families, accept non-negative combined split MDL.  The
         * original strict gate required split>=orig and rejected all attempted
         * chr4 splits. */
        double min_acceptable = (orig_score > 0) ? 0.0 : orig_score - 2.0 * fabs(orig_score);
        res->min_acceptable = min_acceptable;
        if (split_score < min_acceptable) {
            res->decision = "reject_mdl_gate";
            if (w->verbose >= 2)
                fprintf(stderr, "  [split] F%d: rejected: split MDL < min %.1f "
                        "(orig=%.1f, split=%.1f=%.1f+%.1f, n=%d)\n",
                        fam->id, min_acceptable, orig_score, split_score,
                        sub_lo.mdl_score, sub_hi.mdl_score,
                        fam->num_instances);
            free(lo_cons); free(hi_cons);
            free(lo_insts); free(hi_insts);
            continue;
        }

        if (w->verbose >= 2)
            fprintf(stderr, "  [split] F%d: accepted: split into n_lo=%d / n_hi=%d "
                    "at threshold=%.3f (div %.1f%% / %.1f%%), MDL %.1f -> %.1f\n",
                    fam->id, n_lo, n_hi, threshold,
                    mean_lo * 100, mean_hi * 100,
                    orig_score, split_score);

        /* Store accepted result */
        res->valid = 1;
        res->decision = "accept";
        res->lo_insts = lo_insts;
        res->n_lo = n_lo;
        res->lo_cons = lo_cons;
        res->lo_clen = lo_clen;
        res->lo_mdl_score = sub_lo.mdl_score;
        res->lo_model_cost = sub_lo.model_cost;
        res->hi_insts = hi_insts;
        res->n_hi = n_hi;
        res->hi_cons = hi_cons;
        res->hi_clen = hi_clen;
        res->hi_mdl_score = sub_hi.mdl_score;
        res->hi_model_cost = sub_hi.model_cost;
        res->mean_lo = mean_lo;
        res->mean_hi = mean_hi;
        res->orig_score = orig_score;
        res->split_score = split_score;
    }

    return NULL;
}

/* ----------------------------------------------------------------
 * Phase 2 worker: parallel align_refine on split families
 * ---------------------------------------------------------------- */

typedef struct {
    CandidateList     *cl;
    const Genome      *genome;
    const KmerTable   *kt;
    int                k;
    int               *refine_indices;  /* family indices to refine */
    int                n_refine;
    int               *next_refine;     /* shared atomic counter */
} RefineWorkerArgs;

static void *refine_split_worker(void *arg)
{
    RefineWorkerArgs *w = (RefineWorkerArgs *)arg;

    while (1) {
        int ri = __atomic_fetch_add(w->next_refine, 1, __ATOMIC_RELAXED);
        if (ri >= w->n_refine) break;

        int fi = w->refine_indices[ri];
        align_refine_family(&w->cl->families[fi], w->genome, w->kt, w->k,
                            ALIGN_MAX_ITERATIONS);
    }

    return NULL;
}

/* ================================================================
 * refine_split_families — two-phase parallel implementation
 * ================================================================ */

int refine_split_families(CandidateList *cl, const Genome *genome,
                          const KmerTable *kt, int k,
                          glen_t genome_len, int verbose,
                          int num_families, int num_threads)
{
    int orig_n = cl->num_families;

    /* --- Sequential fallback (num_threads <= 1) --- */
    if (num_threads <= 1) {
        int n_splits = 0;
        FILE *audit_fp = split_audit_open();

        for (int fi = 0; fi < orig_n; fi++) {
            CandidateFamily *fam = &cl->families[fi];
            SplitResult audit;
            memset(&audit, 0, sizeof(audit));
            audit.family_idx = fi;
            audit.audited = 1;
            audit.family_id = fam->id;
            audit.consensus_length = fam->consensus_length;
            audit.num_instances = fam->num_instances;
            audit.decision = "not_evaluated";
            audit.valley_pass = 1;

            if (fam->num_instances < REFINE_MIN_SPLIT_INSTANCES) {
                audit.decision = "skip_too_few_instances";
                split_audit_write_result(audit_fp, &audit);
                if (verbose >= 2)
                    fprintf(stderr, "  [split] F%d: skipped: n_instances=%d < "
                            "REFINE_MIN_SPLIT_INSTANCES=%d\n",
                            fam->id, fam->num_instances, REFINE_MIN_SPLIT_INSTANCES);
                continue;
            }

            /* Collect divergence values */
            float *divs = malloc((size_t)fam->num_instances * sizeof(float));
            if (!divs) {
                audit.decision = "error_alloc_divs";
                split_audit_write_result(audit_fp, &audit);
                continue;
            }
            for (int j = 0; j < fam->num_instances; j++)
                divs[j] = fam->instances[j].divergence;

            /* Otsu's method to find split threshold */
            float bimodality;
            int valley_pass;
            float threshold = otsu_threshold(divs, fam->num_instances,
                                             &bimodality, &valley_pass);
            audit.threshold = threshold;
            audit.bimodality = bimodality;
            audit.valley_pass = valley_pass;

            if (bimodality < REFINE_BIMODALITY_THRESH) {
                audit.decision = "reject_bimodality";
                split_audit_write_result(audit_fp, &audit);
                if (verbose >= 2)
                    fprintf(stderr, "  [split] F%d: rejected: bimodality=%.3f < "
                            "threshold=%.3f (n=%d)\n",
                            fam->id, bimodality, REFINE_BIMODALITY_THRESH,
                            fam->num_instances);
                free(divs);
                continue;
            }

            /* Valley depth check */
            if (!valley_pass) {
                audit.decision = "reject_valley";
                split_audit_write_result(audit_fp, &audit);
                if (verbose >= 2)
                    fprintf(stderr, "  [split] F%d: rejected: valley not deep enough "
                            "(bimodality=%.3f, n=%d)\n",
                            fam->id, bimodality, fam->num_instances);
                free(divs);
                continue;
            }

            /* Split instances into two groups */
            int n_lo = 0, n_hi = 0;
            for (int j = 0; j < fam->num_instances; j++) {
                if (divs[j] <= threshold) n_lo++;
                else n_hi++;
            }
            audit.n_lo = n_lo;
            audit.n_hi = n_hi;

            if (n_lo < REFINE_MIN_CLUSTER_SIZE ||
                n_hi < REFINE_MIN_CLUSTER_SIZE) {
                audit.decision = "reject_cluster_size";
                split_audit_write_result(audit_fp, &audit);
                if (verbose >= 2)
                    fprintf(stderr, "  [split] F%d: rejected: n_lo=%d or n_hi=%d < "
                            "REFINE_MIN_CLUSTER_SIZE=%d (threshold=%.3f, n=%d)\n",
                            fam->id, n_lo, n_hi, REFINE_MIN_CLUSTER_SIZE,
                            threshold, fam->num_instances);
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
            audit.mean_lo = mean_lo;
            audit.mean_hi = mean_hi;
            audit.div_gap = mean_hi - mean_lo;

            if (mean_hi - mean_lo < REFINE_MIN_DIV_GAP) {
                audit.decision = "reject_div_gap";
                split_audit_write_result(audit_fp, &audit);
                if (verbose >= 2)
                    fprintf(stderr, "  [split] F%d: rejected: div_gap=%.3f < "
                            "REFINE_MIN_DIV_GAP=%.3f (mean_lo=%.3f mean_hi=%.3f, n=%d)\n",
                            fam->id, mean_hi - mean_lo, REFINE_MIN_DIV_GAP,
                            mean_lo, mean_hi, fam->num_instances);
                free(divs);
                continue;
            }

            /* Allocate instance arrays for each group */
            Instance *lo_insts = malloc((size_t)n_lo * sizeof(Instance));
            Instance *hi_insts = malloc((size_t)n_hi * sizeof(Instance));
            if (!lo_insts || !hi_insts) {
                audit.decision = "error_alloc_split_instances";
                split_audit_write_result(audit_fp, &audit);
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
                audit.decision = "error_subset_consensus";
                audit.lo_clen = lo_clen;
                audit.hi_clen = hi_clen;
                split_audit_write_result(audit_fp, &audit);
                free(lo_insts); free(hi_insts);
                free(lo_cons); free(hi_cons);
                continue;
            }

            /* Score original family */
            mdl_score_family(fam, genome_len, num_families);
            double orig_score = fam->mdl_score;
            audit.orig_score = orig_score;

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

            mdl_score_family(&sub_lo, genome_len, num_families + 1);
            mdl_score_family(&sub_hi, genome_len, num_families + 1);
            double split_score = sub_lo.mdl_score + sub_hi.mdl_score;
            audit.lo_mdl_score = sub_lo.mdl_score;
            audit.hi_mdl_score = sub_hi.mdl_score;
            audit.lo_clen = lo_clen;
            audit.hi_clen = hi_clen;
            audit.split_score = split_score;

            /* Relaxed gate (chr4 90×80 experiment 2026-05-02): for positive
             * original families, accept non-negative combined split MDL. */
            double min_acceptable = (orig_score > 0) ? 0.0 : orig_score - 2.0 * fabs(orig_score);
            audit.min_acceptable = min_acceptable;
            if (split_score < min_acceptable) {
                audit.decision = "reject_mdl_gate";
                split_audit_write_result(audit_fp, &audit);
                if (verbose >= 2)
                    fprintf(stderr, "  [split] F%d: rejected: split MDL < min %.1f "
                            "(orig=%.1f, split=%.1f=%.1f+%.1f, n=%d)\n",
                            fam->id, min_acceptable, orig_score, split_score,
                            sub_lo.mdl_score, sub_hi.mdl_score,
                            fam->num_instances);
                free(lo_cons); free(hi_cons);
                free(lo_insts); free(hi_insts);
                continue;
            }

            /* Accept split */
            audit.decision = "accept";
            audit.valid = 1;
            split_audit_write_result(audit_fp, &audit);
            if (verbose >= 2)
                fprintf(stderr, "  [split] F%d: accepted: split into n_lo=%d / n_hi=%d "
                        "at threshold=%.3f (div %.1f%% / %.1f%%), MDL %.1f -> %.1f\n",
                        fam->id, n_lo, n_hi, threshold,
                        mean_lo * 100, mean_hi * 100,
                        orig_score, split_score);
            else if (verbose)
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
                fam = &cl->families[fi];
            }

            CandidateFamily *new_fam = &cl->families[cl->num_families];
            memset(new_fam, 0, sizeof(CandidateFamily));
            new_fam->id = (mdl_uid_t)cl->num_families;
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
            align_refine_family(&cl->families[cl->num_families - 1], genome,
                                kt, k, ALIGN_MAX_ITERATIONS);

            n_splits++;
        }

        if (audit_fp) fclose(audit_fp);
        return n_splits;
    }

    /* ================================================================
     * Parallel path (num_threads > 1): two-phase approach
     * ================================================================ */

    /* --- Phase 1: parallel split analysis (read-only on cl) --- */

    SplitResult *results = calloc((size_t)orig_n, sizeof(SplitResult));
    if (!results) return 0;

    int next_idx = 0;
    int n_workers = num_threads < orig_n ? num_threads : orig_n;
    if (n_workers < 1) n_workers = 1;

    SplitAnalysisArgs sa_args;
    sa_args.families     = cl->families;
    sa_args.genome       = genome;
    sa_args.genome_len   = genome_len;
    sa_args.num_families = num_families;
    sa_args.verbose      = verbose;
    sa_args.orig_n       = orig_n;
    sa_args.next_idx     = &next_idx;
    sa_args.results      = results;

    if (n_workers > 1) {
        pthread_t *threads = malloc((size_t)(n_workers - 1) * sizeof(pthread_t));
        if (!threads) { free(results); return 0; }

        for (int t = 0; t < n_workers - 1; t++)
            pthread_create(&threads[t], NULL, split_analysis_worker, &sa_args);

        /* Main thread also participates */
        split_analysis_worker(&sa_args);

        for (int t = 0; t < n_workers - 1; t++)
            pthread_join(threads[t], NULL);

        free(threads);
    } else {
        split_analysis_worker(&sa_args);
    }

    {
        FILE *audit_fp = split_audit_open();
        if (audit_fp) {
            for (int fi = 0; fi < orig_n; fi++)
                split_audit_write_result(audit_fp, &results[fi]);
            fclose(audit_fp);
        }
    }

    /* --- Phase 2a: sequential apply of accepted splits --- */

    int n_splits = 0;

    /* Collect indices of families to refine (both halves of each split) */
    int refine_cap = orig_n;
    int *refine_indices = malloc((size_t)refine_cap * sizeof(int));
    if (!refine_indices) { free(results); return 0; }
    int n_refine = 0;

    for (int fi = 0; fi < orig_n; fi++) {
        SplitResult *res = &results[fi];
        if (!res->valid) continue;

        CandidateFamily *fam = &cl->families[fi];

        if (verbose)
            fprintf(stderr, "  Split F%d: %d instances (div %.1f%%) + "
                    "%d instances (div %.1f%%), MDL %.1f -> %.1f\n",
                    fam->id, res->n_lo, res->mean_lo * 100,
                    res->n_hi, res->mean_hi * 100,
                    res->orig_score, res->split_score);

        /* Replace original family with low-divergence group */
        free(fam->consensus);
        free(fam->instances);
        fam->consensus = res->lo_cons;
        fam->consensus_length = res->lo_clen;
        fam->instances = res->lo_insts;
        fam->num_instances = res->n_lo;
        fam->cap_instances = res->n_lo;
        fam->mdl_score = res->lo_mdl_score;
        fam->model_cost = res->lo_model_cost;

        /* Add high-divergence group as new family */
        if (cl->num_families >= cl->cap_families) {
            int new_cap = cl->cap_families * 2;
            CandidateFamily *tmp = realloc(cl->families,
                                           (size_t)new_cap * sizeof(CandidateFamily));
            if (!tmp) {
                /* Cannot append; free hi resources */
                free(res->hi_cons);
                free(res->hi_insts);
                /* lo already applied, still counts as partial split */
                n_splits++;
                /* Still refine the lo half */
                if (n_refine >= refine_cap) {
                    refine_cap *= 2;
                    refine_indices = realloc(refine_indices,
                                             (size_t)refine_cap * sizeof(int));
                }
                if (refine_indices) refine_indices[n_refine++] = fi;
                continue;
            }
            cl->families = tmp;
            cl->cap_families = new_cap;
            fam = &cl->families[fi]; /* re-acquire after realloc */
        }

        int new_idx = cl->num_families;
        CandidateFamily *new_fam = &cl->families[new_idx];
        memset(new_fam, 0, sizeof(CandidateFamily));
        new_fam->id = (mdl_uid_t)new_idx;
        new_fam->consensus = res->hi_cons;
        new_fam->consensus_length = res->hi_clen;
        new_fam->instances = res->hi_insts;
        new_fam->num_instances = res->n_hi;
        new_fam->cap_instances = res->n_hi;
        new_fam->component_id = fam->component_id;
        new_fam->topology = fam->topology;
        new_fam->estimated_copies = (freq_t)res->n_hi;
        new_fam->mdl_score = res->hi_mdl_score;
        new_fam->model_cost = res->hi_model_cost;
        cl->num_families++;

        /* Grow refine_indices if needed */
        if (n_refine + 2 > refine_cap) {
            refine_cap = refine_cap * 2 + 2;
            int *tmp2 = realloc(refine_indices,
                                (size_t)refine_cap * sizeof(int));
            if (!tmp2) { free(refine_indices); refine_indices = NULL; break; }
            refine_indices = tmp2;
        }

        refine_indices[n_refine++] = fi;       /* lo half */
        refine_indices[n_refine++] = new_idx;  /* hi half */
        n_splits++;
    }

    free(results);

    /* --- Phase 2b: parallel align_refine on all split families --- */

    if (n_refine > 0 && refine_indices) {
        int next_refine = 0;
        int rw = num_threads < n_refine ? num_threads : n_refine;
        if (rw < 1) rw = 1;

        RefineWorkerArgs rw_args;
        rw_args.cl             = cl;
        rw_args.genome         = genome;
        rw_args.kt             = kt;
        rw_args.k              = k;
        rw_args.refine_indices = refine_indices;
        rw_args.n_refine       = n_refine;
        rw_args.next_refine    = &next_refine;

        if (rw > 1) {
            pthread_t *threads = malloc((size_t)(rw - 1) * sizeof(pthread_t));
            if (threads) {
                for (int t = 0; t < rw - 1; t++)
                    pthread_create(&threads[t], NULL, refine_split_worker,
                                   &rw_args);

                /* Main thread participates */
                refine_split_worker(&rw_args);

                for (int t = 0; t < rw - 1; t++)
                    pthread_join(threads[t], NULL);

                free(threads);
            } else {
                /* Fallback: sequential refine */
                refine_split_worker(&rw_args);
            }
        } else {
            refine_split_worker(&rw_args);
        }
    }

    free(refine_indices);

    return n_splits;
}

/* ================================================================
 * Phase 6.4: Pruning marginal families after MDL selection
 * ================================================================ */

/* (start, end, family_idx) interval entry used by the prune sweep-line.
 * end is EXCLUSIVE (start + aligned_length).  family_idx points back
 * into cl->families so the sweep can skip self-coverage and treat
 * dynamically pruned families as no longer covering. */
typedef struct {
    gpos_t start;
    gpos_t end;
    int    fam_idx;
} PruneInterval;

static int cmp_prune_interval(const void *x, const void *y)
{
    gpos_t sa = ((const PruneInterval *)x)->start;
    gpos_t sb = ((const PruneInterval *)y)->start;
    if (sa < sb) return -1;
    if (sa > sb) return  1;
    return 0;
}

static int cmp_idxscore_score_asc(const void *x, const void *y);
typedef struct { int idx; double score; } PruneIdxScore;
static int cmp_idxscore_score_asc(const void *x, const void *y)
{
    /* Sort weakest (lowest score) first.  Stable secondary key on idx
     * keeps the order deterministic across equal-score families. */
    double sa = ((const PruneIdxScore *)x)->score;
    double sb = ((const PruneIdxScore *)y)->score;
    if (sa < sb) return -1;
    if (sa > sb) return  1;
    int ia = ((const PruneIdxScore *)x)->idx;
    int ib = ((const PruneIdxScore *)y)->idx;
    return (ia > ib) - (ia < ib);
}

/* Lower-bound search: returns smallest index i such that
 * intervals[i].start >= key, or n if no such i exists. */
static int prune_lower_bound_start(const PruneInterval *iv, int n, gpos_t key)
{
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = lo + ((hi - lo) >> 1);
        if (iv[mid].start < key) lo = mid + 1;
        else                     hi = mid;
    }
    return lo;
}

/*
 * Count bases of [q_start, q_end) covered by ANY interval whose
 * fam_idx is NOT excl_fam AND whose family has not been pruned
 * (pruned[fam_idx] == 0).  Implemented as a bounded sweep over
 * sorted intervals: scan from the first interval whose start could
 * still overlap q_start (binary searched via max_end_before lookup
 * is approximated by walking back a fixed window — simpler approach:
 * scan from first start > q_start - max_len, where max_len is the
 * largest interval end-start seen).  We instead binary-search to
 * the first interval with start >= q_start, then walk backwards
 * while previous starts are still within max_len of q_start.
 *
 * Returns the cardinality of the union of (other-family) overlap
 * with [q_start, q_end), in bp.
 */
static int64_t prune_count_other_coverage(const PruneInterval *iv, int n,
                                          gpos_t q_start, gpos_t q_end,
                                          int excl_fam,
                                          const uint8_t *pruned,
                                          gpos_t max_interval_len)
{
    if (q_end <= q_start || n == 0) return 0;

    /* Find first index whose start could overlap [q_start, q_end).
     * An interval [s,e) overlaps [q_start, q_end) iff s < q_end and
     * e > q_start.  Earliest possible overlapping index has
     * s >= q_start - max_interval_len (since e > q_start ⇒ s > q_start - len).
     * Binary-search for the first such start. */
    gpos_t scan_from = q_start - max_interval_len;
    int i = prune_lower_bound_start(iv, n, scan_from);

    /* Walk forward, sweep-merging overlapping intervals clipped to
     * [q_start, q_end).  Track current merged interval [cur_s, cur_e)
     * and accumulate union length. */
    int64_t covered = 0;
    gpos_t cur_s = -1, cur_e = -1; /* -1 sentinel = no current run */

    for (; i < n; i++) {
        if (iv[i].start >= q_end) break;          /* sorted, all later miss */
        if (iv[i].fam_idx == excl_fam) continue;  /* skip self */
        if (pruned[iv[i].fam_idx]) continue;      /* skip pruned families */
        gpos_t s = iv[i].start;
        gpos_t e = iv[i].end;
        if (e <= q_start) continue;               /* entirely before query */
        if (s < q_start) s = q_start;
        if (e > q_end)   e = q_end;
        if (e <= s)      continue;

        if (cur_s < 0) {
            cur_s = s; cur_e = e;
        } else if (s <= cur_e) {
            if (e > cur_e) cur_e = e;             /* extend */
        } else {
            covered += (int64_t)(cur_e - cur_s);
            cur_s = s; cur_e = e;
        }
    }
    if (cur_s >= 0) covered += (int64_t)(cur_e - cur_s);
    return covered;
}

/*
 * Public entry point for the sweep-line core.  Exposed (file-local)
 * for the sweep-line unit test.  Operates on a CandidateList;
 * computes per-family "exclusive bp per instance" using sorted
 * intervals + binary search rather than an O(genome_len) bitmap.
 *
 * The function returns the number of pruned families (and marks their
 * acceptance state as pruned, while zeroing the legacy mdl_score alias).
 * Memory footprint is O(num_intervals),
 * NOT O(genome_len) — so it scales to multi-Gbp genomes such as
 * wheat (17 Gb) without a single fatal allocation.
 */
int refine_prune_families_sweepline(CandidateList *cl, glen_t genome_len,
                                    int verbose, int num_families)
{
    /* Count accepted families */
    int n_accepted = 0;
    for (int i = 0; i < cl->num_families; i++)
        if (candidate_is_accepted(&cl->families[i])) n_accepted++;

    if (n_accepted <= 1) return 0;

    /* Build (idx, score) pairs for accepted families and sort weakest first.
     * ENG-N10: replaces the previous O(n^2) selection sort with qsort
     * (O(n log n)) — same numerical result, vastly cheaper at n_accepted
     * in the tens of thousands. */
    PruneIdxScore *order = malloc((size_t)n_accepted * sizeof(PruneIdxScore));
    if (!order) return 0;

    int oi = 0;
    for (int i = 0; i < cl->num_families; i++) {
        if (candidate_is_accepted(&cl->families[i])) {
            order[oi].idx = i;
            order[oi].score = candidate_report_score(&cl->families[i]);
            oi++;
        }
    }
    qsort(order, (size_t)n_accepted, sizeof(PruneIdxScore),
          cmp_idxscore_score_asc);

    /* ENG-N2: replace `uint8_t cov[genome_len]` (1 byte / base; 17 GB on
     * wheat) with a sorted array of intervals.  Memory is O(num_intervals)
     * which is bounded by total accepted-family instances, NOT genome
     * size.  The "exclusive coverage" query (how many bases of an
     * instance are covered ONLY by the candidate family) is implemented
     * by counting OTHER-family coverage of the instance via a binary-
     * searched local sweep, then subtracting from instance length.
     *
     * Pruning a family does NOT physically remove its intervals; we
     * just flip pruned[fam_idx] = 1 and the sweep skips its contribution
     * for subsequent candidates. */
    int64_t total_intervals = 0;
    for (int i = 0; i < cl->num_families; i++) {
        if (!candidate_is_accepted(&cl->families[i])) continue;
        total_intervals += cl->families[i].num_instances;
    }

    if (total_intervals == 0) { free(order); return 0; }

    PruneInterval *intervals =
        malloc((size_t)total_intervals * sizeof(PruneInterval));
    if (!intervals) { free(order); return 0; }

    /* pruned flag is per-family-in-the-CandidateList, indexed by idx
     * (NOT by position in order[]).  Allocated per-family, not per-base,
     * so it's O(cl->num_families) bytes — negligible. */
    uint8_t *pruned = calloc((size_t)cl->num_families, 1);
    if (!pruned) { free(intervals); free(order); return 0; }

    int64_t ii = 0;
    gpos_t  max_len = 0;
    for (int i = 0; i < cl->num_families; i++) {
        if (!candidate_is_accepted(&cl->families[i])) continue;
        CandidateFamily *f = &cl->families[i];
        for (int j = 0; j < f->num_instances; j++) {
            gpos_t s = f->instances[j].position;
            gpos_t e = s + (gpos_t)f->instances[j].aligned_length;
            /* Clip to genome bounds (matches original behavior, which
             * silently ignored out-of-range positions). */
            if (s < 0) s = 0;
            if (e > genome_len) e = genome_len;
            if (e <= s) continue;
            intervals[ii].start   = s;
            intervals[ii].end     = e;
            intervals[ii].fam_idx = i;
            if (e - s > max_len) max_len = e - s;
            ii++;
        }
    }
    int64_t n_intervals = ii;
    qsort(intervals, (size_t)n_intervals, sizeof(PruneInterval),
          cmp_prune_interval);

    int n_pruned = 0;

    for (int k = 0; k < n_accepted; k++) {
        int fi = order[k].idx;
        CandidateFamily *f = &cl->families[fi];
        if (!candidate_is_accepted(f)) continue; /* defensive */
        if (pruned[fi]) continue;        /* defensive */

        double exclusive_savings = 0;
        int exclusive_instances = 0;

        for (int j = 0; j < f->num_instances; j++) {
            Instance *inst = &f->instances[j];
            gpos_t s = inst->position;
            gpos_t e = s + (gpos_t)inst->aligned_length;
            if (s < 0) s = 0;
            if (e > genome_len) e = genome_len;
            int alen = (int)inst->aligned_length;
            if (e <= s || alen <= 0) continue;

            int64_t other = prune_count_other_coverage(
                intervals, (int)n_intervals, s, e, fi, pruned, max_len);
            int64_t excl  = (int64_t)(e - s) - other;
            if (excl < 0) excl = 0;             /* defensive */
            int excl_bases = (int)excl;

            /* Skip instances with <25% exclusive coverage (matches
             * original threshold). */
            if (excl_bases < alen / 4) continue;

            int edits = (int)(inst->divergence * excl_bases + 0.5f);
            double lit = 2.0 * excl_bases;
            double enc = mdl_instance_cost_full(excl_bases, edits,
                                                f->consensus_length,
                                                num_families);
            exclusive_savings += (lit - enc);
            exclusive_instances++;
        }

        double exclusive_score = exclusive_savings - candidate_model_cost(f);

        if (exclusive_instances == 0) {
            /* No instance has >=25% exclusive coverage — purely redundant.
             * Mark family pruned; subsequent sweep queries will skip its
             * intervals, exactly mirroring the original cov[]-- behavior. */
            if (verbose)
                fprintf(stderr, "  Pruned F%d: excl_score=%.1f "
                        "(model=%.1f, excl_inst=%d)\n",
                        f->id, exclusive_score, candidate_model_cost(f),
                        exclusive_instances);
            pruned[fi] = 1;
            f->mdl_score = 0;
            f->mdl.report_score = 0.0;
            f->mdl.accept_state = CAND_ACCEPT_PRUNED;
            f->mdl.quality_tier = CAND_TIER_REJECT;
            f->mdl.quality_flags |= CAND_QF_PRUNED_REDUNDANT;
            n_pruned++;
        }
    }

    free(pruned);
    free(intervals);
    free(order);
    return n_pruned;
}

int refine_prune_families(CandidateList *cl, glen_t genome_len, int verbose,
                          int num_families)
{
    /* ENG-N2 + ENG-N10 (QUALITY_PROPOSAL_v6 Tier 1.5b): the previous
     * implementation allocated a uint8_t[genome_len] coverage counter
     * (17 GB on wheat) and used an O(n^2) selection sort over n_accepted.
     * Both are replaced by the sweep-line variant in
     * refine_prune_families_sweepline, which is memory-bound by the
     * total number of instance intervals (not genome length) and uses
     * qsort for ordering.  The numerical decision (prune iff every
     * instance has <25% exclusive coverage from currently-accepted
     * families) is preserved exactly. */
    return refine_prune_families_sweepline(cl, genome_len, verbose,
                                           num_families);
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
    int    seq_index;     /* ENG-N9: chromosome / sequence identifier so the
                           * sweep does not pair instances across chr boundaries.
                           * Used as the primary sort key in cmp_instance_entry. */
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
    /* ENG-N9: primary key seq_index, secondary key start (both ascending).
     * The sweep loop below relies on sorted order to exit early via a
     * `break` once entries[j].start - entries[i].start > D.  Sorting by
     * (seq_index, start) makes that exit fire naturally at every chromo-
     * some boundary: as soon as j moves into a higher seq_index,
     * entries[j].start may be small (back to 0) but entries[j].seq_index
     * differs from entries[i].seq_index, so the cooccurrence loop skips
     * those and the natural "start - start > D" break catches the very
     * next same-chr j (or ends the inner loop entirely if i was the last
     * entry on its chromosome). */
    int ka = ((const InstanceEntry *)x)->seq_index;
    int kb = ((const InstanceEntry *)y)->seq_index;
    if (ka < kb) return -1;
    if (ka > kb) return  1;
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

/* Hash a canonicalized (fam_a, fam_b) pair into a power-of-two table. */
static inline unsigned coocc_pair_hash(int a, int b, int ht_size)
{
    return (((unsigned)a * 2654435761u) ^ ((unsigned)b * 40503u))
           & (unsigned)(ht_size - 1);
}

/* Double the co-occurrence hash table and re-insert every existing pair.
 * Returns 0 on success, -1 on allocation failure (caller aborts the sweep).
 * Growing keeps the load factor bounded so probe chains stay O(1).  Without
 * it the fixed table saturates on repeat-dense genomes (distinct family pairs
 * far exceed the slot count): every lookup then scans the whole table and
 * pairs beyond capacity are silently dropped — the root cause of the
 * multi-hour assembly spin. */
static int coocc_ht_grow(int **ht_p, int *ht_size_p,
                         const CooccPair *pairs, int num_pairs)
{
    int new_size = *ht_size_p * 2;
    int *nht = malloc((size_t)new_size * sizeof(int));
    if (!nht) return -1;
    for (int i = 0; i < new_size; i++) nht[i] = -1;
    for (int p = 0; p < num_pairs; p++) {
        unsigned h = coocc_pair_hash(pairs[p].fam_a, pairs[p].fam_b, new_size);
        for (int probe = 0; probe < new_size; probe++) {
            unsigned idx = (h + (unsigned)probe) & (unsigned)(new_size - 1);
            if (nht[idx] == -1) { nht[idx] = p; break; }
        }
    }
    free(*ht_p);
    *ht_p = nht;
    *ht_size_p = new_size;
    return 0;
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
            entries[ei].seq_index = inst->seq_index;  /* ENG-N9 */
            ei++;
        }
    }
    entries_count = ei;

    qsort(entries, (size_t)entries_count, sizeof(InstanceEntry),
          cmp_instance_entry);

    /* ---- Step 2: Compute proximity distance D ---- */
    /* D = median_consensus_length * 4, clamped to [500, 30000].
     * Was [500, 5000] but observed (chr4 RM benchmark): long LTR-like
     * truth intervals (>10 kb) often fragment into 4-6 family pieces
     * spanning 15-30 kb.  A 5 kb proximity cap couldn't see them as
     * co-occurring; bumping to 30 kb lets the synteny detection see
     * longer-range patterns. */
    int64_t cons_len_sum = 0;
    for (int i = 0; i < n; i++)
        cons_len_sum += cl->families[i].consensus_length;
    int median_cons = (int)(cons_len_sum / n);
    int D = median_cons * 4;
    if (D < 500)   D = 500;
    if (D > 30000) D = 30000;

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

    /* Hash table for pair lookup: map (fam_a, fam_b) -> index in pairs[].
     * Open-addressing, power-of-two sized, and GROWN as pairs accumulate so
     * the load factor stays bounded (see coocc_ht_grow for why a fixed table
     * is the multi-hour-spin bug). */
    int ht_size = 16 * 1024;  /* must stay a power of two */
    int *ht = malloc((size_t)ht_size * sizeof(int));
    if (!ht) { free(entries); free(pairs); return 0; }
    for (int t = 0; t < ht_size; t++) ht[t] = -1;

    if (verbose)
        fprintf(stderr, "  Fragment assembly: co-occurrence sweep over %d "
                "instances (D=%d bp)...\n", entries_count, D);

    /* Work budget: the sweep is O(entries x window-density).  On pathological
     * inputs (very many dense instances) the comparison count can explode;
     * assembly is an optional MDL-gated enhancement, so past the budget we
     * skip it with a diagnostic instead of spinning indefinitely. */
    const int64_t SWEEP_BUDGET = 4000000000LL;
    int64_t sweep_work = 0;
    int sweep_aborted = 0;

    for (int i = 0; i < entries_count; i++) {
        int fa = entries[i].family_idx;

        /* Look forward within distance D */
        for (int j = i + 1; j < entries_count; j++) {
            if (entries[j].seq_index != entries[i].seq_index) break;
            if (entries[j].start - entries[i].start > D) break;

            if (++sweep_work > SWEEP_BUDGET) {
                fprintf(stderr,
                        "WARNING: fragment-assembly co-occurrence sweep exceeded "
                        "its work budget (%lld comparisons over %d instances); "
                        "skipping assembly for this run. Lower family/instance "
                        "counts (e.g. -max-instances) if assembly is required.\n",
                        (long long)SWEEP_BUDGET, entries_count);
                sweep_aborted = 1;
                goto sweep_done;
            }

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

            /* Keep the table below ~70% load so probe chains stay short and
             * an empty slot always exists (no silent drops). */
            if ((int64_t)(num_pairs + 1) * 10 >= (int64_t)ht_size * 7) {
                if (coocc_ht_grow(&ht, &ht_size, pairs, num_pairs) != 0) {
                    sweep_aborted = 1;   /* OOM — abandon assembly safely */
                    goto sweep_done;
                }
            }

            /* Find or create pair in hash table */
            unsigned h = coocc_pair_hash(pa, pb, ht_size);
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

sweep_done:
    free(ht);

    if (sweep_aborted) {
        free(entries);
        free(pairs);
        return 0;
    }

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

        /* Co-occurrence FRACTION guard (#3 chimera): the absolute >=3 floor is
         * trivially met by chance for high-copy families, which can fuse
         * unrelated-but-adjacent TEs into a chimeric consensus.  Require a
         * meaningful fraction of the smaller family's copies to co-occur. */
        if (min_inst > 0 &&
            (float)p->count < REFINE_ASSEMBLE_MIN_COOCC_FRAC * (float)min_inst) {
            if (verbose >= 2)
                fprintf(stderr, "  Assembly skip (low co-occ fraction): F%d/F%d "
                        "count=%d / min_inst=%d (%.0f%% < %.0f%%)\n",
                        fa->id, fb->id, p->count, min_inst,
                        100.0f * (float)p->count / (float)min_inst,
                        100.0f * REFINE_ASSEMBLE_MIN_COOCC_FRAC);
            continue;
        }

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
                if (entries[j].seq_index != entries[i].seq_index) break;
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
                if (entries[j].seq_index != entries[i].seq_index) break;
            if (entries[j].start - entries[i].start > D) break;
                if (entries[j].family_idx == fb_idx) a_first_count++;
            }
        }
        for (int i = 0; i < entries_count; i++) {
            if (entries[i].family_idx != fb_idx) continue;
            for (int j = i + 1; j < entries_count; j++) {
                if (entries[j].seq_index != entries[i].seq_index) break;
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

        /* Mean divergence of the merged instances (chimera-fit measure).
         * Each merged instance carries its divergence against its own part
         * consensus; a copy-weighted mean over all merged instances flags
         * joins whose constituent copies fit poorly (chimera signature). */
        double assembled_div_sum = 0.0;
        for (int j = 0; j < merged_n; j++)
            assembled_div_sum += merged_inst[j].divergence;
        double assembled_mean_div =
            (merged_n > 0) ? assembled_div_sum / (double)merged_n : 1.0;

        /* Accept gate: the join must (a) beat the sum of parts in MDL,
         * (b) yield a positive MDL score (so we never waste a slot on a join
         * that selection would drop at mdl_score<=0), and (c) keep merged-copy
         * divergence within the real-TE band — the volume-driven MDL gate alone
         * passes high-copy chimeras with large positive scores. */
        if (assembled.mdl_score <= sum_individual ||
            assembled.mdl_score <= 0.0 ||
            assembled_mean_div > REFINE_ASSEMBLE_MAX_DIV) {
            if (verbose && assembled.mdl_score > sum_individual)
                fprintf(stderr,
                        "  Assembly REJECT F%d+F%d len=%d n=%d MDL=%.1f "
                        "div=%.3f (gate: mdl>0 & div<=%.2f)\n",
                        first_fam->id, second_fam->id, new_cons_len, merged_n,
                        assembled.mdl_score, assembled_mean_div,
                        REFINE_ASSEMBLE_MAX_DIV);
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
        if (verbose)
            fprintf(stderr, "    (accepted div=%.3f)\n", assembled_mean_div);

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
                cl->families[write_idx].id = (mdl_uid_t)write_idx;
                write_idx++;
            }
        }
        cl->num_families = write_idx;
    }

    free(uf_parent);
    free(uf_rnk);

    return n_assemblies;
}

/* ================================================================
 * Phase 6.5: Tandem-instance coalescing (post-selection reporting)
 * ================================================================
 *
 * Walks each accepted family's instance list and merges consecutive
 * same-strand instances whose gap < coalesce_factor × consensus_length.
 * Operates AFTER mdl_select_library + refine_prune_families: this is a
 * pure reporting transform, NOT a model decision.
 *
 * Why: empirical chr4 RM benchmark shows long truth intervals (>10 kb)
 * are tandem arrays of a single underlying repeat unit.  RM merges them
 * into one annotation; mdl-repeat naturally emits each copy as its own
 * instance.  The mismatch tanks interval-level recall even when the bp-
 * level coverage is fine.  Coalescing aligns reporting with truth
 * convention.
 */

typedef struct { gpos_t start; int idx; } InstByStart;
static int cmp_inst_by_start(const void *a, const void *b)
{
    gpos_t sa = ((const InstByStart *)a)->start;
    gpos_t sb = ((const InstByStart *)b)->start;
    if (sa < sb) return -1;
    if (sa > sb) return 1;
    return 0;
}

int refine_coalesce_tandem_instances(CandidateList *cl,
                                     float coalesce_factor,
                                     int verbose)
{
    int total_coalesced = 0;
    int total_families_affected = 0;

    for (int fi = 0; fi < cl->num_families; fi++) {
        CandidateFamily *f = &cl->families[fi];
        if (!candidate_is_accepted(f)) continue;
        if (f->num_instances <= 1)  continue;

        int gap_threshold = (int)(coalesce_factor * (float)f->consensus_length);
        if (gap_threshold < 50) gap_threshold = 50;  /* min sanity */

        /* Sort instances by genome position to enable single-pass scan. */
        InstByStart *order = malloc((size_t)f->num_instances * sizeof(InstByStart));
        if (!order) continue;
        for (int j = 0; j < f->num_instances; j++) {
            order[j].start = f->instances[j].position;
            order[j].idx   = j;
        }
        qsort(order, (size_t)f->num_instances, sizeof(InstByStart),
              cmp_inst_by_start);

        /* Walk sorted instances; coalesce successive pairs into the
         * leftmost.  Mark merged ones with aligned_length = -1 to remove
         * after the pass. */
        int family_coalesced = 0;
        int  active_idx = order[0].idx;
        Instance *active = &f->instances[active_idx];

        for (int k = 1; k < f->num_instances; k++) {
            int    cur_idx = order[k].idx;
            Instance *cur  = &f->instances[cur_idx];

            /* ENG-N8 cross-chromosome guard.  Per R3, the runtime guard
             * alone (no sort-key change) is sufficient: when we cross
             * a chromosome boundary just advance the active pointer.
             * Without this, a tandem-coalesce can fuse instances from
             * two different sequences whose gpos_t happen to be close,
             * producing a BED interval that spans a chromosome boundary. */
            if (cur->seq_index != active->seq_index) {
                active_idx = cur_idx;
                active     = cur;
                continue;
            }

            gpos_t active_end = active->position + active->aligned_length;
            int    gap        = (int)(cur->position - active_end);

            if (cur->strand == active->strand &&
                gap >= -10 &&            /* tiny overlap OK (TSD) */
                gap <= gap_threshold) {
                /* Coalesce cur into active */
                gpos_t new_end = cur->position + cur->aligned_length;
                if (new_end > active_end) {
                    glen_t old_alen = active->aligned_length;
                    active->aligned_length = (glen_t)(new_end - active->position);

                    /* Sum edit counts; keep score additive. */
                    active->num_edits += cur->num_edits;
                    /* Add gap as edit-equivalent (insertions in the
                     * tandem-array spacer are real edits relative to
                     * a tandem-coalesced consensus interpretation). */
                    if (gap > 0) active->num_edits += gap;

                    /* Recompute divergence as edits / new aligned_length. */
                    if (active->aligned_length > 0) {
                        active->divergence =
                            (float)active->num_edits /
                            (float)active->aligned_length;
                        if (active->divergence > 1.0f) active->divergence = 1.0f;
                    }
                    active->cons_end   = active->cons_end > cur->cons_end
                                       ? active->cons_end : cur->cons_end;
                    active->score     += cur->score;
                    (void)old_alen;
                }
                /* Mark cur as removed */
                cur->aligned_length = -1;
                family_coalesced++;
            } else {
                /* Move active forward */
                active_idx = cur_idx;
                active     = cur;
            }
        }

        free(order);

        /* Compact: remove instances with aligned_length == -1. */
        if (family_coalesced > 0) {
            int wi = 0;
            for (int j = 0; j < f->num_instances; j++) {
                if (f->instances[j].aligned_length < 0) continue;
                if (wi != j) f->instances[wi] = f->instances[j];
                wi++;
            }
            f->num_instances = wi;
            total_coalesced += family_coalesced;
            total_families_affected++;
        }
    }

    if (verbose && total_coalesced > 0) {
        fprintf(stderr,
                "Tandem coalescing: merged %d instance pairs across %d families\n",
                total_coalesced, total_families_affected);
    }
    return total_coalesced;
}

/* ================================================================
 * Drop chimeric / over-extended long families (post-selection)
 * ================================================================ */

int refine_drop_chimeric_long(CandidateList *cl, int verbose)
{
    int dropped = 0;

    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *f = &cl->families[i];
        if (!candidate_is_accepted(f)) continue;
        if (f->consensus_length < REFINE_CHIMERA_MIN_LEN) continue;
        if (f->num_instances <= 0) continue;

        /* Mean per-copy divergence against the consensus. For a long chimera or
         * over-extension the copies match only a fragment, so the unmatched
         * bulk scores as edits and mean divergence is high; a genuine long
         * element is conserved and stays low. */
        double sum = 0.0;
        for (int j = 0; j < f->num_instances; j++)
            sum += (double)f->instances[j].divergence;
        double mean_div = sum / (double)f->num_instances;

        if (mean_div > REFINE_CHIMERA_MAX_DIV) {
            f->mdl.accept_state = CAND_ACCEPT_REJECTED;
            dropped++;
            if (verbose)
                fprintf(stderr,
                        "  drop chimeric-long F%d len=%d copies=%d mean_div=%.3f\n",
                        f->id, f->consensus_length, f->num_instances, mean_div);
        }
    }

    if (verbose)
        fprintf(stderr,
                "Chimera filter: dropped %d long high-divergence families "
                "(len>=%d, div>%.2f)\n",
                dropped, REFINE_CHIMERA_MIN_LEN, (double)REFINE_CHIMERA_MAX_DIV);

    return dropped;
}
