#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mdl.h"

/* Global MDL mode — default: exact lgamma binomial */
MdlMode g_mdl_mode = MDL_MODE_EXACT;

/* ================================================================
 * Rissanen's universal integer code
 * ================================================================ */

/* log2(c0) where c0 = 2.865064... is the normalizing constant */
static const double LOG2_C0 = 1.5179605508986484;

double L_int(int64_t n)
{
    /* L_int is only defined for positive integers (Rissanen 1978).
     * Returning 0.0 silently let callers propagate degenerate inputs
     * into MDL scoring as if they were free; INFINITY makes any
     * downstream cost-based comparison reject the family unambiguously,
     * surfacing the upstream bug that produced the bad input. */
    if (n <= 0) {
        fprintf(stderr, "L_int: n must be positive (got %" PRId64 ") — "
                "returning INFINITY (call site has invalid input)\n", n);
        return INFINITY;
    }

    double result = LOG2_C0;
    double val = log2((double)n);

    while (val > 0.0) {
        result += val;
        val = log2(val);
    }

    return result;
}

/* ================================================================
 * Per-family MDL scoring
 * ================================================================ */

/*
 * Cost to encode a sequence of length L as a literal (2 bits/base).
 * Accepts int64_t to avoid overflow on multi-gigabase genomes (> 2.1 Gb).
 */
static inline double literal_cost(int64_t length)
{
    return 2.0 * (double)length;
}

/*
 * Per-instance encoding cost.
 *
 * Retains the original sensitivity-preserving formula:
 *   L_int(a_i) + L_int(m_i + 1) + m_i * log2(3) + [position encoding]
 *
 * DESIGN DECISION (Audit §C2): consensus_length and num_families are
 * accepted for API compatibility but (void)-ignored.  The DESIGN_DOC
 * overhead terms (type bit, ceil(log2(R)) family ID, strand bit,
 * consensus pointer) were found too aggressive on real genomes, causing
 * unacceptable sensitivity loss.  Omitting ceil(log2(R)) introduces a
 * ~3-6% bias toward accepting more families — an acceptable trade-off
 * for de novo discovery where sensitivity > specificity.  This also
 * renders the iterative R convergence loop in mdl_select_library()
 * and the post-prune recovery pass in main.c effectively inert (scores
 * do not change when R changes).  Retained as framework for future
 * refinement if per-instance overhead terms are re-introduced.
 */
double mdl_instance_cost_full(int aligned_length, int edits,
                              int consensus_length, int num_families)
{
    (void)consensus_length;
    (void)num_families;

    int a_i = aligned_length;
    int m_i = edits;

    /* Clamp edits to valid range */
    if (m_i < 0) m_i = 0;
    if (m_i > a_i) m_i = a_i;
    if (a_i < 1) a_i = 1;

    double cost = 0.0;

    /* Cost to specify the instance length */
    cost += L_int(a_i);

    /* Cost to specify the number of edits (+1 to avoid L_int(0)) */
    cost += L_int(m_i + 1);

    /* Cost to specify each edit type (3 alternatives at each position) */
    cost += m_i * log2(3.0);

    /* Cost to specify edit positions: log2(C(a_i, m_i)) */
    switch (g_mdl_mode) {
    case MDL_MODE_NONE:
        break;
    case MDL_MODE_EXACT:
        if (m_i > 0 && m_i < a_i)
            cost += (lgamma(a_i + 1.0) - lgamma(m_i + 1.0)
                     - lgamma(a_i - m_i + 1.0)) / log(2.0);
        break;
    case MDL_MODE_UPPER:
        if (m_i > 0)
            cost += m_i * log2((double)a_i);
        break;
    }

    return cost;
}

/*
 * Legacy wrapper for backward compatibility.
 * Uses R=2 and consensus_length=aligned_length as conservative defaults.
 */
double mdl_instance_cost(int aligned_length, int edits)
{
    return mdl_instance_cost_full(aligned_length, edits, aligned_length, 2);
}

/*
 * Model cost for a family with consensus of length L:
 *   model_cost = L_int(L) + 2 * L
 * (L_int for the length, then 2 bits per base for the consensus itself)
 */
double mdl_model_cost(int consensus_length)
{
    /* Guard: a family with non-positive consensus length is malformed.
     * Returning INFINITY surfaces the bug rather than letting a negative
     * cost flip downstream MDL math. */
    if (consensus_length <= 0) {
        fprintf(stderr, "mdl_model_cost: invalid consensus_length=%d\n",
                consensus_length);
        return INFINITY;
    }
    return L_int(consensus_length) + literal_cost(consensus_length);
}

void mdl_score_family(CandidateFamily *fam, glen_t genome_len, int num_families)
{
    (void)genome_len;

    fam->model_cost = mdl_model_cost(fam->consensus_length);

    double total_savings = 0.0;
    for (int i = 0; i < fam->num_instances; i++) {
        Instance *inst = &fam->instances[i];
        int a_i = (int)inst->aligned_length;
        /* num_edits comes from align_collect_instances (align.c, banded DP
         * with ALIGN_MATCH / ALIGN_MISMATCH / g_align_gap = 1 / -1 / -5).
         * It is NOT the consensus-vs-consensus alignment used by
         * refine.c's merge stage (REFINE_MATCH/MISMATCH/GAP = 2/-3/-2),
         * which produces identity % only and never feeds MDL. */
        int m_i = inst->num_edits;
        if (m_i < 0) {
            fprintf(stderr, "WARNING: instance with negative num_edits=%d, "
                    "falling back to divergence estimate\n", m_i);
            m_i = (int)(inst->divergence * a_i + 0.5f);
        }
        if (m_i > a_i) m_i = a_i;

        double lit = literal_cost(a_i);
        double enc = mdl_instance_cost_full(a_i, m_i,
                                            fam->consensus_length,
                                            num_families);
        double saving = lit - enc;

        total_savings += saving;
    }

    fam->mdl_score = total_savings - fam->model_cost;
}

/* ================================================================
 * Greedy MDL library selection
 * ================================================================ */

/* Comparison for sorting families by MDL score (descending) */
static int cmp_mdl_desc(const void *a, const void *b)
{
    const CandidateFamily *fa = (const CandidateFamily *)a;
    const CandidateFamily *fb = (const CandidateFamily *)b;
    if (fa->mdl_score > fb->mdl_score) return -1;
    if (fa->mdl_score < fb->mdl_score) return 1;
    return 0;
}

/* (start, end) interval — file-scope so qsort comparator below can take it.
 * Used by the ENG-N12 sweep-line replacement of the covered-bit bitmap.
 * end is EXCLUSIVE.  Memory is O(num_intervals), not O(genome_len). */
typedef struct {
    gpos_t start;
    gpos_t end;
} MdlInterval;

/* qsort comparator: sort MdlInterval by start ascending.  External linkage
 * (no `static`) because the embedded `extern` declaration inside
 * mdl_select_library is plain `extern int cmp_mdl_interval_start(...)`. */
int cmp_mdl_interval_start(const void *x, const void *y)
{
    gpos_t sa = ((const MdlInterval *)x)->start;
    gpos_t sb = ((const MdlInterval *)y)->start;
    if (sa < sb) return -1;
    if (sa > sb) return  1;
    return 0;
}

MDLResult mdl_select_library(CandidateList *cl, glen_t genome_len)
{
    MDLResult result;
    memset(&result, 0, sizeof(result));

    /* Single-pass scoring.  The previous design called for an iterative
     * R-convergence loop (per-instance cost depends on log2(R), R depends
     * on which families are accepted); that loop was removed because the
     * current per-instance cost formula does not depend on R (see
     * mdl_instance_cost_full).  R is passed downstream solely as an API
     * parameter for future overhead-term reintroduction. */
    int R_estimate = cl->num_families;
    if (R_estimate < 2) R_estimate = 2;

    for (int i = 0; i < cl->num_families; i++)
        mdl_score_family(&cl->families[i], genome_len, R_estimate);

    /* Sort by MDL score descending */
    qsort(cl->families, (size_t)cl->num_families, sizeof(CandidateFamily),
          cmp_mdl_desc);

    /*
     * Greedy unique-coverage selection with standalone fallback.
     *
     * The MDL two-part code assumes each genome position is encoded by at
     * most ONE family.  Naively summing per-instance savings across families
     * with overlapping instances double-counts those positions and produces
     * total_savings > 2N, driving compression_ratio negative.
     *
     * Algorithm: walk families in mdl_score-desc order.  For each family,
     * count exclusive (not-yet-claimed) positions per instance and re-evaluate
     * its contribution using only the exclusive bases.  Accept the family if
     * EITHER:
     *   (a) exclusive_savings > model_cost (marginal pass — rigorous),  OR
     *   (b) standalone_score > 0 AND consensus_length >= 50
     *       AND num_instances >= 3 (Stage B fix — preserves valid families
     *       whose territory was claimed by a higher-scoring overlapping
     *       family but who are themselves valid TE families in their own
     *       right; see REFINE_TRACE_REPORT bottleneck #1).
     * Mark accepted instances' positions as claimed before considering the
     * next family.
     *
     * On accepted families, fam->mdl_score is rewritten to the exclusive
     * score (when (a) applies) or kept at the standalone score (when only
     * (b) applies), so downstream consumers (output_fasta gate, prune)
     * always see a positive value.  Rejected families have mdl_score zeroed.
     * Only exclusive_savings are accumulated into result.total_savings to
     * keep the two-part DL valid (no overlap double-counting).
     */
    if (genome_len <= 0) {
        fprintf(stderr, "mdl_select: invalid genome_len=%" PRId64 "\n",
                (int64_t)genome_len);
        return result;
    }

    /* ENG-N12 (QUALITY_PROPOSAL_v6 Tier 1.5b): replace the
     * O(genome_len/8) covered-bit bitmap with a sorted-merged dynamic
     * array of "covered" intervals.  Memory is O(num_accepted_intervals)
     * — bounded by total instances across accepted families, NOT genome
     * length.  This unblocks wheat/maize-scale runs (where the previous
     * 2 GB / 287 MB single allocation either OOM'd or silently failed
     * via the calloc-NULL path that returned an empty result).
     *
     * Data structure: `covered_iv[]` is sorted by start AND merged
     * (no overlaps, no adjacencies — i.e. covered_iv[k].end < covered_iv[k+1].start).
     *
     * Operations needed:
     *   1. count_covered(s, e): how many bases of [s, e) are already in covered_iv.
     *      Implemented by binary-searching for the first interval with end > s,
     *      then sweeping forward while interval.start < e.
     *   2. add_intervals(family): merge a family's instance intervals
     *      into covered_iv after acceptance.  Implemented as a streaming
     *      two-pointer merge with a side buffer + buffer swap.
     */

    /* MdlInterval is declared at file scope (above mdl_select_library) so
     * that cmp_mdl_interval_start (also file-scope, used as a qsort
     * comparator below) can refer to it. */

    /* Helper: lower_bound on .end (returns smallest i s.t. iv[i].end > key)
     * via binary search.  Equivalent to `first interval that could overlap
     * a query starting at key`. */
    #define MDL_LB_END_GT(iv, n, key, out) do {                              \
        int _lo = 0, _hi = (n);                                              \
        while (_lo < _hi) {                                                  \
            int _mid = _lo + ((_hi - _lo) >> 1);                             \
            if ((iv)[_mid].end <= (key)) _lo = _mid + 1;                     \
            else                         _hi = _mid;                        \
        }                                                                    \
        (out) = _lo;                                                         \
    } while (0)

    MdlInterval *covered_iv = NULL;
    int64_t covered_n = 0;
    /* covered_cap tracks the malloc'd capacity of covered_iv so the swap
     * below correctly migrates the (buf, cap) pair between covered_iv
     * and merge_buf without leaking. */
    int64_t covered_cap = 0;

    /* Standalone-fallback gate (Stage B fix 1).
     *
     * A family with positive standalone MDL score, ≥50 bp consensus, and
     * ≥3 instances is a valid TE family in its own right.  The pure
     * unique-coverage greedy rejects it whenever a higher-scoring family
     * already claimed its territory — destroying ~28% of families and
     * 49% of covered bases on chr4 (REFINE_TRACE_REPORT bottleneck #1).
     *
     * To preserve such families, accept them via condition (b) below
     * even when their exclusive contribution alone is non-positive.  The
     * three thresholds (positive standalone score, length ≥ 50, ≥ 3
     * instances) gate out spurious noise while readmitting overlapped-
     * but-real families. */
    #define STANDALONE_MIN_LEN 50
    #define STANDALONE_MIN_INST 3

    /* Per-family scratch buffers, grown as needed across iterations.
     * `fam_iv` holds the family's instance intervals (clipped to genome,
     * sorted-merged) for both the exclusive-bp query and the merge step. */
    MdlInterval *fam_iv = NULL;
    int64_t fam_iv_cap = 0;
    /* `merge_buf` is the streaming-merge output buffer used by accept(). */
    MdlInterval *merge_buf = NULL;
    int64_t merge_buf_cap = 0;

    /* Comparator for qsort on MdlInterval.start ascending. */
    /* GCC nested-function would be cleaner but is non-portable;
     * use a file-scope static via function-pointer.  Inline here via
     * a local struct + a compare-by-start pattern: since MdlInterval has
     * (start, end), we can sort using a wrapper.  To keep the patch
     * confined to mdl.c without polluting the file scope, declare a
     * static comparator function (forward-declared at file scope). */

    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *fam = &cl->families[i];

        /* Pre-screen: families with non-positive raw score or <2 instances
         * cannot produce a positive exclusive score either.  Reject and zero.
         * (Note: standalone_score is exactly fam->mdl_score at this point,
         * before any rewrite below.) */
        if (fam->mdl_score <= 0 || fam->num_instances < 2) {
            fam->mdl_score = 0.0;
            continue;
        }

        /* Capture standalone score (set by mdl_score_family above) so it
         * remains available after the per-family mdl_score field is
         * rewritten with the exclusive score on accept. */
        double standalone_score = fam->mdl_score;

        /* --- Build the family's clipped, sorted-merged interval list.
         * We compute this regardless of accept/reject because it's the
         * substrate for both the exclusive-bp query and the merge step,
         * and sorting once amortizes well. */
        if ((int64_t)fam->num_instances > fam_iv_cap) {
            int64_t new_cap = fam_iv_cap ? fam_iv_cap * 2 : 64;
            while (new_cap < fam->num_instances) new_cap *= 2;
            MdlInterval *tmp = realloc(fam_iv,
                                       (size_t)new_cap * sizeof(MdlInterval));
            if (!tmp) {
                fprintf(stderr, "mdl_select: OOM growing fam_iv to %" PRId64
                        " entries\n", new_cap);
                fam->mdl_score = 0.0;
                continue;
            }
            fam_iv = tmp;
            fam_iv_cap = new_cap;
        }

        int64_t fam_n = 0;
        for (int j = 0; j < fam->num_instances; j++) {
            Instance *inst = &fam->instances[j];
            gpos_t s = inst->position;
            gpos_t e = s + (gpos_t)inst->aligned_length;
            if (s < 0) s = 0;
            if (e > genome_len) e = genome_len;
            if (e <= s) continue;
            fam_iv[fam_n].start = s;
            fam_iv[fam_n].end   = e;
            fam_n++;
        }
        if (fam_n == 0) {
            fam->mdl_score = 0.0;
            continue;
        }

        /* Sort family intervals by start.  qsort needs a function
         * pointer — defined at file scope below mdl_select_library
         * (cmp_mdl_interval_start). */
        extern int cmp_mdl_interval_start(const void *, const void *);
        qsort(fam_iv, (size_t)fam_n, sizeof(MdlInterval),
              cmp_mdl_interval_start);

        /* In-place merge of overlapping/adjacent family intervals.
         * After this, fam_iv[0..fam_merged_n) is sorted-merged. */
        int64_t fam_merged_n = 0;
        {
            gpos_t cs = fam_iv[0].start;
            gpos_t ce = fam_iv[0].end;
            for (int64_t t = 1; t < fam_n; t++) {
                if (fam_iv[t].start <= ce) {
                    if (fam_iv[t].end > ce) ce = fam_iv[t].end;
                } else {
                    fam_iv[fam_merged_n].start = cs;
                    fam_iv[fam_merged_n].end   = ce;
                    fam_merged_n++;
                    cs = fam_iv[t].start;
                    ce = fam_iv[t].end;
                }
            }
            fam_iv[fam_merged_n].start = cs;
            fam_iv[fam_merged_n].end   = ce;
            fam_merged_n++;
        }

        /* Compute exclusive savings: walk the family's *original*
         * (un-merged) instance list — preserves the per-instance
         * accounting (one literal_cost / mdl_instance_cost_full call
         * per instance) — and for each instance, count uncovered bases
         * via binary search on covered_iv. */
        double exclusive_savings = 0.0;
        int    has_any_exclusive = 0;

        for (int j = 0; j < fam->num_instances; j++) {
            Instance *inst = &fam->instances[j];
            gpos_t start = inst->position;
            int    alen  = (int)inst->aligned_length;
            if (alen <= 0) continue;

            gpos_t qs = start;
            gpos_t qe = start + (gpos_t)alen;
            if (qs < 0) qs = 0;
            if (qe > genome_len) qe = genome_len;
            if (qe <= qs) continue;

            int64_t qlen = (int64_t)(qe - qs);

            /* Count overlap with covered_iv via local sweep starting at
             * the first covered interval whose end > qs. */
            int idx;
            MDL_LB_END_GT(covered_iv, covered_n, qs, idx);
            int64_t overlap = 0;
            for (int64_t t = idx; t < covered_n; t++) {
                if (covered_iv[t].start >= qe) break;
                gpos_t os = covered_iv[t].start > qs ? covered_iv[t].start : qs;
                gpos_t oe = covered_iv[t].end   < qe ? covered_iv[t].end   : qe;
                if (oe > os) overlap += (int64_t)(oe - os);
            }

            int64_t excl = qlen - overlap;
            if (excl < 0) excl = 0; /* defensive */
            int excl_bases = (int)excl;
            if (excl_bases == 0) continue;

            int excl_edits = (int)(inst->divergence * (double)excl_bases + 0.5);
            if (excl_edits < 0) excl_edits = 0;
            if (excl_edits > excl_bases) excl_edits = excl_bases;

            double lit = literal_cost(excl_bases);
            double enc = mdl_instance_cost_full(excl_bases, excl_edits,
                                                fam->consensus_length,
                                                R_estimate);
            exclusive_savings += (lit - enc);
            has_any_exclusive = 1;
        }

        double exclusive_score = exclusive_savings - fam->model_cost;
        int marginal_pass = (has_any_exclusive && exclusive_score > 0.0);

        /* Standalone fallback (Stage B fix 1) — see comment above */
        int standalone_pass = (standalone_score > 0.0 &&
                               fam->consensus_length >= STANDALONE_MIN_LEN &&
                               fam->num_instances    >= STANDALONE_MIN_INST);

        if (!marginal_pass && !standalone_pass) {
            /* Reject: model cost not justified, even allowing for the
             * standalone fallback. */
            fam->mdl_score = 0.0;
            continue;
        }

        /* Accept: choose the score to report and the savings to accumulate.
         *
         * If the family qualifies on its own (exclusive_score > 0), prefer
         * it — that's the rigorous unique-coverage contribution.  Otherwise
         * report the standalone score so downstream `mdl_score > 0` filters
         * still see this family as accepted; accumulate only the exclusive
         * savings (zero if no exclusive bases) so the two-part code total
         * remains a valid bound (no double counting of overlapped bases). */
        if (marginal_pass) {
            fam->mdl_score = exclusive_score;
            result.total_savings += exclusive_savings;
        } else {
            fam->mdl_score = standalone_score;
            if (has_any_exclusive)
                result.total_savings += exclusive_savings;
        }

        result.num_accepted++;
        result.total_model_cost += fam->model_cost;

        /* Streaming merge: covered_iv (sorted-merged) + fam_iv[0..fam_merged_n)
         * (sorted-merged) → merge_buf (sorted-merged).  Two-pointer linear
         * scan; output is the union, count newly-covered bases on the fly. */
        int64_t need_cap = covered_n + fam_merged_n + 4;
        if (need_cap > merge_buf_cap) {
            int64_t new_cap = merge_buf_cap ? merge_buf_cap * 2 : 256;
            while (new_cap < need_cap) new_cap *= 2;
            MdlInterval *tmp = realloc(merge_buf,
                                       (size_t)new_cap * sizeof(MdlInterval));
            if (!tmp) {
                fprintf(stderr, "mdl_select: OOM growing merge_buf to %" PRId64
                        " entries\n", new_cap);
                /* Roll back acceptance — leave covered_iv unchanged.
                 * The family was already counted in num_accepted above; for
                 * correctness, undo it. */
                result.num_accepted--;
                result.total_model_cost -= fam->model_cost;
                if (marginal_pass) result.total_savings -= exclusive_savings;
                else if (has_any_exclusive)
                                   result.total_savings -= exclusive_savings;
                fam->mdl_score = 0.0;
                continue;
            }
            merge_buf = tmp;
            merge_buf_cap = new_cap;
        }

        int64_t a = 0, b = 0, m = 0;
        gpos_t cur_s = 0, cur_e = 0;
        int    have_cur = 0;

        /* Helper macro: append [s,e) to merge_buf, merging if it touches
         * the current run. */
        #define MDL_EMIT_RUN(_s, _e) do {                                    \
            if (!have_cur) { cur_s = (_s); cur_e = (_e); have_cur = 1; }    \
            else if ((_s) <= cur_e) { if ((_e) > cur_e) cur_e = (_e); }     \
            else {                                                          \
                merge_buf[m].start = cur_s;                                  \
                merge_buf[m].end   = cur_e;                                  \
                m++;                                                         \
                cur_s = (_s); cur_e = (_e);                                  \
            }                                                                \
        } while (0)

        while (a < covered_n && b < fam_merged_n) {
            if (covered_iv[a].start <= fam_iv[b].start) {
                MDL_EMIT_RUN(covered_iv[a].start, covered_iv[a].end);
                a++;
            } else {
                MDL_EMIT_RUN(fam_iv[b].start, fam_iv[b].end);
                b++;
            }
        }
        while (a < covered_n)    { MDL_EMIT_RUN(covered_iv[a].start,
                                                 covered_iv[a].end); a++; }
        while (b < fam_merged_n) { MDL_EMIT_RUN(fam_iv[b].start,
                                                 fam_iv[b].end);     b++; }
        if (have_cur) {
            merge_buf[m].start = cur_s;
            merge_buf[m].end   = cur_e;
            m++;
        }
        #undef MDL_EMIT_RUN

        /* Compute newly-covered bp = sum(merged) - sum(prev). */
        int64_t prev_total = 0;
        for (int64_t t = 0; t < covered_n; t++)
            prev_total += covered_iv[t].end - covered_iv[t].start;
        int64_t new_total = 0;
        for (int64_t t = 0; t < m; t++)
            new_total += merge_buf[t].end - merge_buf[t].start;
        result.bases_covered += (new_total - prev_total);

        /* Swap merge_buf into covered_iv (no copy). */
        MdlInterval *swp_iv  = covered_iv;
        int64_t      swp_cap = covered_cap;
        covered_iv  = merge_buf;
        covered_cap = merge_buf_cap;
        covered_n   = m;
        merge_buf     = swp_iv;
        merge_buf_cap = swp_cap;
    }

    #undef MDL_LB_END_GT
    #undef STANDALONE_MIN_LEN
    #undef STANDALONE_MIN_INST

    free(fam_iv);
    free(merge_buf);
    free(covered_iv);

    /* Compute total DL — include L_int(R) for number of families */
    int R_final = result.num_accepted;
    if (R_final < 1) R_final = 1;
    result.dl_library = L_int(R_final) + result.total_model_cost;
    result.dl_total = literal_cost(genome_len) - result.total_savings +
                      result.dl_library;
    /* dl_total is bounded below by 0 in theory; clamp defensively against
     * floating-point round-off and report only meaningful ratios. */
    if (result.dl_total < 0.0) result.dl_total = 0.0;
    result.compression_ratio = result.dl_total / literal_cost(genome_len);

    return result;
}
