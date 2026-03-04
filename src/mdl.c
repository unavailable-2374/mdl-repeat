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
    if (n <= 0) {
        fprintf(stderr, "L_int: n must be positive (got %" PRId64 ")\n", n);
        return 0.0;
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
 */
static inline double literal_cost(int length)
{
    return 2.0 * length;
}

/*
 * Per-instance encoding cost.
 *
 * Retains the original sensitivity-preserving formula:
 *   L_int(a_i) + L_int(m_i + 1) + m_i * log2(3) + [position encoding]
 *
 * The consensus_length and num_families parameters are accepted for API
 * compatibility but not used in the cost formula — the DESIGN_DOC overhead
 * terms (type bit, family ID, strand, consensus pointer) were found to be
 * too aggressive, causing unacceptable sensitivity loss on real genomes.
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
        int m_i = inst->num_edits;  /* exact edits from alignment */
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

MDLResult mdl_select_library(CandidateList *cl, glen_t genome_len)
{
    MDLResult result;
    memset(&result, 0, sizeof(result));

    /* Iterative R convergence: R and scores are circularly dependent.
     * Start with R = num_families (upper bound), iterate up to 3 times. */
    int R_estimate = cl->num_families;
    if (R_estimate < 2) R_estimate = 2;

    for (int iter = 0; iter < 3; iter++) {
        /* Score all families with current R estimate */
        for (int i = 0; i < cl->num_families; i++) {
            mdl_score_family(&cl->families[i], genome_len, R_estimate);
        }

        /* Count accepted */
        int new_R = 0;
        for (int i = 0; i < cl->num_families; i++) {
            if (cl->families[i].mdl_score > 0 &&
                cl->families[i].num_instances >= 2)
                new_R++;
        }
        if (new_R < 2) new_R = 2;

        if (new_R == R_estimate)
            break; /* converged */

        R_estimate = new_R;
    }

    /* Sort by MDL score descending */
    qsort(cl->families, (size_t)cl->num_families, sizeof(CandidateFamily),
          cmp_mdl_desc);

    /* Accept all families with positive raw MDL score and ≥2 instances */
    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *fam = &cl->families[i];

        if (fam->mdl_score <= 0)
            break; /* sorted: all remaining are non-positive */

        if (fam->num_instances < 2)
            continue;

        result.num_accepted++;
        result.total_model_cost += fam->model_cost;

        /* Accumulate raw savings */
        double fam_savings = fam->mdl_score + fam->model_cost; /* undo subtraction */
        result.total_savings += fam_savings;
    }

    /* Count unique coverage via bitmap (for reporting only) */
    size_t bitmap_bytes = ((size_t)genome_len + 7) / 8;
    uint8_t *covered = calloc(bitmap_bytes, 1);
    if (!covered) {
        fprintf(stderr, "mdl_select: could not allocate %" PRId64 " bytes for bitmap\n",
                (int64_t)bitmap_bytes);
        return result;
    }

    #define BIT_SET(pos) (covered[(pos) >> 3] |= (1 << ((pos) & 7)))
    #define BIT_GET(pos) (covered[(pos) >> 3] & (1 << ((pos) & 7)))

    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *fam = &cl->families[i];
        if (fam->mdl_score <= 0) break;
        if (fam->num_instances < 2) continue;

        for (int j = 0; j < fam->num_instances; j++) {
            Instance *inst = &fam->instances[j];
            gpos_t start = inst->position;
            int alen = (int)inst->aligned_length;

            for (int p = 0; p < alen; p++) {
                gpos_t gp = start + p;
                if (gp >= 0 && gp < genome_len) {
                    if (!BIT_GET(gp))
                        result.bases_covered++;
                    BIT_SET(gp);
                }
            }
        }
    }

    #undef BIT_SET
    #undef BIT_GET
    free(covered);

    /* Compute total DL — include L_int(R) for number of families */
    int R_final = result.num_accepted;
    if (R_final < 1) R_final = 1;
    result.dl_library = L_int(R_final) + result.total_model_cost;
    result.dl_total = literal_cost((int)genome_len) - result.total_savings +
                      result.dl_library;
    result.compression_ratio = result.dl_total / literal_cost((int)genome_len);

    return result;
}
