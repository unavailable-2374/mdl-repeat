#ifndef MDL_MDL_H
#define MDL_MDL_H

#include "types.h"
#include "candidates.h"
#include "genome.h"

/*
 * MDL (Minimum Description Length) scoring for repeat families.
 *
 * Uses Rissanen's universal integer code and two-part coding scheme:
 *   DL(total) = DL(Library) + DL(Genome | Library)
 *
 * A family is accepted if its MDL score (savings - model_cost) > 0.
 */

/* --- MDL mode for position encoding --- */
typedef enum {
    MDL_MODE_NONE  = 0,  /* no position encoding (original behavior) */
    MDL_MODE_EXACT = 1,  /* lgamma exact binomial C(a_i, m_i) (default) */
    MDL_MODE_UPPER = 2   /* upper bound: m_i * log2(a_i) */
} MdlMode;

extern MdlMode g_mdl_mode;

/* --- Rissanen's universal integer code --- */

/*
 * L_int(n): encode a positive integer n using Rissanen's universal code.
 * L_int(n) = log2*(n) + log2(c0)
 * where log2*(n) = log2(n) + log2(log2(n)) + ... (sum of positive terms)
 * and c0 ≈ 2.865064 (so log2(c0) ≈ 1.5180...)
 *
 * Returns: number of bits needed to encode n.
 */
double L_int(int64_t n);

/* --- MDL cost helpers (used by refine module) --- */

/*
 * Per-instance encoding cost (extended API).
 *
 * Uses the original sensitivity-preserving formula:
 *   L_int(a_i) + L_int(m_i+1) + m_i*log2(3) + [position encoding]
 *
 * The consensus_length and num_families parameters are accepted for API
 * compatibility (used by iterative R convergence and fragment assembly)
 * but do not affect the per-instance cost formula.
 */
double mdl_instance_cost_full(int aligned_length, int edits,
                              int consensus_length, int num_families);

/*
 * Legacy wrapper (backward compat for modes that don't track R/consensus_length).
 * Uses R=2, consensus_length=aligned_length as conservative defaults.
 */
double mdl_instance_cost(int aligned_length, int edits);

/*
 * Model cost for a family: L_int(consensus_length) + 2*consensus_length
 */
double mdl_model_cost(int consensus_length);

/* --- Per-family MDL scoring --- */

/*
 * Compute MDL score for a single candidate family.
 * Fills fam->mdl_score and fam->model_cost.
 *
 * Parameters:
 *   fam          - candidate family with consensus and instances
 *   genome_len   - total genome length (for literal cost estimation)
 *   num_families - R, number of accepted families (for per-instance overhead)
 */
void mdl_score_family(CandidateFamily *fam, glen_t genome_len, int num_families);

/* --- Library-level selection --- */

typedef struct {
    double total_model_cost;    /* sum of accepted family model costs */
    double total_savings;       /* sum of accepted family savings */
    double dl_library;          /* total DL of the repeat library */
    double dl_genome_given_lib; /* DL of genome given library */
    double dl_total;            /* dl_library + dl_genome_given_lib */
    double compression_ratio;   /* dl_total / (2 * genome_length) */
    int    num_accepted;        /* number of accepted families */
    int64_t bases_covered;      /* total genome bases covered by instances */
} MDLResult;

/*
 * MDL-based library selection with iterative R convergence.
 *
 * R (number of accepted families) and per-instance costs are circularly
 * dependent: scores depend on ceil(log2(R)), but R depends on scores.
 * We iterate up to 3 times until R converges.
 *
 * 1. Score all families with current R estimate
 * 2. Count accepted (score > 0 and >= 2 instances)
 * 3. If accepted count changed, re-score with new R; else converge
 * 4. Compute unique coverage via bitmap (reporting only)
 *
 * Library cost includes L_int(R) for encoding the number of families.
 *
 * Modifies candidates in place (reorders, updates scores).
 * Returns the MDL result summary.
 */
MDLResult mdl_select_library(CandidateList *cl, glen_t genome_len);

#endif /* MDL_MDL_H */
