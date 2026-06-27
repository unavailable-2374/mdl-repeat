/*
 * Unit tests for mdl.c — Rissanen integer code, per-instance cost,
 * model cost, and library selection (including unique-coverage gating
 * introduced in M1#1).
 *
 * Build: see tests/run_tests.sh.  Tests are self-contained — no
 * external deps beyond mdl.c + libm.  Each block prints PASS / FAIL
 * lines and the program exits 0 iff every assertion passes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../src/mdl.h"

static int g_passed = 0;
static int g_failed = 0;

#define CHECK(cond, fmt, ...) do {                                        \
    if (cond) {                                                            \
        printf("  PASS: " fmt "\n", ##__VA_ARGS__);                        \
        g_passed++;                                                        \
    } else {                                                               \
        printf("  FAIL: " fmt " [%s:%d]\n", ##__VA_ARGS__, __FILE__,       \
               __LINE__);                                                  \
        g_failed++;                                                        \
    }                                                                      \
} while (0)

#define CHECK_NEAR(actual, expected, tol, fmt, ...) do {                  \
    double _a = (actual), _e = (expected);                                 \
    int _ok = fabs(_a - _e) < (tol);                                       \
    if (_ok) {                                                             \
        printf("  PASS: " fmt " (got %.4f, expected %.4f)\n",              \
               ##__VA_ARGS__, _a, _e);                                     \
        g_passed++;                                                        \
    } else {                                                               \
        printf("  FAIL: " fmt " (got %.4f, expected %.4f, tol=%g) "        \
               "[%s:%d]\n",                                                \
               ##__VA_ARGS__, _a, _e, (double)(tol),                       \
               __FILE__, __LINE__);                                        \
        g_failed++;                                                        \
    }                                                                      \
} while (0)

/* ============================================================
 * 1. Rissanen L_int reference values
 * ============================================================ */
static void test_L_int(void)
{
    printf("\n[1] L_int reference values\n");
    struct { int64_t n; double expected; } cases[] = {
        {1,    1.518},
        {2,    2.518},
        {3,    3.768},
        {10,   7.364},
        {100, 12.880},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        CHECK_NEAR(L_int(cases[i].n), cases[i].expected, 0.01,
                   "L_int(%lld)", (long long)cases[i].n);
    }
    /* Monotonicity: L_int strictly grows for n > 1 */
    CHECK(L_int(2) > L_int(1), "L_int monotonic 2 > 1");
    CHECK(L_int(1000) > L_int(100), "L_int monotonic 1000 > 100");
    /* Boundary: L_int(<= 0) returns INFINITY (not 0) per M3#9 */
    CHECK(isinf(L_int(0)),  "L_int(0) returns INFINITY");
    CHECK(isinf(L_int(-5)), "L_int(-5) returns INFINITY");
}

/* ============================================================
 * 2. mdl_model_cost — monotonicity + invalid input
 * ============================================================ */
static void test_model_cost(void)
{
    printf("\n[2] mdl_model_cost\n");
    double c10  = mdl_model_cost(10);
    double c100 = mdl_model_cost(100);
    double c1k  = mdl_model_cost(1000);
    CHECK(c10 > 0, "model_cost(10) positive (%.2f)", c10);
    CHECK(c100 > c10, "model_cost monotonic: c100 > c10 (%.2f > %.2f)",
          c100, c10);
    CHECK(c1k > c100, "model_cost monotonic: c1k > c100");
    /* Asymptotically dominated by 2*L term */
    CHECK_NEAR(c100 / 100.0, 2.0, 0.5,
               "model_cost(L)/L approaches 2 bits/base for large L");
    /* Invalid input returns INFINITY (M3#9) */
    CHECK(isinf(mdl_model_cost(0)),  "model_cost(0) returns INFINITY");
    CHECK(isinf(mdl_model_cost(-1)), "model_cost(-1) returns INFINITY");
}

/* ============================================================
 * 3. mdl_instance_cost_full — boundary m, three modes
 * ============================================================ */
static void test_instance_cost(void)
{
    printf("\n[3] mdl_instance_cost_full — three modes × boundary edits\n");

    /* Reference: a=200, m=0 — perfect match */
    g_mdl_mode = MDL_MODE_NONE;
    double e_none_m0  = mdl_instance_cost_full(200, 0, 200, 100);
    g_mdl_mode = MDL_MODE_EXACT;
    double e_exact_m0 = mdl_instance_cost_full(200, 0, 200, 100);
    g_mdl_mode = MDL_MODE_UPPER;
    double e_upper_m0 = mdl_instance_cost_full(200, 0, 200, 100);

    /* m=0 means no position info needed: all three modes equal */
    CHECK_NEAR(e_none_m0,  e_exact_m0, 0.01, "m=0: NONE == EXACT");
    CHECK_NEAR(e_exact_m0, e_upper_m0, 0.01, "m=0: EXACT == UPPER");
    /* perfect match cost much less than 2*a=400 (literal) */
    CHECK(e_exact_m0 < 50.0, "m=0 perfect match cost (%.2f) << literal 400",
          e_exact_m0);

    /* a=200, m=20 (10% divergence): UPPER should >= EXACT, both >= NONE */
    g_mdl_mode = MDL_MODE_NONE;
    double e_none_m20  = mdl_instance_cost_full(200, 20, 200, 100);
    g_mdl_mode = MDL_MODE_EXACT;
    double e_exact_m20 = mdl_instance_cost_full(200, 20, 200, 100);
    g_mdl_mode = MDL_MODE_UPPER;
    double e_upper_m20 = mdl_instance_cost_full(200, 20, 200, 100);

    CHECK(e_upper_m20 >= e_exact_m20 - 0.01,
          "10%% div: UPPER (%.2f) >= EXACT (%.2f) [m*log2(a) >= log2(C(a,m))]",
          e_upper_m20, e_exact_m20);
    CHECK(e_exact_m20 > e_none_m20,
          "10%% div: EXACT (%.2f) > NONE (%.2f) [adds position term]",
          e_exact_m20, e_none_m20);

    /* a=200, m=200 (everything edited — degenerate): EXACT and UPPER
     * still finite, NONE = m*log2(3) + log terms */
    g_mdl_mode = MDL_MODE_EXACT;
    double e_exact_full = mdl_instance_cost_full(200, 200, 200, 100);
    CHECK(isfinite(e_exact_full), "m==a EXACT finite (%.2f)", e_exact_full);

    /* Restore default mode for downstream tests */
    g_mdl_mode = MDL_MODE_EXACT;

    /* Monotonicity in m: cost grows with edits up to m = a/2 (for EXACT) */
    double e_m5  = mdl_instance_cost_full(200,  5, 200, 100);
    double e_m20_chk = mdl_instance_cost_full(200, 20, 200, 100);
    CHECK(e_m20_chk > e_m5, "EXACT cost monotonic in m for m << a/2");

    /* Negative m / m > a should be clamped — function must not crash */
    double e_neg = mdl_instance_cost_full(200, -5, 200, 100);
    double e_overflow = mdl_instance_cost_full(200, 999, 200, 100);
    CHECK(isfinite(e_neg),      "negative m clamped, finite cost");
    CHECK(isfinite(e_overflow), "m > a clamped, finite cost");
}

/* ============================================================
 * 4. mdl_select_library — unique-coverage exclusivity (M1#1 fix)
 * ============================================================ */

/* Helper: build a fresh CandidateFamily on the heap, attach instances. */
static void make_family(CandidateFamily *fam, int id, int cons_len,
                        int n_inst, gpos_t starts[], int aligned_lens[],
                        int edits[])
{
    memset(fam, 0, sizeof(*fam));
    fam->id = (mdl_uid_t)id;
    fam->consensus_length = cons_len;
    fam->consensus = calloc((size_t)cons_len, 1);  /* zeros are valid 'A' */
    fam->num_instances = n_inst;
    fam->cap_instances = n_inst;
    fam->instances = calloc((size_t)n_inst, sizeof(Instance));
    for (int i = 0; i < n_inst; i++) {
        fam->instances[i].position       = starts[i];
        fam->instances[i].aligned_length = aligned_lens[i];
        fam->instances[i].num_edits      = edits[i];
        fam->instances[i].divergence     = (float)edits[i] / aligned_lens[i];
        fam->instances[i].strand         = 1;
        fam->instances[i].seq_index      = 0;
    }
}

static void free_family(CandidateFamily *fam)
{
    free(fam->consensus);
    free(fam->instances);
}

static void test_select_unique_coverage(void)
{
    printf("\n[4] mdl_select_library — unique-coverage exclusivity "
           "(standalone fallback)\n");

    /* Two families, identical instances at the same coordinates.
     *
     * Naive (pre-M1#1) summing would double-count the savings,
     * potentially dropping compression_ratio < 0.  After the M1#1
     * unique-coverage gate, the loser had mdl_score zeroed.  After Stage B
     * fix 1 (standalone fallback for cons >= 50bp + >= 3 instances),
     * both families pass the gate — each is a "valid TE family in its
     * own right" — but bases_covered must still count each unique
     * position exactly once and compression_ratio must remain in [0,1].
     */
    CandidateList cl;
    memset(&cl, 0, sizeof(cl));
    cl.num_families = 2;
    cl.cap_families = 2;
    cl.families = calloc(2, sizeof(CandidateFamily));

    gpos_t starts[3]  = {1000, 5000, 9000};
    int    alens[3]   = {200,  200,  200};
    int    edits_a[3] = {2,    3,    1};   /* high-quality family */
    int    edits_b[3] = {15,   12,   18};  /* lower-quality, same positions */

    make_family(&cl.families[0], 0, 200, 3, starts, alens, edits_a);
    make_family(&cl.families[1], 1, 200, 3, starts, alens, edits_b);

    glen_t genome_len = 100000;
    MDLResult r = mdl_select_library(&cl, genome_len);

    /* Stage B: both families are valid (cons >= 50, n >= 3, standalone > 0)
     * so both should be accepted via the fallback gate. */
    CHECK(r.num_accepted == 2,
          "two families on identical positions, both valid standalone → "
          "both accepted via fallback gate (got %d)", r.num_accepted);
    CHECK(r.compression_ratio >= 0.0 && r.compression_ratio <= 1.0,
          "compression_ratio in [0,1] (got %.4f)", r.compression_ratio);
    /* No family should have its score zeroed: both pass standalone gate. */
    int n_pos = 0, n_zero = 0, n_exclusive = 0, n_standalone = 0;
    int n_warn = 0, n_fallback_flag = 0, n_no_exclusive_flag = 0;
    for (int i = 0; i < cl.num_families; i++) {
        if (cl.families[i].mdl_score >  0.0) n_pos++;
        if (cl.families[i].mdl_score == 0.0) n_zero++;
        if (cl.families[i].mdl.accept_state == CAND_ACCEPT_EXCLUSIVE) n_exclusive++;
        if (cl.families[i].mdl.accept_state == CAND_ACCEPT_STANDALONE) n_standalone++;
        if (cl.families[i].mdl.quality_tier == CAND_TIER_WARN) n_warn++;
        if (cl.families[i].mdl.quality_flags & CAND_QF_STANDALONE_FALLBACK)
            n_fallback_flag++;
        if (cl.families[i].mdl.quality_flags & CAND_QF_NO_EXCLUSIVE_BASES)
            n_no_exclusive_flag++;
    }
    CHECK(n_pos == 2 && n_zero == 0,
          "score: 2 positive + 0 zero (got %d / %d)", n_pos, n_zero);
    CHECK(n_exclusive == 1 && n_standalone == 1,
          "acceptance provenance: 1 exclusive + 1 standalone (got %d / %d)",
          n_exclusive, n_standalone);
    CHECK(n_warn == 1 && n_fallback_flag == 1 && n_no_exclusive_flag == 1,
          "quality provenance: standalone fallback marked warn with no-exclusive flag");

    /* bases_covered must NOT double-count overlap.  Total unique
     * positions = 3 * 200 = 600 (the two families occupy identical
     * coordinates). */
    CHECK(r.bases_covered == 600,
          "bases_covered == 600 (no double-counting; got %lld)",
          (long long)r.bases_covered);

    free_family(&cl.families[0]);
    free_family(&cl.families[1]);
    free(cl.families);
}

/*
 * Verify the unique-coverage gate STILL rejects families that fail BOTH
 * (a) marginal coverage AND (b) the standalone fallback (cons >= 50 +
 * n_inst >= 3).  Use n_instances = 2 to disqualify condition (b).
 */
static void test_select_unique_coverage_no_fallback(void)
{
    printf("\n[4b] mdl_select_library — unique-coverage rejects when "
           "standalone fallback also fails\n");

    CandidateList cl;
    memset(&cl, 0, sizeof(cl));
    cl.num_families = 2;
    cl.cap_families = 2;
    cl.families = calloc(2, sizeof(CandidateFamily));

    /* Two instances per family — below the standalone-fallback threshold
     * of 3, so the only path to acceptance is via marginal exclusive
     * coverage.  Identical positions → loser must be rejected. */
    gpos_t starts[2]  = {1000, 5000};
    int    alens[2]   = {200,  200};
    int    edits_a[2] = {2,    3};
    int    edits_b[2] = {15,   12};

    make_family(&cl.families[0], 0, 200, 2, starts, alens, edits_a);
    make_family(&cl.families[1], 1, 200, 2, starts, alens, edits_b);

    MDLResult r = mdl_select_library(&cl, 100000);

    CHECK(r.num_accepted == 1,
          "n_inst=2 (below standalone-fallback threshold), "
          "identical positions → exactly 1 accepted (got %d)", r.num_accepted);
    int n_pos = 0, n_zero = 0, n_exclusive = 0, n_rejected = 0;
    int rejected_preserved_standalone = 0;
    for (int i = 0; i < cl.num_families; i++) {
        if (cl.families[i].mdl_score >  0.0) n_pos++;
        if (cl.families[i].mdl_score == 0.0) n_zero++;
        if (cl.families[i].mdl.accept_state == CAND_ACCEPT_EXCLUSIVE) n_exclusive++;
        if (cl.families[i].mdl.accept_state == CAND_ACCEPT_REJECTED) {
            n_rejected++;
            if (cl.families[i].mdl.standalone_score > 0.0)
                rejected_preserved_standalone++;
        }
    }
    CHECK(n_pos == 1 && n_zero == 1,
          "score rewrite: 1 positive + 1 zero (got %d / %d)", n_pos, n_zero);
    CHECK(n_exclusive == 1 && n_rejected == 1,
          "acceptance provenance: 1 exclusive + 1 rejected (got %d / %d)",
          n_exclusive, n_rejected);
    CHECK(rejected_preserved_standalone == 1,
          "rejected overlapped family preserves positive standalone_score");
    CHECK(r.bases_covered == 400,
          "bases_covered == 400 (got %lld)", (long long)r.bases_covered);

    free_family(&cl.families[0]);
    free_family(&cl.families[1]);
    free(cl.families);
}

static void test_select_disjoint_families(void)
{
    printf("\n[5] mdl_select_library — both families accepted when disjoint\n");

    CandidateList cl;
    memset(&cl, 0, sizeof(cl));
    cl.num_families = 2;
    cl.cap_families = 2;
    cl.families = calloc(2, sizeof(CandidateFamily));

    gpos_t starts_a[3] = {1000, 5000, 9000};
    gpos_t starts_b[3] = {20000, 30000, 40000};   /* disjoint */
    int alens[3]  = {200, 200, 200};
    int edits[3]  = {2, 3, 1};

    make_family(&cl.families[0], 0, 200, 3, starts_a, alens, edits);
    make_family(&cl.families[1], 1, 200, 3, starts_b, alens, edits);

    MDLResult r = mdl_select_library(&cl, 100000);

    CHECK(r.num_accepted == 2,
          "disjoint families: both accepted (got %d)", r.num_accepted);
    CHECK(r.bases_covered == 1200,
          "bases_covered == 1200 (got %lld)", (long long)r.bases_covered);
    CHECK(r.compression_ratio >= 0.0 && r.compression_ratio <= 1.0,
          "compression_ratio in [0,1] (got %.4f)", r.compression_ratio);
    CHECK(cl.families[0].mdl.quality_tier == CAND_TIER_CORE &&
          cl.families[1].mdl.quality_tier == CAND_TIER_CORE,
          "disjoint high-quality families remain core tier");
    CHECK(cl.families[0].mdl.quality_flags == CAND_QF_NONE &&
          cl.families[1].mdl.quality_flags == CAND_QF_NONE,
          "disjoint high-quality families have no quality flags");

    free_family(&cl.families[0]);
    free_family(&cl.families[1]);
    free(cl.families);
}

/* ============================================================ */

int main(void)
{
    printf("=== mdl.c unit tests (M3#10) ===\n");

    test_L_int();
    test_model_cost();
    test_instance_cost();
    test_select_unique_coverage();
    test_select_unique_coverage_no_fallback();
    test_select_disjoint_families();

    printf("\n%d passed, %d failed\n", g_passed, g_failed);
    return (g_failed == 0) ? 0 : 1;
}
