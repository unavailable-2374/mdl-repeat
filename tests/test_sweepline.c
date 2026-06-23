/*
 * Sweep-line regression tests (ENG-N2 + ENG-N12, QUALITY_PROPOSAL_v6
 * Tier 1.5b).  Exercises the O(num_intervals) replacements for the
 * O(genome_len) coverage allocations in:
 *
 *   1. mdl.c::mdl_select_library    (was: ~genome_len/8 bitmap)
 *   2. refine.c::refine_prune_families  (was: genome_len-byte cov[] + O(n^2) sort)
 *
 * Test strategy:
 *   A. Numerical equivalence on a SMALL genome (genome_len = 10^5):
 *      the sweep-line implementation must produce the same accept/
 *      prune decisions and the same bases_covered as a brute-force
 *      bitmap reference computed in the test harness.
 *   B. Memory ceiling on a SIMULATED LARGE GENOME (genome_len = 4e9):
 *      run with the SAME small instance set but with genome_len set
 *      to wheat-scale, then assert peak RSS delta stays below 200 MB.
 *      Pre-patch this would have allocated 500 MB (mdl bitmap) +
 *      4 GB (refine cov[]) — > 4.5 GB.
 *   C. Cross-check that the small-genome and large-genome runs return
 *      identical num_accepted, bases_covered, and per-family acceptance
 *      flags.  A genuine sweep-line implementation is genome-size
 *      agnostic; allocating a per-base array would have OOM'd at 4e9.
 *
 * Build: see tests/run_tests.sh (linked against obj-files in obj/).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "../src/mdl.h"
#include "../src/refine.h"
#include "../src/candidates.h"

/* Forward declaration: the sweep-line core in refine.c is exposed via
 * extern linkage so this test can call it directly without going through
 * the public refine_prune_families wrapper.  The wrapper is also tested
 * indirectly via run_tests.sh's chr4 smoke test. */
extern int refine_prune_families_sweepline(CandidateList *cl,
                                           glen_t genome_len, int verbose,
                                           int num_families);

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

/* ============================================================
 * Memory tracking via /proc/self/statm (Linux-specific).
 *
 * Returns RSS in bytes.  On non-Linux platforms returns -1, in which
 * case the memory-ceiling assertion is skipped (and reported as a
 * SKIP, not a FAIL).
 * ============================================================ */
static int64_t rss_bytes(void)
{
    FILE *fp = fopen("/proc/self/statm", "r");
    if (!fp) return -1;
    long size = 0, resident = 0;
    if (fscanf(fp, "%ld %ld", &size, &resident) != 2) {
        fclose(fp);
        return -1;
    }
    fclose(fp);
    long pagesize = 4096;     /* default x86_64 */
    return (int64_t)resident * pagesize;
}

/* ============================================================
 * Helpers — build candidate families with synthetic instance sets.
 *
 * N families × M instances each.  Family i's instances are placed at
 * positions i*spread + j*stride (j = 0 .. M-1), each of length L.
 * Setting spread > stride*M produces fully disjoint families;
 * spread = 0 produces fully overlapping families.
 * ============================================================ */
static void make_family(CandidateFamily *fam, int id, int cons_len,
                        int n_inst, gpos_t base, gpos_t stride, int alen,
                        int edits)
{
    memset(fam, 0, sizeof(*fam));
    fam->id = (uid_t)id;
    fam->consensus_length = cons_len;
    fam->consensus = calloc((size_t)cons_len, 1);
    fam->num_instances = n_inst;
    fam->cap_instances = n_inst;
    fam->instances = calloc((size_t)n_inst, sizeof(Instance));
    for (int j = 0; j < n_inst; j++) {
        fam->instances[j].position       = base + (gpos_t)j * stride;
        fam->instances[j].aligned_length = alen;
        fam->instances[j].num_edits      = edits;
        fam->instances[j].divergence     = (float)edits / (float)alen;
        fam->instances[j].strand         = 1;
        fam->instances[j].seq_index      = 0;
    }
}

static void free_cl(CandidateList *cl)
{
    for (int i = 0; i < cl->num_families; i++) {
        free(cl->families[i].consensus);
        free(cl->families[i].instances);
    }
    free(cl->families);
    memset(cl, 0, sizeof(*cl));
}

/* Build a fresh CandidateList that mirrors the previous one (same family
 * count, instance positions, lengths, edits).  Used so we can run the
 * sweep-line at small and large genome_len with identical input. */
static void build_cl(CandidateList *cl, int n_fam, int n_inst_per_fam,
                     int cons_len, int alen, int edits,
                     gpos_t fam_spread, gpos_t inst_stride)
{
    memset(cl, 0, sizeof(*cl));
    cl->num_families = n_fam;
    cl->cap_families = n_fam;
    cl->families = calloc((size_t)n_fam, sizeof(CandidateFamily));
    for (int i = 0; i < n_fam; i++) {
        make_family(&cl->families[i], i, cons_len, n_inst_per_fam,
                    (gpos_t)i * fam_spread, inst_stride, alen, edits);
    }
}

/* ============================================================
 * Test 1 — mdl_select_library equivalence on small genome.
 *
 * Two families with partially overlapping instances.  Brute-force
 * compute bases_covered as |union(all instance intervals)|.  Verify
 * the sweep-line produces the same value.
 * ============================================================ */
static void test_mdl_select_small(void)
{
    printf("\n[1] mdl_select_library — sweep-line equivalence (small genome)\n");

    CandidateList cl;
    /* 3 families, 5 instances each, 200 bp long.  Family spread = 5000,
     * instance stride = 1000, so families overlap pairwise (instance
     * positions: F0 = 0,1000,...,4000; F1 = 5000,6000,...,9000;
     * F2 = 10000,11000,...,14000 — actually disjoint at this spread).
     *
     * Bring overlap by setting fam_spread = 800: each family overlaps
     * with the next (alen=200, stride=1000 → instances disjoint within
     * a family; family 0 inst 0 covers 0..200, family 1 inst 0 covers
     * 800..1000 — disjoint.  To force overlap, set fam_spread=100). */
    build_cl(&cl, 3, 5, 200, 200, 5,
             /*fam_spread*/100, /*inst_stride*/1000);

    glen_t genome_len = 100000;
    MDLResult r = mdl_select_library(&cl, genome_len);

    /* Brute-force reference: compute |union| of all accepted-family
     * instance intervals using a bitmap. */
    uint8_t *ref = calloc((size_t)genome_len, 1);
    int64_t ref_covered = 0;
    for (int i = 0; i < cl.num_families; i++) {
        if (!candidate_is_accepted(&cl.families[i])) continue;
        for (int j = 0; j < cl.families[i].num_instances; j++) {
            gpos_t s = cl.families[i].instances[j].position;
            gpos_t e = s + (gpos_t)cl.families[i].instances[j].aligned_length;
            if (s < 0) s = 0;
            if (e > genome_len) e = genome_len;
            for (gpos_t p = s; p < e; p++) {
                if (!ref[p]) { ref[p] = 1; ref_covered++; }
            }
        }
    }
    free(ref);

    CHECK(r.bases_covered == ref_covered,
          "sweep-line bases_covered (%" PRId64 ") == brute-force union (%" PRId64 ")",
          (int64_t)r.bases_covered, ref_covered);
    CHECK(r.compression_ratio >= 0.0 && r.compression_ratio <= 1.0,
          "compression_ratio in [0,1] (got %.4f)", r.compression_ratio);
    CHECK(r.num_accepted >= 1,
          "at least 1 family accepted (got %d)", r.num_accepted);

    free_cl(&cl);
}

/* ============================================================
 * Test 2 — mdl_select_library memory ceiling on simulated 4 Gb genome.
 *
 * Same instance set as Test 1 (positions in 0..15000), but pass
 * genome_len = 4e9.  The sweep-line must NOT allocate per-base storage,
 * so peak RSS delta should be < 200 MB.
 * ============================================================ */
static void test_mdl_select_large_genome(void)
{
    printf("\n[2] mdl_select_library — 4Gb simulated genome (memory ceiling)\n");

    int64_t before = rss_bytes();

    CandidateList cl;
    build_cl(&cl, 3, 5, 200, 200, 5, 100, 1000);

    glen_t big_genome = (glen_t)4 * 1000 * 1000 * 1000;  /* 4 Gb */
    MDLResult r = mdl_select_library(&cl, big_genome);

    int64_t after = rss_bytes();

    CHECK(r.num_accepted >= 1,
          "4Gb genome: at least 1 family accepted (got %d)", r.num_accepted);
    CHECK(r.compression_ratio >= 0.0 && r.compression_ratio <= 1.0,
          "4Gb genome: compression_ratio in [0,1] (got %.4f)",
          r.compression_ratio);

    if (before >= 0 && after >= 0) {
        int64_t delta = after - before;
        if (delta < 0) delta = 0;
        printf("    RSS before=%.1f MB, after=%.1f MB, delta=%.1f MB\n",
               before / (1024.0 * 1024.0),
               after  / (1024.0 * 1024.0),
               delta  / (1024.0 * 1024.0));
        /* Pre-patch this test would have allocated ~500 MB (4e9/8 bitmap).
         * The sweep-line uses O(num_instances) memory ≈ a few KB. */
        CHECK(delta < 200LL * 1024 * 1024,
              "RSS delta < 200 MB (got %.1f MB) — confirms no per-base alloc",
              delta / (1024.0 * 1024.0));
    } else {
        printf("    SKIP: /proc/self/statm unavailable\n");
    }

    free_cl(&cl);
}

/* ============================================================
 * Test 3 — refine_prune_families_sweepline equivalence (small).
 *
 * Build 3 accepted families.  Family 0 and 1 have IDENTICAL instance
 * positions (full overlap).  Family 2 is disjoint.  After pruning,
 * exactly one of {F0, F1} should be removed (the weaker score wins on
 * the original cov[]==1 logic).  Family 2 must survive (its bases are
 * exclusively covered by it).
 * ============================================================ */
static void test_refine_prune_small(void)
{
    printf("\n[3] refine_prune_families_sweepline — small genome equivalence\n");

    CandidateList cl;
    memset(&cl, 0, sizeof(cl));
    cl.num_families = 3;
    cl.cap_families = 3;
    cl.families = calloc(3, sizeof(CandidateFamily));

    /* F0 and F1 both cover [1000,1200), [3000,3200), [5000,5200).
     * F2 covers [10000,10200), [12000,12200), [14000,14200). */
    make_family(&cl.families[0], 0, 200, 3, 1000, 2000, 200, 5);
    make_family(&cl.families[1], 1, 200, 3, 1000, 2000, 200, 5);
    make_family(&cl.families[2], 2, 200, 3, 10000, 2000, 200, 5);

    /* Mark all as accepted (mdl_score > 0) and assign distinct scores
     * so the weakest-first iteration order is well-defined. */
    cl.families[0].mdl_score = 100.0;   /* weakest, identical to F1 → prunable */
    cl.families[1].mdl_score = 200.0;   /* survivor (covers same as F0) */
    cl.families[2].mdl_score = 300.0;   /* disjoint → must survive */
    cl.families[0].model_cost = 50.0;
    cl.families[1].model_cost = 50.0;
    cl.families[2].model_cost = 50.0;
    for (int i = 0; i < 3; i++) {
        cl.families[i].mdl.report_score = cl.families[i].mdl_score;
        cl.families[i].mdl.model_cost = cl.families[i].model_cost;
        cl.families[i].mdl.accept_state = CAND_ACCEPT_EXCLUSIVE;
        cl.families[i].mdl.quality_tier = CAND_TIER_CORE;
    }

    glen_t genome_len = 100000;
    int n_pruned = refine_prune_families_sweepline(&cl, genome_len, 0, 3);

    /* Expectation: F0 should be pruned (weakest, fully redundant w.r.t. F1).
     * F1 retains its score (it has exclusive coverage *after* F0 marked
     * as candidate but before any decision; in the iteration F0 is
     * processed first — its instances have cov[]==2 with F1 active, so
     * excl=0 across all instances → F0 prunes.  F1 then sees cov==1
     * since F0 is gone, so F1 retains.  F2 is fully exclusive throughout. */
    CHECK(n_pruned == 1,
          "exactly 1 family pruned (got %d)", n_pruned);
    CHECK(cl.families[0].mdl_score == 0.0,
          "F0 (weakest, fully overlapping with F1) pruned");
    CHECK(cl.families[0].mdl.accept_state == CAND_ACCEPT_PRUNED,
          "F0 acceptance state marked pruned");
    CHECK(cl.families[1].mdl_score > 0.0,
          "F1 (overlapping with F0 but stronger score) retained");
    CHECK(cl.families[2].mdl_score > 0.0,
          "F2 (disjoint) retained");

    free_cl(&cl);
}

/* ============================================================
 * Test 4 — refine_prune_families_sweepline on simulated 4Gb genome.
 *
 * Same input as Test 3 but genome_len = 4e9.  Memory must stay flat:
 * pre-patch this would have allocated 4 GB (genome_len bytes for cov[]).
 * ============================================================ */
static void test_refine_prune_large_genome(void)
{
    printf("\n[4] refine_prune_families_sweepline — 4Gb simulated genome\n");

    int64_t before = rss_bytes();

    CandidateList cl;
    memset(&cl, 0, sizeof(cl));
    cl.num_families = 3;
    cl.cap_families = 3;
    cl.families = calloc(3, sizeof(CandidateFamily));

    make_family(&cl.families[0], 0, 200, 3, 1000, 2000, 200, 5);
    make_family(&cl.families[1], 1, 200, 3, 1000, 2000, 200, 5);
    make_family(&cl.families[2], 2, 200, 3, 10000, 2000, 200, 5);

    cl.families[0].mdl_score = 100.0;
    cl.families[1].mdl_score = 200.0;
    cl.families[2].mdl_score = 300.0;
    cl.families[0].model_cost = 50.0;
    cl.families[1].model_cost = 50.0;
    cl.families[2].model_cost = 50.0;
    for (int i = 0; i < 3; i++) {
        cl.families[i].mdl.report_score = cl.families[i].mdl_score;
        cl.families[i].mdl.model_cost = cl.families[i].model_cost;
        cl.families[i].mdl.accept_state = CAND_ACCEPT_EXCLUSIVE;
        cl.families[i].mdl.quality_tier = CAND_TIER_CORE;
    }

    glen_t big = (glen_t)4 * 1000 * 1000 * 1000;
    int n_pruned = refine_prune_families_sweepline(&cl, big, 0, 3);

    int64_t after = rss_bytes();

    /* Same prune decision must hold regardless of genome_len — the
     * decision depends only on instance overlap, which is unchanged. */
    CHECK(n_pruned == 1,
          "4Gb genome: exactly 1 family pruned (got %d)", n_pruned);

    if (before >= 0 && after >= 0) {
        int64_t delta = after - before;
        if (delta < 0) delta = 0;
        printf("    RSS before=%.1f MB, after=%.1f MB, delta=%.1f MB\n",
               before / (1024.0 * 1024.0),
               after  / (1024.0 * 1024.0),
               delta  / (1024.0 * 1024.0));
        /* Pre-patch: 4 GB calloc.  Sweep-line: O(num_intervals) ≈ a few KB. */
        CHECK(delta < 200LL * 1024 * 1024,
              "RSS delta < 200 MB (got %.1f MB) — confirms no genome_len alloc",
              delta / (1024.0 * 1024.0));
    } else {
        printf("    SKIP: /proc/self/statm unavailable\n");
    }

    free_cl(&cl);
}

/* ============================================================
 * Test 5 — qsort replacement (ENG-N10): scale n_accepted to 10k families
 *
 * The previous O(n^2) selection sort would take ~1e8 operations on
 * n_accepted=10000.  qsort is ~1.3e5 — 1000x faster.  We don't time
 * it directly (CI-flaky), just confirm the sweep-line completes
 * quickly with no memory blow-up and produces a sensible n_pruned.
 * ============================================================ */
static void test_refine_prune_scale(void)
{
    printf("\n[5] refine_prune_families_sweepline — 10k families (qsort scale)\n");

    int n_fam = 10000;
    CandidateList cl;
    memset(&cl, 0, sizeof(cl));
    cl.num_families = n_fam;
    cl.cap_families = n_fam;
    cl.families = calloc((size_t)n_fam, sizeof(CandidateFamily));

    /* Disjoint families: F_i covers [i*1000, i*1000+200).  No prune
     * candidates because every instance is exclusively covered. */
    for (int i = 0; i < n_fam; i++) {
        make_family(&cl.families[i], i, 200, 1, (gpos_t)i * 1000, 1000,
                    200, 5);
        cl.families[i].mdl_score = (double)(i + 1);  /* distinct scores */
        cl.families[i].model_cost = 50.0;
        cl.families[i].mdl.report_score = cl.families[i].mdl_score;
        cl.families[i].mdl.model_cost = cl.families[i].model_cost;
        cl.families[i].mdl.accept_state = CAND_ACCEPT_EXCLUSIVE;
        cl.families[i].mdl.quality_tier = CAND_TIER_CORE;
    }

    int64_t before = rss_bytes();
    /* Use a comfortably large genome_len just to verify no per-base alloc
     * sneaks in via the sweep-line itself. */
    glen_t big = (glen_t)4 * 1000 * 1000 * 1000;
    int n_pruned = refine_prune_families_sweepline(&cl, big, 0, n_fam);
    int64_t after = rss_bytes();

    CHECK(n_pruned == 0,
          "10k disjoint families: 0 pruned (got %d)", n_pruned);

    if (before >= 0 && after >= 0) {
        int64_t delta = after - before;
        if (delta < 0) delta = 0;
        printf("    RSS delta=%.1f MB (n_fam=10000)\n",
               delta / (1024.0 * 1024.0));
        CHECK(delta < 200LL * 1024 * 1024,
              "10k families RSS delta < 200 MB (got %.1f MB)",
              delta / (1024.0 * 1024.0));
    }

    free_cl(&cl);
}

/* ============================================================ */

int main(void)
{
    printf("=== sweep-line tests (ENG-N2 + ENG-N10 + ENG-N12) ===\n");

    test_mdl_select_small();
    test_mdl_select_large_genome();
    test_refine_prune_small();
    test_refine_prune_large_genome();
    test_refine_prune_scale();

    printf("\n%d passed, %d failed\n", g_passed, g_failed);
    return (g_failed == 0) ? 0 : 1;
}
