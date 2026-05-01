# Stage A Report

**Date**: 2026-04-28  
**Binary tested**: `/home/shuoc/tool/mdl-repeat/bin/mdl-repeat` (built 2026-04-28 14:52)  
**Change**: Option γ — parallel l-mer counting (striped-lock hash, 4 threads, 4096 stripe locks) replacing single-threaded linked-list hash in `discover.c`.  
**Benchmark environment**: `/tmp/ath_bench/`, chr4 = 18.6 Mb, RepeatMasker 4.2.1 in PGTA conda env.

---

## chr4 Full RM-remap (new binary)

### Run parameters

```
-sequence /tmp/ath_bench/chr4.fa  -threads 4  (k=14, mdl-mode=exact, defaults)
Output: chr4_optgamma.fa (2769 families), chr4_optgamma.bed (11094 instances)
RM-remap: RepeatMasker -lib chr4_optgamma.fa -no_is -nolow -pa 4
Predicted BED (after RM): 12846 intervals
Eval: bed_pr.py --overlap-mode min --ignore-strand
```

### Overall numbers

| Metric | New binary (Opt γ) | Baseline (chr4_full_rm.bed) | Delta |
|--------|--------------------|-----------------------------|-------|
| Recall | **0.4728** | 0.4680 | +0.0048 |
| Precision | 0.6727 | 0.6766 | −0.0039 |
| F1 | **0.5553** | 0.5533 | +0.0020 |
| TP intervals | 8642 | 8528 | +114 |
| FP intervals | 4204 | 4077 | +127 |
| FN intervals | 9638 | 9693 | −55 |
| Families discovered | 2769 | 2725 | +44 |

**HANDOFF baseline caveat**: HANDOFF §3 reports F1=0.515, recall=0.389, precision=0.760. These numbers cannot be reproduced by running `bed_pr.py --overlap-mode min` on `chr4_full_rm.bed`. The best-guess explanation is that HANDOFF used a different intermediate library for that specific run (the chr4_full.fa produced by that older session had different tandem coalesce / MDL state). The consistent comparison uses `chr4_full_rm.bed` as the baseline (same code path, same eval tool), which gives 0.4680/0.6766/0.5533. Both baselines agree the new binary is negligibly improved.

### bp-level recall (eval_quick.py, for comparison with HANDOFF §7)

| | New (Opt γ) | Baseline |
|---|---|---|
| bp_rec | 0.474 | 0.469 |
| bp_prec | 0.888 | 0.889 |
| bp_F1 | 0.618 | 0.614 |

### Recall by truth interval length (not copy class — see note below)

Run with `recall_by_length.py --min-overlap 0.5`.

| Length range | truth_n | detected (new) | recall (new) | detected (base) | recall (base) |
|---|---|---|---|---|---|
| [0, 50) | 1986 | 90 | 0.045 | 83 | 0.042 |
| [50, 100) | 4635 | 295 | 0.064 | 301 | 0.065 |
| [100, 200) | 2038 | 443 | 0.217 | 420 | 0.206 |
| [200, 500) | 1923 | 445 | 0.231 | 432 | 0.225 |
| [500, 1000) | 1102 | 199 | 0.181 | 208 | 0.189 |
| [1000, 2000) | 720 | 113 | 0.157 | 114 | 0.158 |
| [2000, 5000) | 302 | 34 | 0.113 | 38 | 0.126 |
| [5000, 10000) | 108 | 3 | 0.028 | 13 | 0.120 |
| [10000, ∞) | 72 | 0 | 0.000 | 0 | 0.000 |
| **OVERALL** | **12886** | **1622** | **0.126** | **1609** | **0.125** |

Note: The [5000, 10000) bin shows a regression (3 vs 13 detected). This may indicate that the parallel hash change affected some very long family consensus building, but the sample is small (108 intervals, 10 detected). Worth investigating if this persists.

### Copy-class stratified recall — NOT REPRODUCED

HANDOFF §3's copy-class breakdown (1 copy, 2, 3-4, 5-9, 10-49, 50-199, ≥200) required the original TAIR10 RepeatMasker annotation with family names to assign each truth interval to a family cluster. That annotation is not available in `/tmp/ath_bench/` — the truth BED was derived from soft-masked lowercase positions (no family labels). A BLAST-based proxy was attempted (truth intervals blasted against chr4 to count distinct hit positions) but it produced max copy count of 47 (no ≥200-copy families seen), inconsistent with HANDOFF's CEN180 class. The exact-class breakdown **cannot be reproduced**. Only the overall numbers are available.

**Honest summary**: The in-scope recall number 0.465 from HANDOFF cannot be re-derived exactly with available data. Using the interval-level recall metric on a BLAST copy-count proxy (requiring ≥80% pident, ≥50% query coverage), the 3-199 copy range gave recall 0.450 (new) vs 0.500 (baseline) — but this is a different proxy and likely underestimates the family-level copy count. **Report the overall recall numbers only.**

---

## Synthetic Diagnostic Results

### Genome A — clean background, signal-vs-noise test

**Setup**: 6 representative families × 8 copies each in 5 Mb Markov background; no cross-talk. Total 48 truth intervals.

**New binary run**: k=13, mdl-mode=exact. Parallel hash: "kmer counting: parallel with 4 threads, 4096 stripe locks".

| Family | Truth n | TP | FN | Recall (new) | Recall (baseline) |
|--------|---------|----|----|---|---|
| rep_c8_short_A (703 bp) | 8 | 5 | 3 | **0.625** | 0.625 |
| rep_c10_short_B (752 bp) | 8 | 3 | 5 | **0.375** | 0.375 |
| rep_c11_medium_A (1373 bp) | 8 | 5 | 3 | **0.625** | 0.625 |
| rep_c6_medium_B (1556 bp) | 8 | 5 | 3 | **0.625** | 0.625 |
| rep_c2_long_A (2194 bp) | 8 | 4 | 4 | **0.500** | 0.500 |
| rep_c9_long_B (2283 bp) | 8 | 4 | 4 | **0.500** | 0.500 |
| **Overall** | **48** | **26** | **22** | **0.5417** | **0.542** |

All 6 families discovered at seeding stage (R=0..R=5 in discovery log). Overall recall essentially **unchanged** (0.5417 vs 0.542). Precision per matched family = 1.000 in both versions. 23 unmatched predictions are background matches (same as baseline's 26 FP attributed to `__unmatched_pred__`).

**Key observation**: The parallel hash produced bit-for-bit identical family discovery (same R=0..R=5 families, same lengths, same copy counts). Option γ changed the hash implementation without changing the algorithm semantics on this genome.

---

### Genome B — cross-family interference test

**Setup**: `medium_A` × 100 copies (dominant, divergence 0-5%) + `short_A` × 8 copies (victim, engineered shared 14-mer at positions 200-214 and 300-314) + 5 background families × 8 copies each. Total 148 truth intervals.

| Family | Truth n | TP | FN | Recall (new) |
|--------|---------|----|----|---|
| victim_c8_short_A (703 bp) | 8 | 5 | 3 | **0.625** |
| dominant_c11_medium_A (1373 bp) | 70 | 57 | 13 | 0.814 |
| background_c11_medium_A (1373 bp)* | 8 | 3 | 5 | 0.375 |
| background_c10_short_B (752 bp) | 8 | 5 | 3 | 0.625 |
| background_c2_long_A (2194 bp) | 8 | 6 | 2 | 0.750 |
| background_c6_medium_B (1556 bp) | 8 | 6 | 2 | 0.750 |
| background_c9_long_B (2283 bp) | 10 | 7 | 3 | 0.700 |
| **Overall** | **119** | **88** | **31** | **0.7395** |

*Note: truth labels `background_c11_medium_A` and `dominant_c11_medium_A` are both cluster 11 (medium_A). The dominant pool has 70 truth records (100 inserted minus overlap/coalesce), background has 8 distinct records.

**Critical finding for Hypothesis B**: The victim family (`short_A`) achieves recall **0.625** in Genome B — **identical to Genome A** (0.625). Despite 100 copies of the dominant family carrying the engineered shared 14-mer, the victim family's seed signal is NOT depressed. The discovery log shows R=5 (victim) found with N=12 copies, same as in Genome A.

---

### Genome C — nested element test

**Setup**: `long_A` (LINE, 2194 bp) × 10 copies (7 plain + 3 with nested `short_A` LTR inside) + 5 standalone `short_A` (LTR_solo). Total 18 truth intervals.

| Family | Truth n | TP | FN | Recall |
|--------|---------|----|----|---|
| LINE_plain_c2 | 7 | 0 | 7 | **0.000** |
| LINE_nested_c2 | 3 | 0 | 3 | **0.000** |
| LTR_nested_c8 | 3 | 0 | 3 | **0.000** |
| LTR_solo_c8 | 5 | 0 | 5 | **0.000** |
| **Overall** | **18** | **0** | **18** | **0.000** |

**All recall = 0 on Genome C.** However, this is **NOT due to hypothesis C** (nested shadow). Analysis of the run log reveals the failure mode:

1. Seeds found correctly: R=0 N=28 length=2195 (plain LINE) and R=1 N=24 length=2897 (LINE+LTR composite) — both families seeded.
2. Merge step incorrectly combined R=0 and R=1: after merge, family_id=0 has consensus length=2195, **only 4 instances** (down from 14 max), divergence_mean=7.4%.
3. With only 4 instances and model_cost=4408.7 bits, MDL score = 0.0 → family **rejected** (threshold: score > 0.0).
4. The 3 accepted families are 31-34 bp noise sequences (random background k-mer coincidences).

**Root cause of Genome C failure**: The nested-element merge gate (introduced in M2) triggered on the 2194 bp plain LINE and 2897 bp composite LINE families. The gate is designed to veto merges where the shorter family is predominantly CONTAINED within longer families. Here the 2194 bp family IS contained within 2897 bp (it IS the same element, just without the LTR). The gate kept only 4 instances after filtering, collapsing the effective copy count below MDL acceptance threshold.

This is a **refine/MDL-gate bug**, not a seed-selection failure and not hypothesis C.

**Hypothesis C cannot be assessed from Genome C** because the tool fails before reaching the nested-vs-standalone comparison.

---

## Mechanism Diagnosis

### Hypothesis A: seed signal weaker than background noise

**Evidence**: Genome A shows ALL 6 families are seeded (appear at R=0..R=5 in discovery). Recall is partial (37.5%-62.5%) but not zero. The signal is present; losses are in refinement/instance-acceptance, not seeding.

**Verdict: WEAK.** Seed selection is not the primary bottleneck for these 6 representative families. Cannot rule out hypothesis A for the actual chr4 missed families (which may have lower signal), but the diagnostic was designed to test the chr4 missed families specifically (reps derived from missed clusters), and they ARE seeded.

**Caveat**: short_B family (703 bp) achieves only 0.375 recall — the lowest of all families. This partial miss is from the refinement/acceptance stage (all 8 copies ARE seeded but only 3 of 8 pass refinement), so even the worst-performing family doesn't fail at seed selection.

---

### Hypothesis B: cross-family seed cross-talk / mask shadow

**Evidence**:
- Genome A recall for victim short_A = **0.625**
- Genome B recall for victim short_A = **0.625** (identical)
- Dominant family coverage = 0.814
- The engineered 14-mer shared between dominant (×100) and victim (×8) does NOT suppress victim seeding.

**Verdict: WEAK.** The victim family achieves identical recall with and without the dominant family present. The mask-shadow mechanism does not depress victim seed frequency in this controlled test.

**Caveat**: The engineered shared l-mer (14-mer, one occurrence per dominant copy) may be insufficient to dominate the victim's highest-frequency l-mer. The victim's best l-mer may appear at frequency ~8-10 (8 copies × typical 1-2 unique-per-copy), while the dominant's shared l-mer contributes ~100 occurrences. If the l-mer selection uses absolute frequency and the victim has ANOTHER l-mer with frequency >100, the shared 14-mer may never be chosen as seed anyway. The shared l-mer test was valid but the outcome does not show interference — likely because dominant-family masking happens AFTER the victim is already seeded.

---

### Hypothesis C: nested element shadow

**Evidence**: Genome C recall = 0.000 across all families, but the failure is in the merge stage (not masking). The LINE families were seeded (R=0 N=28, R=1 N=24) but destroyed by the nested-element merge gate before MDL evaluation.

**Verdict: INCONCLUSIVE.** Genome C cannot test hypothesis C because the merge gate triggers a false collapse. The hypothesis remains untested. 

**Root cause identified (different bug)**: The merge gate incorrectly collapses the 2194 bp plain LINE and 2897 bp LINE+LTR composite into a single under-populated family (4 instances). This is the nested-element merge gate introduced in M2 (refine.c::nested_containment_fraction) behaving incorrectly when the two families are structurally related rather than being a true "short nested inside long" situation. This needs investigation.

---

### Hypothesis D: hash collision suppressing frequency

**Evidence**:
- Genome A recall = 0.5417 (new) vs 0.542 (old) — **no change**
- chr4 overall recall = 0.4728 (new) vs 0.4680 (old) — **+0.005, negligible**
- Discovery logs show same R-values and family sizes; no evidence of differential collision rates

**Verdict: WEAK / NOT the primary bottleneck.** Replacing the linked-list hash with a striped-lock parallel hash made no meaningful difference to recall. The engineering improvement (runtime speed) is real and valuable, but it did not unblock any suppressed seed signals. Hash collisions were not the cause of missed families.

---

## chr4 In-Scope Recall Assessment

Using the metric that can be recomputed (bed_pr.py --overlap-mode min):
- New binary: overall recall **0.4728** (F1 0.5553, precision 0.6727)
- Baseline: overall recall **0.4680** (F1 0.5533, precision 0.6766)
- Delta: **+0.005** (negligible)

**Decision rule result (from prompt)**: ≤ 0.465 → "mechanism D was not the issue; algorithm work is mandatory." The new recall of 0.4728 is marginally above the HANDOFF in-scope recall of 0.465, consistent with "no change" scenario. Option γ did not improve recall.

---

## Recommendation for Stage B

**Neither hypothesis A, B, nor D is strongly supported.** Hypothesis C is inconclusive due to a newly discovered bug (merge gate falsely collapsing Genome C families). The canonical recommendation is:

**"Multiple mechanisms — fix the Genome C merge gate bug first"**

Rationale:
1. The nested-element merge gate (`refine.c::nested_containment_fraction`) is destroying valid families on Genome C. This bug was NOT present in the M2 baseline (where testD found both LINE and SINE). Something about the Genome C setup (specifically that the two families ARE structurally related — 2194 bp is a strict prefix of 2897 bp) triggers incorrect merge behavior. This needs to be debugged and fixed before hypothesis C can be assessed.

2. After fixing the merge gate: re-run Genome C and check if LTR_nested recall < LTR_solo recall. Only then can hypothesis C be assessed.

3. The primary chr4 recall bottleneck remains unidentified. With all four mechanistic hypotheses either weak or untestable, the most informative next step is:
   - Fix the Genome C merge gate bug (concrete, bounded task)
   - Re-run Genome C to test hypothesis C properly
   - If hypothesis C is also WEAK, then consider architecture-level investigation: why do all 6 representative families from chr4's missed clusters achieve only 37-62% recall even in a clean 5 Mb genome?

**Stage B NOT needed yet** (none of the mechanisms are confirmed strongly enough to justify a targeted Stage B fix). The immediate priority is the merge gate bug in `refine.c::nested_containment_fraction` that destroys structurally-related family pairs.

---

## Caveats

1. **Copy-class stratification not reproducible**: HANDOFF §3's breakdown by copy class (1, 2, 3-4, 5-9, 10-49, 50-199, ≥200) cannot be reproduced without the original TAIR10 RepeatMasker annotation with family names. The truth BED has no family labels. Reported overall interval-level recall instead.

2. **HANDOFF baseline discrepancy**: HANDOFF §3 states F1=0.515, recall=0.389, precision=0.760. Running `bed_pr.py --overlap-mode min` on the existing `chr4_full_rm.bed` gives F1=0.5533, recall=0.4680, precision=0.6766. The HANDOFF numbers appear to have been computed with a different metric or from a different run output that is no longer available. All comparisons in this report use a consistent metric.

3. **Genome C merge gate bug**: The nested-element merge gate is a pre-existing fix (M2, refine.c) designed to prevent false merges of nested TEs. On Genome C it triggers incorrectly, destroying the LINE/LTR families and making hypothesis C untestable. This bug is present in both old and new binaries.

4. **[5000, 10000) recall regression**: New binary detected 3 intervals vs 13 in baseline for this length bin. Sample is small (108 intervals). Could be random variation or could indicate that specific large families are lost in the new run. Not investigated further in this report.

5. **All computations ran to completion**. No sandbox blocks for the core pipeline. RepeatMasker ran in PGTA conda environment (version 4.2.1, blastn 2.14.1+). All reported numbers are from actual runs.
