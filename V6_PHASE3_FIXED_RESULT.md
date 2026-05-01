# V6 Phase 3 — Fixed Binary Evaluation Report
**Date:** 2026-04-28  
**Binary:** `/home/shuoc/tool/mdl-repeat/bin/mdl-repeat` (rebuilt from working tree, 2026-04-28)  
**Two regression fixes applied to working tree (not yet committed to git):**
1. `ALIGN_MAX_EXTENSION` reverted from 20,000 → 10,000 in `src/align.c` line 21
2. `remap_instance_coordinates()` added per-chunk in `discover_chunked` worker thread (`src/main.c` ~line 492–493)

**Previous broken result:** `/home/shuoc/tool/mdl-repeat/V6_PHASE3_RESULT.md` (not overwritten — evidence of regressions)

---

## EXECUTIVE SUMMARY

| Test | v6 broken (#18) | v6_fixed (this run) | Status |
|------|----------------|---------------------|--------|
| BIO-N2 ATHILA (FL 12kb + soloLTR 600bp) | FAIL — FL missing, only soloLTR found | FAIL — same identical output | **STILL FAILING** |
| chr4 BED negative coords | PASS (single chr, no chunking) | PASS | No change |
| chr4 family recall 80×50 | 120/147 = 0.8163 | 120/147 = 0.8163 | Unchanged |
| chr4 family recall 80×80 | 112/147 = 0.7619 | 112/147 = 0.7619 | Unchanged |
| chr4 family recall 90×80 | 87/147 = 0.5918 | 87/147 = 0.5918 | Unchanged |
| TAIR10 BED negative coords | 5904 negative rows | **1 negative row** (pre-existing edge case) | **MOSTLY FIXED** |
| TAIR10 no-cross-chr boundary | — | PASS (0 over-boundary rows) | PASS |
| TAIR10 family recall 80×50 | 1228/1517 = 0.8095 | **1293/1517 = 0.8523** | **+65 families (+5.3%)** |
| TAIR10 family recall 80×80 | 1154/1517 = 0.7607 | **1234/1517 = 0.8134** | **+80 families (+6.9%)** |
| TAIR10 family recall 90×80 | 908/1517 = 0.5985 | **999/1517 = 0.6585** | **+91 families (+9.8%)** |

**Coordinate fix outcome:** 5903 of 5904 bad TAIR10 BED rows are fixed. The 1 remaining negative row (`1 -2857 12629`) is a pre-existing edge case: an alignment that extends into the 5000-bp PAD region before chr1's sequence start. This is unrelated to the cross-chromosome remap bug and was present as `unknown -2857 12629` even in the committed HEAD.

**TAIR10 recall surge:** The coordinate fix also restored 16,736 correctly-mapped families (vs 11,401 in v6 broken). Many families were previously under-represented in recall because their instances had nonsense coordinates → affected MDL scoring and instance counting → caused those families to fall below thresholds.

**BIO-N2 unresolved:** The two claimed fixes do NOT address the BIO-N2 failure mode. The FL ATHILA family remains missing. Root cause is a destructive merge-chimera bug in the v6 extension algorithm (separate from ALIGN_MAX_EXTENSION value). Details below.

---

## 1. BIO-N2 ATHILA Synthetic — STILL FAILING

### Result: unchanged from v6 broken

| Family | Length | Copies | MDL score | Status |
|--------|--------|--------|-----------|--------|
| R=0    | 600 bp | 64     | 71,782    | soloLTR — found correctly |
| R=3    | 33 bp  | 3      | 5.7       | junk |
| R=8    | 34 bp  | 2      | 5.5       | junk |
| R=4    | 35 bp  | 3      | 4.2       | junk |
| R=6    | 31 bp  | 3      | 3.1       | junk |
| R=2    | 30 bp  | 3      | 0.6       | junk |

Post-MDL stats (including rejected):
```
family_id  consensus_length  num_instances  mdl_score
1          17201             1              0.0   ← FL chimera, rejected
```

**Expected:** Two families — 12,200bp FL ATHILA (N≈20) AND 600bp soloLTR (N≈80)

### Why ALIGN_MAX_EXTENSION revert did NOT fix BIO-N2

The chimera grows 12,201bp → 17,201bp = +5,000bp, which is within the 10,000 cap (per-direction per-iteration). The cap value is irrelevant at this scale.

### Confirmed: BIO-N2 failure is a v6 working-tree regression

Cross-check against committed HEAD (git stash):
- HEAD binary (committed, no v6 changes): finds BOTH 12,201bp FL (N=85, MDL=112,689) AND 600bp soloLTR (N=76, MDL=71,782) — **PASS**
- v6 working tree: same destructive merge produces 17,201bp chimera — **FAIL**

This confirms the BIO-N2 failure was introduced by the v6 working-tree changes and is NOT fixed by the ALIGN_MAX_EXTENSION revert.

### Root cause (confirmed)

1. Discovery correctly finds R=1 (12,201bp, N=40) and R=2 (12,201bp, N=14) — both representing FL ATHILA seeded from different internal-region l-mers.
2. Merge stage: length ratio = 1.0, passes 80-80-80. `estimate_merge_score` approves. Union-find merges them.
3. Re-refinement (`align_refine_family`) on the merged family (N=27) calls the v6 `extend_direction()` function which uses a new `ext_compute_score` / `ExtDP` algorithm (completely rewritten from HEAD).
4. The new extension algorithm does not terminate correctly at the FL element boundary. Consensus grows +5,000bp per iteration → 17,201bp chimera.
5. Only 1 genome position matches the 17,201bp chimera (the specific FL copy whose flanking sequence was incorporated). MDL score = 0.0, family rejected.
6. The original R=1, R=2 FL families are gone (union-find consumed them) — cannot be recovered.

**The bug is in the v6 `extend_direction()` rewrite, not in the ALIGN_MAX_EXTENSION constant.**

### Recommended fix (not implemented)

In `align_refine_family()` (or immediately after re-refinement of a merged family), cap the final consensus length to `max(len_a, len_b) * 1.5`. For two 12,201bp families, this would cap at 18,302bp, which would not prevent the chimera in this case (17,201 < 18,302). A tighter cap of `max(len_a, len_b) * 1.2` = 14,641bp would reject the 17,201bp result. Alternatively, fix the stopping criterion in the new `extend_direction` implementation to match the HEAD behavior.

---

## 2. chr4 v6_fixed vs Baselines — Family-Level Recall

**Run command:** `bin/mdl-repeat -sequence chr4.fa -output chr4_v6_fixed.fa -threads 4 -v`  
**Wall clock:** 7m 48s  
**Output:** 2,912 families, 12,480 instances

### chr4 family-level recall (all 3 criteria, truth clusters: 5,458 total, 147 in-scope ≥3)

| Criterion | v2 (Stage B) | K (Stage B+K) | v6 broken | **v6_fixed** | Δ fixed vs v2 |
|-----------|-------------|---------------|-----------|--------------|---------------|
| **80×50** IN-SCOPE | 123/147 = 0.8367 | 121/147 = 0.8231 | 120/147 = 0.8163 | **120/147 = 0.8163** | **−3** |
| **80×80** IN-SCOPE | 114/147 = 0.7755 | 113/147 = 0.7687 | 112/147 = 0.7619 | **112/147 = 0.7619** | **−2** |
| **90×80** IN-SCOPE | 89/147 = 0.6054  | 88/147 = 0.5986  | 87/147 = 0.5918  | **87/147 = 0.5918**  | **−2** |

**chr4 result: no change.** The two fixes (ALIGN_MAX_EXTENSION revert, chunked remap) do not affect single-chromosome runs. chr4 recall is identical between v6 broken and v6_fixed — as expected, since chr4 uses neither chunking nor the v6 extension algorithm in a way that differs from v6 broken for its specific elements.

### chr4 multi-chr boundary check: PASS (single chromosome, trivially clean)
- All 12,480 BED rows are on `chr4`, no negative coordinates.

---

## 3. TAIR10 v6_fixed vs Baselines — Family-Level Recall

**Run command:** `bin/mdl-repeat -sequence tair10_nuclear.fa -output tair10_v6_fixed.fa -threads 4 -chunk-size 30 -v`  
**Wall clock:** 25m 38s  
**Output:** 16,736 families (+5,335 vs v6 broken), 178,120 instances (+10,534 vs v6 broken)

### TAIR10 family-level recall (truth clusters: 31,206 total, 1,517 in-scope ≥3)

| Criterion | v2 (Stage B) | v6 broken | **v6_fixed** | Δ fixed vs v2 | Δ fixed vs v6 |
|-----------|-------------|-----------|--------------|---------------|---------------|
| **80×50** IN-SCOPE | 1240/1517 = 0.8174 | 1228/1517 = 0.8095 | **1293/1517 = 0.8523** | **+53 (+4.3%)** | **+65 (+5.3%)** |
| **80×80** IN-SCOPE | 1180/1517 = 0.7779 | 1154/1517 = 0.7607 | **1234/1517 = 0.8134** | **+54 (+4.6%)** | **+80 (+6.9%)** |
| **90×80** IN-SCOPE | 924/1517 = 0.6091  | 908/1517 = 0.5985  | **999/1517 = 0.6585**  | **+75 (+7.8%)** | **+91 (+9.8%)** |

**TAIR10 result: substantial improvement** from the coordinate remap fix. v6_fixed is the best result on TAIR10 across all criteria, beating both v2 and v6 broken by large margins.

### Why the remap fix improves recall so dramatically

The v6 broken run had 5,904 instances with invalid coordinates on chr2. When `align_collect_instances` or downstream refinement scanned these corrupted positions, it would:
- Fail to find valid genome sequence → miss re-refinement opportunities
- Score instances incorrectly → MDL underestimates savings → fewer families pass MDL
- Corrupt BED output but also propagate bad coordinates through consensus rebuild

With correct coordinates, all 4 chunks can properly refine their families. The additional 5,335 families represent real repeat families from chr2 (and cross-chr interactions) that were corrupted/lost in the broken run.

### TAIR10 per-class breakdown at 80×50

| Class | v2 covered/total | v6 broken covered | v6_fixed covered |
|-------|-----------------|------------------|-----------------|
| 1     | —/26,760        | 4845             | **5584**        |
| 2     | —/2,929         | 1531             | **1742**        |
| 3-9   | —/1,430         | 1144/1430 = 0.800 | **1208/1430 = 0.845** |
| 10-99 | —/87            | 84/87 = 0.966    | **85/87 = 0.977** |
| 100+  | 0/0             | 0/0              | 0/0             |

---

## 4. Multi-Chr Boundary Check — TAIR10 v6_fixed

### Negative-start coordinates

**v6 broken:** 5,904 rows with negative start (3.5% of all 167,586 rows), all on chr2  
**v6_fixed:** **1 row** with negative start (0.0006% of 178,120 rows), on chr1

The remaining row: `1  -2857  12629  R=52  222  -`

This is NOT a cross-chromosome remap bug. Root cause: a banded alignment hit extends into the 5000-bp PADLENGTH buffer before the first nucleotide of chr1. The same row appeared as `unknown  -2857  12629` in the committed HEAD (before any v6 changes), confirming it is a pre-existing boundary-alignment edge case unrelated to the chunked remap fix.

The `PADLENGTH` check in `align_collect_instances` (`if (ai.genome_start < PADLENGTH) continue;`) should have filtered this, but the banded alignment can report a start before PADLENGTH if the anchor is near position 5000. This is a very minor issue affecting 1 row out of 178,120.

### Over-chromosome-boundary check: PASS
Zero rows where `end > chr_length` for any chromosome. No instance spans a chromosome boundary.

---

## 5. Side-by-Side Comparison Table

### chr4 (18.6 Mb, single chromosome)

| Metric | v2 | K | v6 broken | v6_fixed | Comment |
|--------|----|----|-----------|----------|---------|
| Families output | 2,932 | 2,920 | 2,911 | 2,912 | Stable |
| Instances output | — | — | 12,200 | 12,480 | +280 |
| 80×50 in-scope | **123/147** | 121/147 | 120/147 | 120/147 | −3 vs v2 |
| 80×80 in-scope | **114/147** | 113/147 | 112/147 | 112/147 | −2 vs v2 |
| 90×80 in-scope | **89/147** | 88/147 | 87/147 | 87/147 | −2 vs v2 |
| Negative BED rows | 0 | 0 | 0 | 0 | No chunking on chr4 |
| Wall clock | ~10m | ~10m | — | 7m48s | Parallel hash speedup |

### TAIR10 nuclear (119 Mb, 5 chromosomes, 4 chunks)

| Metric | v2 | v6 broken | v6_fixed | Comment |
|--------|-----|-----------|----------|---------|
| Families output | 11,462 | 11,401 | **16,736** | +5,335 from remap fix |
| Instances output | — | 167,586 | **178,120** | +10,534 |
| Negative BED rows | 0 | **5,904** | **1** | Residual: boundary-alignment edge case |
| Over-chr-boundary rows | 0 | ? | **0** | Confirmed PASS |
| 80×50 in-scope | 1240/1517 | 1228/1517 | **1293/1517** | +53 vs v2, +65 vs v6 |
| 80×80 in-scope | 1180/1517 | 1154/1517 | **1234/1517** | +54 vs v2, +80 vs v6 |
| 90×80 in-scope | 924/1517 | 908/1517 | **999/1517** | +75 vs v2, +91 vs v6 |
| Wall clock | ~40m | ~33m | 25m38s | 4 parallel chunks |

---

## 6. BIO-N2 Diagnosis: Which v6 Change Broke It

**Hypothesis tested:** ALIGN_MAX_EXTENSION 20k → 10k revert would fix BIO-N2.  
**Result:** NO. The revert is irrelevant because the chimera grows by only 5,000bp (within the 10k cap).

**Finding:** The BIO-N2 failure is caused by a behavioral change in `extend_direction()` in the v6 working tree. The v6 version completely rewrote this function from the HEAD version (which uses the RepeatScout-style `compute_totalbestscore_right`/left). The new `ext_compute_score`/`ExtDP` algorithm permits larger consensus extensions when given 27 instances (which have high aggregate score), producing a 17,201bp chimera instead of stopping at 12,201bp.

**The ALIGN_MAX_EXTENSION revert is correct** (prevents future overshoots beyond 10k), but it doesn't fix the current ~5k overshoot in BIO-N2.

---

## 7. Files Produced

| File | Description |
|------|-------------|
| `/tmp/ath_bench/v6_phase3/bio_n2_fixed/athila_lib.fa` | v6_fixed BIO-N2 library (identical to v6 broken: only soloLTR) |
| `/tmp/ath_bench/v6_phase3/bio_n2_fixed/bio_n2.log` | v6_fixed BIO-N2 run log |
| `/tmp/ath_bench/v6_phase3/bio_n2_fixed/bio_n2_stats.tsv` | Per-family stats (17201bp chimera rejected) |
| `/tmp/ath_bench/v6_phase3/bio_n2_fixed/athila_lib_head.fa` | HEAD (committed) BIO-N2 library — BOTH FL and soloLTR found |
| `/tmp/ath_bench/v6_phase3/bio_n2_fixed/bio_n2_head.log` | HEAD BIO-N2 run log (shows correct recovery) |
| `/tmp/ath_bench/v6_phase3/chr4_v6_fixed.fa` | chr4 v6_fixed library (2,912 families) |
| `/tmp/ath_bench/v6_phase3/chr4_v6_fixed.bed` | chr4 v6_fixed instances (12,480, all valid) |
| `/tmp/ath_bench/v6_phase3/chr4_v6_fixed.log` | chr4 v6_fixed run log |
| `/tmp/ath_bench/v6_phase3/chr4_v6_fixed.tsv` | chr4 v6_fixed per-family stats |
| `/tmp/ath_bench/v6_phase3/tair10_v6_fixed.fa` | TAIR10 v6_fixed library (16,736 families) |
| `/tmp/ath_bench/v6_phase3/tair10_v6_fixed.bed` | TAIR10 v6_fixed instances (178,120 rows; 1 negative-start edge case) |
| `/tmp/ath_bench/v6_phase3/tair10_v6_fixed.log` | TAIR10 v6_fixed run log |
| `/tmp/ath_bench/v6_phase3/tair10_v6_fixed.tsv` | TAIR10 v6_fixed per-family stats |

---

## 8. Next Steps Before Phase 4

### Blocker 1: BIO-N2 still failing (FL ATHILA not recovered)

**Required fix:** In the v6 `extend_direction()` rewrite, identify why the stopping criterion permits 5,000bp overshoot past the true element boundary when N=27. Either:
- (a) Add a post-merge consensus length cap: `if len > max(len_a, len_b) × 1.2: revert to pre-merge`
- (b) Fix the `ext_compute_score` / WHEN_TO_STOP termination to match HEAD behavior

Test: `bio_n2_fixed/athila_lib.fa` must contain a family with length ≥ 10,000bp after fix.

### Blocker 2: 1 residual negative BED row on chr1 (minor)

The row `1  -2857  12629` is a boundary-alignment hitting the PAD region before chr1. Fix: strengthen the PADLENGTH check in `align_collect_instances` to reject any hit where `genome_start < PADLENGTH + k` (ensure the alignment genuinely starts in the genome content, not in padding). This is a minor issue but should be fixed before claiming "clean" BED output.

### Non-blocker: chr4 recall -3 vs v2 baseline

The 3 families lost relative to v2 (all 6-copy, 562-1916bp) appear to be a consequence of v6 algorithm changes unrelated to the two fixes. Use `-trace-dir` to identify whether they are being merged into other families or pruned by MDL.

### Can we proceed to Phase 4?

**TAIR10 coordinate fix works** (5903/5904 bad rows fixed, family recall improved dramatically). The TAIR10 output is now suitable for downstream use.

**BIO-N2 must be fixed first** if Phase 4 requires correct handling of large LTR retrotransposons (>10kb). The v6 extension algorithm has a confirmed bug for elements of that size class.
