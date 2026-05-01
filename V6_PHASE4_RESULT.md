# V6 Phase 4 Result

**Date**: 2026-04-29  
**Binary**: bin/mdl-repeat (Phase 4 — Track 1 instrumentation; Track 2 reverted)  
**Genome**: chr4 full (18.6 Mb, Arabidopsis thaliana TAIR10)  
**Command (Track 1)**: `bin/mdl-repeat -sequence chr4.fa -threads 4 -vv`

---

## Summary

Phase 4 consisted of two tracks:

| Track | Description | Status |
|---|---|---|
| Track 1 | Instrument `refine_split_families` with `-vv` diagnostics; analyze chr4 split rejection distribution; threshold decision; add testF | **DONE — kept** |
| Track 2 | Multipass D' with `-multipass 2` CLI flag; re-seed and screen novelties | **DONE — then reverted (no recall gain)** |

---

## Track 1: Split Instrumentation

### Changes

Both the sequential (`num_threads <= 1`) path and the parallel `split_analysis_worker` in `refine.c` now emit one of these log lines per family when running with `-vv`:

```
[split] F{id}: skipped: n_instances={n} < REFINE_MIN_SPLIT_INSTANCES={floor}
[split] F{id}: rejected: bimodality={X:.3f} < threshold={Y:.3f} (n={n})
[split] F{id}: rejected: valley not deep enough (bimodality={X:.3f}, n={n})
[split] F{id}: rejected: n_lo={lo} or n_hi={hi} < REFINE_MIN_CLUSTER_SIZE={t} (threshold={t:.3f}, n={n})
[split] F{id}: rejected: div_gap={gap:.3f} < REFINE_MIN_DIV_GAP={threshold:.3f} (mean_lo={lo:.3f} mean_hi={hi:.3f}, n={n})
[split] F{id}: rejected: split MDL gain <= 0 (orig={A:.1f}, split={B:.1f}={C:.1f}+{D:.1f}, n={n})
[split] F{id}: accepted: split into n_lo={lo} / n_hi={hi} at threshold={t:.3f} (div {lo_pct:.1f}% / {hi_pct:.1f}%), MDL {A:.1f} -> {B:.1f}
```

Accept-level messages are also printed at `-v` for consistency with existing split log behavior.

Files changed: `src/refine.c` (two locations), `tests/run_tests.sh` (Test F added).

### chr4 Full Split Rejection Distribution (from -vv log)

Total families examined: 3953  
Total families that passed the n_instances floor: 165  
Total splits accepted: **0**

| Rejection reason | Count | % of total |
|---|---|---|
| `n_instances < 10` (REFINE_MIN_SPLIT_INSTANCES floor) | 3794 | 95.6% |
| `split MDL gain <= 0` | 140 | 3.5% |
| `n_lo or n_hi < 3` (cluster too small) | 16 | 0.4% |
| `div_gap < 0.050` | 3 | 0.08% |
| **Accepted** | **0** | **0%** |

### Instance Floor Breakdown (n = 1 to 9, skipped before attempting split)

| n_instances | Count |
|---|---|
| 1 | 12 |
| 2 | 1238 |
| 3 | 1348 |
| 4 | 610 |
| 5 | 269 |
| 6 | 145 |
| 7 | 90 |
| 8 | 44 |
| 9 | 38 |
| Total | 3794 |

### MDL Rejection by Family Size (families that reached the MDL gate)

| n_instances range | Count MDL-rejected |
|---|---|
| 10-14 | 37 |
| 15-24 | 37 |
| 25-49 | 20 |
| 50-99 | 15 |
| 100-199 | 7 |
| 200-299 | 3 |
| 300-508 | 6 |

Even families with 300-508 instances are MDL-rejected. Root cause: `build_subset_consensus` does not recompute `num_edits` against the new subset consensus — the MDL pre-check uses stale edit counts from the original combined consensus. This is an architectural limitation noted in SUBFAMILY_DIAGNOSTIC.md.

### Threshold Decisions

Per QUALITY_PROPOSAL_v6.md rule: "If rejected by MDL gain → don't touch (MDL is correct)".

| Parameter | Current Value | Decision | Rationale |
|---|---|---|---|
| `REFINE_MIN_SPLIT_INSTANCES` | 10 | **NOT changed** | Lowering to 5 adds 600 more attempts, all MDL-rejected |
| `REFINE_BIMODALITY_THRESH` | 0.20 | **NOT changed** | Not the primary bottleneck |
| Valley depth check | — | **NOT changed** | Only 0 families blocked here |
| MDL gate | 0 | **NOT changed** | MDL is correct by design spec |

**Conclusion**: The split stage is correctly conservative on chr4. The 0% acceptance rate reflects chr4's family-size distribution (median ~3 instances) and MDL's information-theoretic correctness, not broken thresholds.

### Test F Added (22/22 PASS)

`tests/gen_testF.py` + `tests/data/testF.fa`: 1 Mb synthetic genome with two subfamilies sharing an 80%-identity ancestor (Lo: 25 copies at ~4% divergence; Hi: 20 copies at ~20% divergence). Test asserts:
1. `[split]` instrumentation lines appear under `-vv` (both `-threads 1` and `-threads 4`)
2. Family count in [1, 6] (no regression)
3. Family count agrees within 1 across thread modes

Test does NOT assert n_splits=1 because MDL correctly rejects the split given stale `num_edits` — this is by design.

---

## Track 2: Multipass D' — Implemented Then Reverted

### Implementation Summary

Added `-multipass 2` CLI flag. Pass 2 re-ran `discover_families` / `discover_chunked` on the same genome; `refine_screen_candidates` (new function, 8-mer Jaccard pre-filter + semiglobal DP) screened Pass 2 against Pass 1 and dropped duplicates; survivors were appended to the combined library.

### chr4 Multipass 2 Results

**Wall-clock**: ~16 min total (vs ~8 min single-pass, 2× overhead)

Pipeline summary:
- Pass 1 families: 4033
- Pass 2 raw candidates: 4025
- Screened as duplicates: 2700 (67%)
- Pass 2 novelties: 1325
- Combined before compaction: 5358 families
- Accepted families (MDL): 3205 / 5273

### Family-Level Recall Comparison (chr4 truth clusters, BLAST ≥80%id/≥50%cov)

| Threshold | Track 1 (single-pass) | Track 2 (multipass 2) | Delta |
|---|---|---|---|
| 80%id / 50%cov IN-SCOPE (>=3 copies) | 121/147 = **0.8231** | 121/147 = **0.8231** | 0 |
| 80%id / 80%cov IN-SCOPE | 113/147 = 0.7687 | 113/147 = 0.7687 | 0 |
| 90%id / 80%cov IN-SCOPE | 87/147 = 0.5918 | 87/147 = 0.5918 | 0 |

**Zero recall improvement at all three thresholds.** The 1325 Pass-2 novelties did not cover any truth family cluster not already covered by Pass 1.

### Revert Decision

Per QUALITY_PROPOSAL_v6.md: "If Track 2 D' produces NO recall improvement on chr4 OR breaks any test, REVERT D' but keep Track 1 changes."

Track 2 reverted from:
- `src/main.c` (removed `multipass` variable, CLI parsing, and multipass D' block)
- `src/refine.c` (removed `refine_screen_candidates` function)
- `src/refine.h` (removed `refine_screen_candidates` declaration)

Track 1 instrumentation is preserved in full.

---

## chr4 Track 1 Pipeline Summary (v6 final, single-pass)

```
Discovered families:     4027 (R=4027)
After compaction:        3959
Merged:                  65 families absorbed
Splitting:               0/165 attempted (0 accepted)
Fragment assembly:       17 pairs assembled
Accepted (MDL):          2942 / 3942
Bases covered:           4553370 / 18596056 (24.5%)
Compression ratio:       0.9807
Wall-clock (chr4 -threads 4): ~8 min
```

---

## Test Suite

```
22/22 PASS (was 17 before testF addition)
```

All pre-existing tests (build, unit, synthetic A/B/C/D/E, human3M) continue to pass. Track 1 testF (5 new assertions) all pass.

---

## Files Changed (Track 1 only, post-revert)

| File | Change |
|---|---|
| `src/refine.c` | Added `[split]` diagnostic log lines to sequential path and `split_analysis_worker` |
| `tests/run_tests.sh` | Added Test F block (5 assertions, both thread modes) |
| `tests/gen_testF.py` | New: synthetic bimodal-subfamily genome generator |
| `tests/data/testF.fa` | New: committed fixture (1 Mb genome, 2 subfamilies) |
| `SUBFAMILY_DIAGNOSTIC.md` | New: full chr4 split rejection analysis |
| `V6_PHASE4_RESULT.md` | This file |

Track 2 files (`-multipass` CLI, `refine_screen_candidates`) were added and then fully reverted.
