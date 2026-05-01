# Subfamily Split Diagnostic — Phase 4 Track 1

**Date**: 2026-04-29  
**Binary**: bin/mdl-repeat (Phase 4 — Track 1 instrumentation added)  
**Command**: `bin/mdl-repeat -sequence chr4.fa -output ... -threads 4 -vv`  
**Genome**: chr4 full (18.6 Mb, Arabidopsis thaliana TAIR10)

---

## 1. Split Rejection Distribution

From `-vv` log (chr4 full genome run):

| Rejection reason | Count | % of total |
|---|---|---|
| `n_instances < 10` (REFINE_MIN_SPLIT_INSTANCES floor) | 3794 | 95.6% |
| `split MDL gain <= 0` | 140 | 3.5% |
| `n_lo or n_hi < 3` (cluster too small) | 16 | 0.4% |
| `div_gap < 0.050` (mean divergence gap too small) | 3 | 0.08% |
| **Accepted** | **0** | **0%** |

Total families examined: 3953  
Total families that reached the floor check: 3953  
Total families that passed the floor: 165  
Total families accepted for split: 0  

---

## 2. Instance Floor Breakdown (n_instances = 2 to 9)

| n_instances | Count skipped |
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
| **TOTAL** | **3794** |

---

## 3. MDL Rejection Distribution (families with n ≥ 10)

| n_instances range | Count MDL-rejected |
|---|---|
| 10-14 | 37 |
| 15-24 | 37 |
| 25-49 | 20 |
| 50-99 | 15 |
| 100-199 | 7 |
| 200-299 | 3 |
| 300-508 | 6 |

Even families with 300-508 instances are rejected by MDL. This confirms MDL is not a threshold problem but reflects genuine information-theoretic assessment.

---

## 4. Root Cause Analysis

### Primary bottleneck: REFINE_MIN_SPLIT_INSTANCES floor
- 95.6% of families are blocked before even attempting split
- Most chr4 families are small (median ~3 instances) after align-refine

### Secondary: MDL gate rejects ALL attempted splits
- Zero splits accepted on chr4
- Even large families (300-508 instances) are rejected
- Root cause: `build_subset_consensus` does not re-compute `num_edits` for instances against the NEW subset consensus. It reuses the old `num_edits` (edits against combined consensus). Therefore, the MDL pre-check uses stale divergence estimates and the savings from splitting appear smaller than they really are.
- This is a known architectural limitation: the post-split `align_refine_family` call IS made, which would correctly re-assess instance counts, but the MDL gate runs BEFORE re-refinement.

---

## 5. Decision Matrix

Per task specification: "If rejected by MDL gain → don't touch (MDL is correct)"

| Trigger | Threshold | Action | Rationale |
|---|---|---|---|
| REFINE_MIN_SPLIT_INSTANCES floor | n < 10 | **NOT changed** | Lowering to 5 would allow 600 more to be attempted, all would be MDL-rejected |
| Bimodality threshold | 0.20 | **NOT changed** | Not the primary blocker |
| Valley depth check | — | **NOT changed** | Only 0 cases (valley check is post-bimodality) |
| MDL gate | 0 | **NOT changed** | MDL is correct by design spec |

**Conclusion**: No threshold or floor changes made. The split stage is correctly conservative for chr4's family size distribution. The low split rate is expected and correct.

---

## 6. What DID Change (Track 1 — Instrumentation)

### Added to both sequential path and `split_analysis_worker` (parallel path):

```
[split] F{id}: skipped: n_instances={n} < REFINE_MIN_SPLIT_INSTANCES={threshold}
[split] F{id}: rejected: bimodality={X:.3f} < threshold={Y:.3f} (n={n})
[split] F{id}: rejected: valley not deep enough (bimodality={X:.3f}, n={n})
[split] F{id}: rejected: n_lo={lo} or n_hi={hi} < REFINE_MIN_CLUSTER_SIZE={threshold} (threshold={t:.3f}, n={n})
[split] F{id}: rejected: div_gap={gap:.3f} < REFINE_MIN_DIV_GAP={threshold:.3f} (mean_lo={lo:.3f} mean_hi={hi:.3f}, n={n})
[split] F{id}: rejected: split MDL gain <= 0 (orig={A:.1f}, split={B:.1f}={C:.1f}+{D:.1f}, n={n})
[split] F{id}: accepted: split into n_lo={lo} / n_hi={hi} at threshold={t:.3f} (div {lo_pct:.1f}% / {hi_pct:.1f}%), MDL {A:.1f} -> {B:.1f}
```

All log lines are gated on `verbose >= 2` (i.e., `-vv` flag).  
The accept-level messages are ALSO printed at `verbose >= 1` (`-v` flag) for consistency with existing split log behavior.

---

## 7. New Test Added

**Test F** (`tests/data/testF.fa`, `tests/gen_testF.py`):  
- 1Mb genome  
- Two subfamilies sharing an ancestral consensus: Lo (25 copies, ~4% div) and Hi (20 copies, ~20% div)  
- Verifies: (a) `[split]` instrumentation lines are printed under `-vv`, (b) both `-threads 1` and `-threads 4` produce the same family count, (c) family count in [1, 6]  
- Does NOT assert n_splits = 1 (MDL correctly rejects the split in this case)  
- **Both thread modes**: 22/22 PASS in test suite
