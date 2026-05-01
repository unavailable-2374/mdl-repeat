# Refine Pipeline Trace Report (chr4)

**Date**: 2026-04-28
**Binary**: `bin/mdl-repeat` with -trace-dir instrumentation
**Genome**: /tmp/ath_bench/chr4.fa (18.6 Mb)
**Truth**: /tmp/ath_bench/chr4_full_truth.bed (12,886 intervals, 6,113,883 bp)

## Stage-by-stage recall curve

| Stage | n_fam | n_pred | cov_bp | rec_loose | rec_strict | bp_rec | Δ bp_rec |
|-------|------:|-------:|-------:|----------:|-----------:|-------:|---------:|
| 01 discover  | 4028 | 16099 | 14,650,314 | 0.4387 | 0.3075 | **0.5802** | — |
| 02 compact   | 4027 | 16098 | 14,650,282 | 0.4386 | 0.3075 | 0.5801 | −0.0001 |
| 03 merge     | 3868 | 17228 | 12,292,823 | 0.3963 | 0.2363 | 0.5204 | **−0.0597** |
| 04 split     | 3868 | 17228 | 12,292,823 | 0.3963 | 0.2363 | 0.5204 | 0 |
| 05 assemble  | 3849 | 17228 | 12,292,823 | 0.3963 | 0.2363 | 0.5204 | 0 |
| 06 mdl_select| 2774 | 12169 |  6,309,462 | 0.3030 | 0.1350 | **0.2798** | **−0.2406** |
| 07 prune     | 2769 | 11983 |  2,826,300 | 0.3020 | 0.1338 | 0.2388 | −0.0410 |
| 08 coalesce  | 2769 | 11086 |  7,263,290 | 0.4146 | 0.2906 | **0.3970** | +0.1582 |

(coalesce is a reporting transform — it merges adjacent same-family same-strand instances into longer intervals, raising recall under the strict ≥50%-overlap metric)

## Diagnosis

**Discovery starts strong**: 0.580 bp_rec at stage 01. The instances exist; they just get destroyed downstream.

### Bottleneck #1: MDL select (−0.241 bp_rec)
- Families: 3849 → 2774 (1075 dropped, 28%)
- Instances: 17228 → 12169 (5059 dropped, 29%)
- Covered bp: 12.3M → 6.3M (49% lost)
- Mechanism: `mdl_select_library` uses unique-coverage greedy. Families whose instances are largely overlapped by already-accepted larger families have near-zero marginal savings → rejected. Even when their consensuses are valid, they get pruned because their territory is "stolen" by the greedy first-comer.

### Bottleneck #2: merge (−0.060 bp_rec)
- Families: 4027 → 3868 (155 merged)
- Instances actually went UP: 16098 → 17228 (+1130, the merge re-aligns and finds extras)
- Covered bp: 14.65M → 12.29M (2.36M LOST)
- Mechanism: when family A (len 500) and family B (len 700) merge, the consensus becomes ~600bp; instances that were 500 or 700 now align to 600 → ends trimmed. Net coverage loss despite more instance count.

### Bottleneck #3: prune (−0.041 bp_rec)
- Acting on already-accepted MDL families; removes 5 marginal ones.
- Cov_bp drop 6.3M → 2.8M is largely because the ".bed" dump filters by mdl_score>0 and prune removes entries that no longer satisfy after exclusion-test. **The covered_bp metric here is misleading**; real instance count only changed 12169 → 11983.

## Stage B target

Two changes, in priority order:

### Fix 1 — Loosen MDL select acceptance (high impact, ~25% bp_rec recoverable)
Current: accept iff `marginal_savings > 0` after greedy unique-coverage.
Proposed: accept if EITHER
  - marginal savings > 0, OR
  - standalone savings > 0 AND consensus_length ≥ 50 AND num_instances ≥ 3

Rationale: if a family's standalone MDL is positive, it's a valid TE family in its own right; the greedy unique-coverage shouldn't kill it just because a neighbor was accepted first. The min copy count + min length still gates spurious noise.

### Fix 2 — Tighten merge length-similarity gate (moderate impact, ~6% bp_rec recoverable)
Current: 80-80-80 (80% identity + 80% mutual length coverage).
Add: only merge when `min(len_A, len_B) / max(len_A, len_B) ≥ 0.7` to avoid merging vastly different-sized families that lose coverage at the consensus step.

Rationale: 80% mutual coverage allows e.g. merging 500bp into 700bp (71% mutual), which then loses 200bp at the consensus. A 70% length-ratio guard rejects this case.

## Caveats

- These bp_rec numbers are computed against the trace BED only (no RM-remap). RM-remap with the final library typically inflates recall by 0.05–0.08 because RM finds additional instances. So the absolute "0.580 → 0.397" trajectory understates final library quality.
- The `ge` -based MDL filter at stage 06 displays only mdl_score>0 in BED dumps; the underlying CandidateList retains all families until prune actually removes them.
- Coalesce stage is a reporting transform; it doesn't affect library quality, only how instances are reported.

## Files

- /tmp/ath_bench/trace/01_discover.{fa,bed,tsv} ... 08_coalesce.{fa,bed,tsv}
- /tmp/ath_bench/trace/run.log
- /tmp/ath_bench/trace/stage_recall.py (analysis script)
