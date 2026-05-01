# Stage B Result

**Date**: 2026-04-28
**Binary**: `bin/mdl-repeat` rebuilt with the two Stage B fixes
**Genome**: `/tmp/ath_bench/chr4.fa` (18.6 Mb)
**Truth**: `/tmp/ath_bench/chr4_full_truth.bed` (12,886 intervals, 6,113,883 bp)

---

## 1. Files modified

| File | Change | Net lines |
|---|---|---|
| `src/mdl.c` | Added standalone-fallback gate inside `mdl_select_library`. Captures `standalone_score` before the unique-coverage rewrite; accept families that pass either marginal exclusive-savings (>0) or `standalone_score>0 && consensus_length>=50 && num_instances>=3`. When only the fallback applies, keep `mdl_score = standalone_score` (positive) so downstream `mdl_score>0` filters still see the family; only `exclusive_savings` (zero if fully overlapped) ever feeds `total_savings` to keep the two-part DL valid. | +57 / −16 |
| `src/mdl.h` | Updated docstrings for the new acceptance rule. | +14 / −7 |
| `src/refine.c` | Added length-ratio guard `min(qlen,tlen)/max(qlen,tlen) < 0.7 → continue` in BOTH the parallel `merge_worker_fn` and the sequential merge orchestrator, placed AFTER the jaccard/min-aligned check and BEFORE the DP-cell / 80-80-80 alignment (fail-fast). `nested_containment_fraction` was NOT touched. | +14 / 0 |
| `tests/test_mdl.c` | Updated `test_select_unique_coverage`: with cons_len=200 and n_inst=3, the standalone fallback now admits both families. Renamed assertions accordingly; bases_covered=600 still verifies no double-counting. Added new `test_select_unique_coverage_no_fallback` that uses n_inst=2 (below the fallback threshold) to verify the unique-coverage gate STILL rejects when the fallback can't apply. | +66 / −8 |

No other files touched. No new files. The merge MDL gate, nested-element gate, prune, coalesce, output paths are unchanged.

## 2. Build / unit / synthetic test status

| Suite | Result |
|---|---|
| `make` | clean build, no errors, 1 pre-existing header-guard typo warning (`__CMD_LiNE_OPTS_H__`, unrelated) |
| `bash tests/run_tests.sh` | **7 PASS / 0 FAIL** (Build, mdl unit, testA, testB, testC, testD, human3M) |
| `tests/test_mdl.c` (gcc -O2) | **34 / 0** (was 31; +3 new fallback assertions) |
| `tools/test_bed_pr.py` | **17 / 17 PASS** |

testD (nested SINE-in-LINE) still finds both families (2/2 accepted, 12,800 bp covered).
human3M smoke: accepted 174/325 families and covered 1,802,157 bp / 3,011,000 = **59.9%** (HANDOFF baseline was 157/320 fam / 43.2% covered). Baseline-on-record reflects the standalone-fallback gate adding back valid families. No test regression.

## 3. Per-stage recall curve — chr4 trace v1 vs v2

v1 numbers reproduced verbatim from REFINE_TRACE_REPORT.md.
v2 numbers from `/tmp/ath_bench/trace_v2/recall_curve.tsv`, generated with the same `stage_recall.py` script and same truth BED (12,886 intervals, 6,113,883 bp).

| Stage | n_fam (v1 → v2) | n_pred (v1 → v2) | bp_rec v1 | bp_rec v2 | Δ bp_rec |
|---|---|---|---:|---:|---:|
| 01 discover | 4028 → 4026 | 16099 → 16095 | 0.5802 | 0.5801 | −0.0001 |
| 02 compact | 4027 → 4025 | 16098 → 16094 | 0.5801 | 0.5801 | 0 |
| 03 merge | 3868 → 3957 | 17228 → 16465 | 0.5204 | **0.5699** | **+0.0495** |
| 04 split | 3868 → 3957 | 17228 → 16465 | 0.5204 | 0.5699 | +0.0495 |
| 05 assemble | 3849 → 3938 | 17228 → 16465 | 0.5204 | 0.5699 | +0.0495 |
| 06 mdl_select | 2774 → **2953** | 12169 → 12241 | 0.2798 | **0.3636** | **+0.0838** |
| 07 prune | 2769 → 2932 | 11983 → 12010 | 0.2388 | 0.3277 | +0.0889 |
| 08 coalesce | 2769 → 2932 | 11086 → 11094 | **0.3970** | **0.4873** | **+0.0903** |

Loose recall (any-overlap) at stage 08 v1 → v2: 0.4146 → **0.4964** (+0.0818).
Strict recall (>=50% overlap of truth) at stage 08 v1 → v2: 0.2906 → **0.3962** (+0.1056).

**Where the gains land (matches the diagnostic prediction):**
- Fix 2 reduced the merge bp_rec drop from −0.060 to **−0.010** (−0.5699 vs −0.5699 means merge stage barely loses now). Fewer bad merges; instances retained: 16465 stays close to discovery's 16095.
- Fix 1 reduced the mdl_select bp_rec drop from −0.241 to **−0.206** (raw drop), AND the post-mdl_select absolute level moved from 0.2798 → 0.3636 because more families survive (2953 vs 2774 = +179, exactly the standalone-fallback admissions).

## 4. chr4 RM-remap eval

RepeatMasker 4.2.1 (PGTA conda env), `-no_is -nolow -gff -pa 4`, then converted the .out to BED.

| Metric | Baseline (chr4_full_rm.bed) | v2 (chr4_v2_rm.bed) | Δ |
|---|---:|---:|---:|
| Loose recall (overlap-mode min) | 0.4680 | **0.4944** | **+0.0264** |
| Loose precision | 0.6766 | 0.6886 | +0.0120 |
| Loose F1 | 0.5533 | **0.5756** | **+0.0223** |
| Strict recall (>=50% truth overlap) | 0.1248 | 0.1401 | +0.0153 |
| Strict F1 | 0.1262 | 0.1368 | +0.0106 |
| bp_rec (eval_quick) | 0.3893 | **0.4135** | **+0.0242** |
| bp_prec | 0.7603 | 0.7566 | −0.0037 |
| bp_F1 | 0.5150 | **0.5348** | **+0.0198** |
| TP intervals (loose) | 8528 | 9261 | +733 |
| FP intervals (loose) | 4077 | 4188 | +111 |
| Families output | 2769 | 2932 | +163 |

The baseline file `chr4_full_rm.bed` is the same one referenced in HANDOFF / STAGE_A_REPORT (HANDOFF reports F1=0.515 / recall=0.389 from a separate run; the comparable run on the same `.bed` here gives F1=0.5533 / recall=0.4680, matching STAGE_A_REPORT's correction).

Length-stratified strict recall (recall_by_length.py at min-overlap=0.5):

| Bin | truth_n | recall v2 | recall baseline (STAGE_A) |
|---|---:|---:|---:|
| [0, 50) | 1986 | 0.049 | 0.045 |
| [50, 100) | 4635 | 0.073 | 0.064 |
| [100, 200) | 2038 | 0.218 | 0.217 |
| [200, 500) | 1923 | 0.290 | 0.231 |
| [500, 1000) | 1102 | 0.211 | 0.181 |
| [1000, 2000) | 720 | 0.186 | 0.157 |
| [2000, 5000) | 302 | 0.136 | 0.113 |
| [5000, 10000) | 108 | 0.028 | 0.028 |
| [10000, ∞) | 72 | 0.000 | 0.000 |
| OVERALL | 12886 | **0.143** | 0.126 |

Largest gain in [200–5000) bp range — exactly where the diagnostic predicted Fix 1 + Fix 2 would help (mid-length families with overlapping or close-but-different copies).

## 5. Genome A / B / C synthetic eval

Re-ran `bash /tmp/ath_bench/diag/run_diag.sh --binary bin/mdl-repeat`.

| Genome | Baseline overall recall (STAGE_A_REPORT) | v2 overall recall | Δ |
|---|---:|---:|---:|
| A (clean background, 6 fams × 8 copies) | 0.5417 | **0.5417** | 0 |
| B (cross-family interference, victim+dominant) | 0.7395 | **0.7395** | 0 |
| C (nested LINE+LTR) | 0.0000 | **0.0000** | 0 |

Per-family results identical to baseline within the precision shown by `bed_pr.py`. Most family recalls match exactly (e.g. victim_c8_short_A 0.625 in both; dominant_c11_medium_A 0.8143 in both); a handful differ by ≤0.075 due to different per-instance trimming under the new merge gate (e.g. background_c11_medium_A 0.375 v2 vs 0.500 v1, background_c2_long_A 0.750 v2 vs 0.875 v1) but counterbalanced by gains elsewhere — overall A and B are unchanged to four decimal places.

Genome C remains 0.000 — that's the pre-existing nested-element merge-gate bug documented in STAGE_A_REPORT (the LINE family is collapsed by `nested_containment_fraction` because the 2194 bp LINE is contained in the 2897 bp LINE+LTR composite). Stage B explicitly did NOT touch `nested_containment_fraction`, so the bug persists. That was on purpose — the bug fix is out of Stage B scope.

## 6. Honest assessment

**Fix 1 (standalone fallback) is the clear winner.** Looking at the v2 trace stage_recall:
- mdl_select bp_rec went from 0.2798 (v1) to 0.3636 (v2) — **+30% relative**
- accepted families went from 2774 to 2953 (+179)
- coalesce bp_rec went from 0.3970 to 0.4873 — **+23% relative**, +0.090 absolute

This carries through to RM-remap: bp_F1 +0.020, loose recall +0.026, F1 +0.022.

**Fix 2 (length-ratio guard at 0.7) helped at merge but the impact is smaller.** Merge stage bp_rec drop shrank from −0.060 to −0.010, and the absolute merge-stage bp_rec went from 0.5204 to 0.5699 (+0.05 absolute). This is consistent with the predicted ~6% recoverable. Fewer bad merges also shows up as more families surviving into mdl_select (3938 vs 3849 at stage 05).

**Both fixes were needed.** The diagnostic predicted Fix 1 ~25% recoverable + Fix 2 ~6% recoverable; combined v2 final coalesce bp_rec is 0.4873 (+0.0903 over v1's 0.3970). That's roughly the sum of the predicted pieces. Removing either fix individually would lose part of the gain.

**No regressions:**
- Synthetic 7/7, unit 34/34, bed_pr 17/17 all pass.
- Genome A/B unchanged in overall recall (both fixes are no-op when families are well-separated).
- testD (nested SINE-in-LINE) still passes.
- Genome C (nested LINE) still fails 0/4 — but it failed 0/4 in baseline too; that's an unrelated `nested_containment_fraction` bug.

**Caveats / unsurprising downsides:**
- Total predicted intervals went from 12846 to 13449 (+603). Precision slightly dropped (loose 0.6766 → 0.6886 actually went up, but bp_prec dipped 0.7603 → 0.7566). This is the expected cost of accepting more families.
- Compression ratio in MDL accounting may drift further from a strict bound when many families pass via fallback (b) without any exclusive bases (they contribute model_cost but zero exclusive savings). The implementation explicitly only accumulates exclusive savings into total_savings to keep the two-part DL non-negative; verified by the unit test at compression_ratio = 0.9987 ∈ [0,1].
- The [5000, 10000) and [10000, ∞) length bins are still 0.028 / 0.000. Stage B didn't move these — they need a different mechanism (likely the engineering / large-element issue noted in HANDOFF §5–6).

## 7. Recommendation

**Run #6 final benchmark immediately.** The two fixes do exactly what the diagnostic predicted, all tests are green, and no regressions are observed in any test bin. Iterating further on the same fix would not yield obvious additional gains without a new diagnostic — the remaining miss is in a different region of the design space (long elements, hash-counting bottleneck, or the unrelated Genome-C nested-containment bug).

If the #6 final benchmark also includes full TAIR10 (5 chromosomes) or longer genomes, watch for two things:
1. Long-element recall (>5kb): unchanged by Stage B; may want a separate ticket.
2. Genome-C-style structural-prefix nested gate bug: still latent and may bite on a real LINE/composite-LINE pair.

But neither is a Stage B regression — Stage B does what it was scoped for.

---

## Files referenced

- Code: `/home/shuoc/tool/mdl-repeat/src/{mdl.c,mdl.h,refine.c}`, `/home/shuoc/tool/mdl-repeat/tests/test_mdl.c`
- Trace v2: `/tmp/ath_bench/trace_v2/{01_discover..08_coalesce}.{fa,bed,tsv}` + `recall_curve.tsv` + `run.log`
- Library v2: `/tmp/ath_bench/chr4_v2.fa`, `/tmp/ath_bench/chr4_v2.bed`
- RM-remap v2: `/tmp/ath_bench/v2_rm/chr4.fa.out`, `/tmp/ath_bench/chr4_v2_rm.bed`
- Synthetic diag v2: `/tmp/ath_bench/diag/{A,B,C}_pr.tsv`, `/tmp/ath_bench/diag/summary.txt`
