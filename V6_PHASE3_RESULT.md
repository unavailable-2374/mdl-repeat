# V6 Phase 3 Evaluation Report
**Date:** 2026-04-29  
**Binary:** `/home/shuoc/tool/mdl-repeat/bin/mdl-repeat` (built 2026-04-29 17:36)  
**v6 feature set:** Stage A parallel hash + Stage B MDL fallback + Option K banded DP + ENG-N3 dynamic HASH_SIZE + ENG-N7 `-coalesce-factor` CLI + J' ALIGN_MAX_EXTENSION 20k + B+Q6+ENG-N8/N9/N11 cross-chr + ENG-N2/N10/N12 sweep-line memory

---

## EXECUTIVE SUMMARY

**Multi-chr boundary check: FAIL** — 5,904 BED rows have negative coordinates in TAIR10 v6 output.  
Root cause: `discover_chunked()` returns chunk-local instance coordinates without remapping to genome-global. The `remap_instance_coordinates()` function exists but is only called on the genome-sampling code path (line 1197), not after `discover_chunked`. This is an unfixed ENG-N8 regression in the current binary.

**BIO-N2 ATHILA synthetic: FAIL** — The full-length (12 kb) ATHILA family is missing from the final library. The soloLTR (600 bp) family is correctly recovered. Root cause is a destructive merge-collapse (Bug B regression), not an MDL selection failure.

**chr4 family-level recall: v6 is marginally worse than v2 baseline** (−1 to −3 families across all criteria). No improvement from any v6 patch on chr4.

**TAIR10 family-level recall: v6 is marginally worse than v2 baseline** (−1 to −6 families across all criteria). Library BLAST-based recall is unaffected by the BED coordinate bug.

---

## 1. BIO-N2 ATHILA Synthetic Test

### Setup
- Background: 5 Mb Markov-2 (chr4 dinucleotide composition, seed 42)
- 20 full-length ATHILA copies: `LTR(600bp) + internal(11,000bp) + LTR(600bp) = 12,200bp`, 5% divergence
- 80 solo LTR copies: 600bp each, 5% divergence
- LTR and internal use distinct uniform-random sequences (0 shared 16-mers; verified)
- Truth BED: 20 × ATHILA_FL (12,200bp) + 80 × ATHILA_soloLTR (600bp)
- Run: `bin/mdl-repeat -sequence athila_synth.fa -output athila_lib.fa -instances athila_inst.bed -threads 4 -v`

### Result: **FAIL — FL family absent from final library**

**Final library (6 families):**
| Family | Length | Copies | MDL score |
|--------|--------|--------|-----------|
| R=0    | 600 bp | 64     | 71,782    |
| R=3    | 33 bp  | 3      | 5.7       |
| R=8    | 34 bp  | 2      | 5.5       |
| R=4    | 35 bp  | 3      | 4.2       |
| R=6    | 31 bp  | 3      | 3.1       |
| R=2    | 30 bp  | 3      | 0.6       |

- `R=0` (600bp, 64 copies) = soloLTR family → **found correctly**
- Full-length 12,200bp family → **MISSING**

**Per-family stats after MDL (9 families, pre-prune):**
```
family_id  consensus_length  num_instances  mdl_score
0          600               64             71782.4
3          33                3              5.7
8          34                2              5.5
4          35                3              4.2
6          31                3              3.1
2          30                3              0.6
7          31                3              0.0
5          31                3              0.0
1          17201             1              0.0   ← destroyed FL
```

### Root Cause: Bug B Regression (Destructive Merge)

The pipeline discovered both FL families correctly:
- `R=1`: 12,201bp, 40 copies (seeded from internal region)
- `R=2`: 12,201bp, 14 copies (seeded from different internal-region position)

During the merge stage, R=1 and R=2 were merged (same length → ratio 1.0, passes 80-80-80 identity). Re-refinement after merge produced a **17,201bp chimeric consensus** extending ~5,000bp beyond the FL element boundary (ALIGN_MAX_EXTENSION = 20,000bp overshoots). Only 1 of the 54 original copies aligns to this chimera. MDL correctly rejects the chimera (score = 0.0, model_cost = 34,424 bits) but **the original 40-copy and 14-copy FL families were already consumed by union-find** and cannot be recovered.

**This is NOT an MDL select failure.** The Stage B fallback fires correctly (MDL rejects the chimera). The failure is upstream: the estimate_merge_score gate does not foresee that the post-merge re-refinement will generate a chimera of incorrect length. The nested_containment_fraction gate did not fire because both families have equal length (neither is ≥3× the other).

### Distinction from MDL fallback
- MDL score of the 17,201bp chimera = 0.0 → MDL select correctly rejects it
- But the original FL families (R=1, R=2) are gone before MDL select runs
- Stage B (MDL fallback) cannot resurrect union-find consumed families

**Files saved to:** `/tmp/ath_bench/v6_phase3/bio_n2/`

---

## 2. chr4 v6 vs Baselines — Family-Level Recall

**Run:** `bin/mdl-repeat -sequence /tmp/ath_bench/chr4.fa -output chr4_v6.fa -threads 4 -v`  
**Output:** 2,911 families, 12,200 instances  

### Family-level recall table (chr4, truth clusters: 5,458 total, 147 in-scope ≥3)

| Criterion | v2 (Stage B) | K (Stage B+K) | v6 (current) | Δ v6 vs v2 |
|-----------|-------------|---------------|--------------|------------|
| **80×50** IN-SCOPE | **123/147 = 0.8367** | 121/147 = 0.8231 | 120/147 = 0.8163 | **−3** |
| **80×80** IN-SCOPE | **114/147 = 0.7755** | 113/147 = 0.7687 | 112/147 = 0.7619 | **−2** |
| **90×80** IN-SCOPE | **89/147 = 0.6054**  | 88/147 = 0.5986  | 87/147 = 0.5918  | **−2** |

### Per-class breakdown at 80×50

| Class | v2 total | v2 covered | v6 covered | Δ |
|-------|----------|-----------|-----------|---|
| 3-9   | 142      | 119       | 116       | −3 |
| 10-99 | 5        | 4         | 4         | 0 |
| 100+  | 0        | 0         | 0         | — |
| Overall | 5458   | 1287      | 1258      | −29 |

### Note on corrected denominator
The task references a "corrected denominator" of 129 families (from #14 Phase 1 report, 80×50 = 86.8% = 112/129). This denominator cannot be reproduced from the current `family_eval/truth_clusters.fa.clstr` using any standard length or copy-count filter. No threshold between 50bp and 300bp gives exactly 129 clusters. The raw in-scope denominator of **147** (size ≥ 3) is used throughout this report.

### chr4 multi-chr boundary check: PASS (single chromosome, trivially clean)
All 12,200 BED rows are on `chr4`, no negative coordinates.

---

## 3. TAIR10 v6 vs Baselines — Family-Level Recall

**Run:** `bin/mdl-repeat -sequence tair10_nuclear.fa -output tair10_v6.fa -threads 4 -chunk-size 30 -v`  
**Output:** 11,401 families, 167,586 instances  
**Wall clock:** ~33 minutes (4 chunks: chr1=30.4M, chr2=27.0M, chr3=23.5M, chr4+5=38.3M)

### Family-level recall table (TAIR10, truth clusters: 31,206 total, 1,517 in-scope ≥3)

| Criterion | v2 (Stage B) | K (Stage B+K) | v6 (current) | Δ v6 vs v2 |
|-----------|-------------|---------------|--------------|------------|
| **80×50** IN-SCOPE | **1234/1517 = 0.8134** | 1230/1517 = 0.8108 | 1228/1517 = 0.8095 | **−6** |
| **80×80** IN-SCOPE | **1162/1517 = 0.7660** | 1163/1517 = 0.7666 | 1154/1517 = 0.7607 | **−8** |
| **90×80** IN-SCOPE | **910/1517 = 0.5999**  | 912/1517 = 0.6012  | 908/1517 = 0.5985  | **−2** |

Note: these recall numbers are computed via BLAST (library.fa → truth cluster reps) and are unaffected by the BED coordinate bug.

---

## 4. Multi-Chr Boundary Check: FAIL

### Status: FAIL — 5,904 BED rows with negative/invalid coordinates in TAIR10 v6

**First bad row:**
```
2	-549252	-451639	R=14910	111	-
```

**Scope of the problem:**
- Total bad rows (start < 0): **5,904** out of 167,586 (3.5%)
- Chromosomes affected: chr2 (5,903 rows), plus 1 "unknown" row
- Families affected: **2,476 families** have at least one bad instance coordinate
- Worst case: coordinates as negative as −10,311,559 (chr2)
- The "unknown" chromosome row (`unknown -2857 12629`) indicates seq_index is invalid (out of bounds → defaults to "unknown")

### Root Cause (confirmed by code inspection)

`remap_instance_coordinates()` at line 416 in `main.c` converts chunk-local positions to genome-global positions. This function is only called at line 1197, on the **genome-sampling code path** (`if (sample_segments)`). The **chunked discovery path** (`discover_chunked`, called at line 1016) returns instance positions in chunk-local coordinates and NEVER calls this remapping function.

When `output_bed()` later uses `genome->boundaries[seq_index-1]` as `chr_offset`, it subtracts a genome-global chromosome offset from a chunk-local position, producing large negative numbers.

**The ENG-N8 cross-chr fix in the current binary is incomplete.** The remap function exists but is wired to the wrong code path.

### Fix required (not implemented in current v6)

In `discover_chunked()` (line 495), after combining all chunk results at line ~726 (before `return combined`), call the remap for each chunk's segments:

```c
/* Remap chunk-local instance coords to genome-global coords */
for (int b = 0; b < num_bins; b++) {
    if (chunk_results[b])
        remap_instance_coordinates(chunk_results[b],
                                   bin_segments[b], bin_counts[b]);
}
```

This must happen before `free(bin_segments)`.

---

## 5. Honest Interpretation

### Did v6 patches help?

**No. v6 is marginally worse than v2 on all family-level criteria.**

| Patch | chr4 impact | TAIR10 impact | Assessment |
|-------|-------------|---------------|------------|
| Stage A parallel hash | Speed improvement (no wall-clock data available) | Chunked TAIR10 in ~33 min | ENG-only, correct |
| Stage B MDL fallback | MDL gate fires correctly in merge (vetoed 3 of 723 pairs on TAIR10) | Same | Works as intended |
| Option K banded DP | −2 families vs v2 (K→v6 unchanged at 80×80) | Marginal | No improvement |
| ENG-N3 dynamic HASH_SIZE | Enabled TAIR10 run to complete | Essential | Correct |
| ENG-N7 -coalesce-factor CLI | No change (default 20× used) | — | No impact |
| J' ALIGN_MAX_EXTENSION 20k | Enabled long-element extension but also enables overshoot chimeras (BIO-N2 bug) | — | Neutral/negative |
| ENG-N8 cross-chr fix | **INCOMPLETE** — remap not called for chunked path | **5,904 bad BED rows** | **Bug** |
| ENG-N2/N10/N12 sweep-line memory | Memory management improvement | No change in recall | Correct |

### chr4 regression (v6 vs v2): −3 families at 80×50

The 3 families missed by v6 but present in v2 (from the missed list at 80×80):
- `chr4:7690987-7692903` (size=6)
- `chr4:17820586-17821147` (size=6)  
- `chr4:7690446-7690969` (size=6)

These are short (562-1916bp), 6-copy families that v2 captured but v6 lost. Most likely explanation: the ALIGN_MAX_EXTENSION=20k change or merge-gate changes in v6 altered which consensuses these families merged into or extended past.

### The BIO-N2 bug is a real design issue

The merge of two same-element FL families (seeded from different positions within the 12kb element) producing a chimeric 17kb consensus is not specific to ATHILA. The same will happen with any element >10kb where multiple seeds are drawn from different sub-regions. The `estimate_merge_score` MDL gate is too optimistic because it cannot predict post-merge re-refinement overshoot. A length ceiling in `align_refine_family` proportional to `max(len_a, len_b) × 1.5` would prevent the chimera.

---

## 6. Files

| File | Description |
|------|-------------|
| `/tmp/ath_bench/v6_phase3/bio_n2/athila_synth.fa` | 5 Mb ATHILA synthetic genome |
| `/tmp/ath_bench/v6_phase3/bio_n2/athila_truth.bed` | Ground truth (20 FL + 80 soloLTR) |
| `/tmp/ath_bench/v6_phase3/bio_n2/athila_lib.fa` | v6 output library (6 families, FL missing) |
| `/tmp/ath_bench/v6_phase3/bio_n2/bio_n2.log` | Run log |
| `/tmp/ath_bench/v6_phase3/bio_n2/bio_n2_stats.tsv` | Per-family MDL stats (shows 17201bp chimera) |
| `/tmp/ath_bench/v6_phase3/chr4_v6.fa` | chr4 v6 library (2,911 families) |
| `/tmp/ath_bench/v6_phase3/chr4_v6.bed` | chr4 v6 instances (12,200, all valid) |
| `/tmp/ath_bench/v6_phase3/chr4_v6.log` | chr4 v6 run log |
| `/tmp/ath_bench/v6_phase3/tair10_v6.fa` | TAIR10 v6 library (11,401 families) |
| `/tmp/ath_bench/v6_phase3/tair10_v6.bed` | TAIR10 v6 instances (167,586 rows, **5,904 with negative coords**) |
| `/tmp/ath_bench/v6_phase3/tair10_v6.log` | TAIR10 v6 run log |

---

## 7. Actionable Next Steps

1. **Fix ENG-N8 (blocker):** Move `remap_instance_coordinates()` call inside `discover_chunked()`, before `return combined`, iterating over each chunk's `bin_segments[b]`. One-time remap per chunk after combining. This fixes the 5,904 bad BED rows on TAIR10.

2. **Fix BIO-N2 merge-chimera bug:** In `align_refine_family()` or wherever re-refinement extends a consensus after merge, cap the final consensus length at `max(len_a, len_b) × 1.5`. This prevents a 12kb element being extended to 17kb by boundary overshoot.

3. **Investigate chr4 regression (v6 vs v2 −3 families):** The 3 lost families (all 6-copy, 562-1916bp) need traceback — check if they were merged into other families or pruned by MDL. Use `-trace-dir` after confirming trace-dir works correctly.

4. **Do not report v6 BED output as validated** until ENG-N8 is fixed. The library `.fa` is valid for BLAST-based recall but the instance `.bed` is incorrect for TAIR10.
