# Option K — Banded-DP Refine Extension

Replaces column-by-column majority voting in `src/align.c::extend_direction`
with the same banded dynamic-programming scheme used by
`discover.c::extend_right` (RepeatScout port).  Single instances with small
indels no longer mis-vote at every column past the indel — the band absorbs
shifts up to MAXOFFSET.

## 1. Files modified

| File | Change |
|------|--------|
| `src/align.c` | `extend_direction` body rewritten to banded DP. `align_extend_consensus` and the post-extension realloc/shift code unchanged. New helpers added: `ext_genome_base`, `ext_compute_score`, `ExtDP` struct. |

Only `src/align.c` was modified.  No header signatures changed.  No
constants introduced.  Reuses existing `g_align_gap`, `g_align_maxoffset`,
`ALIGN_MATCH`, `ALIGN_MISMATCH`, `ALIGN_CAPPENALTY`, `ALIGN_WHEN_TO_STOP`,
`EXTENSION_SLACK`.

Net diff in `extend_direction` region (lines 651..975 in the new file):
- Old column-vote loop body removed (~100 lines)
- New banded-DP body added (~200 lines, including helpers)
- Post-extension consensus realloc/shift (lines 945..end) unchanged

## 2. Algorithm summary

For each new column `y` past the current consensus edge:

```
for a in {A, C, G, T}:
    total_a = 0
    for n in 0..N-1:                       # contributing instances
        floor = max(0, bestbest[n] + CAPPENALTY)
        bs = floor
        for offset in [-M, +M]:
            bs = max(bs, score(y, n, offset, a))   # 3-case DP
        total_a += bs
besta = argmax(total_a)
commit score[y%2][n][*] for besta
update bestbest[n] and besttotalbestscore -> besty
if y - besty >= adaptive_WHEN_TO_STOP: break
```

`score(y, n, offset, a)` mirrors `discover.c::compute_score_right` with
the three-case structure (gap-in-sequence, diagonal, gap-in-consensus).
Genome-position mapping retains the existing forward/reverse formula
from `align_rebuild_consensus`, augmented with `+ offset` (forward) /
`- offset` (reverse) for the band.

Segdup-prevention checks (`< 2`, `> 1000 && < 3`, `> 5000 && < 5`) are
preserved — they now gate on `nrepeatocc` (count of instances with
positive score), the natural banded-DP analogue of the column-vote
`total < N` count.

The adaptive `WHEN_TO_STOP` (longer extensions earn longer quiet
windows, `max(WHEN_TO_STOP, besty/10)`) mirrors discover.c.

## 3. Build

```
$ make 2>&1 | tail -3
gcc -O3 -Wall -Wextra -std=c11 -pthread -march=native -flto -MMD -MP -c src/align.c -o obj/align.o
gcc -O3 -Wall -Wextra -std=c11 -pthread -march=native -flto -o bin/mdl-repeat ...
lto-wrapper: warning: using serial compilation of 2 LTRANS jobs
```

Zero compiler warnings with `-Wall -Wextra -O3`.

## 4. Test status

### 4a. `tests/run_tests.sh` — 7/7 PASS

```
=== MDL-Repeat Test Suite ===
--- Building ---
  PASS: Build succeeds with zero warnings
--- mdl.c unit tests ---
  34 assertions passed, 0 failed
  PASS: mdl.c unit tests
--- Test A: Multiple families (1MB genome) ---
  Accepted families:    3 / 3
  Bases covered:        32150 / 1011000 (3.2%)
  PASS
--- Test B: Tandem array (MVP limitation) ---
  INFO: Test B found 0 families
--- Test C: Detection limit ---
  Accepted families:    1 / 1
  Bases covered:        1200 / 211000 (0.6%)
  PASS
--- Test D: Nested TE (SINE inside LINE) ---
  Accepted families:    2 / 2
  Bases covered:        12800 / 511000 (2.5%)
  PASS
--- Human3M smoke test ---
  Accepted families:    174 / 325
  Bases covered:        1802157 / 3011000 (59.9%)
  Compression ratio:    0.9803
  PASS
===========================
Results: 7 PASS, 0 FAIL
===========================
```

testD finds both LINE and SINE — exactly the case the spec called out
as a sanity check.  Banded DP did not regress on the synthetic suite.

### 4b. `tests/test_mdl.c` — 34/34 PASS

```
$ gcc -O2 -Wall -Wextra -std=c11 -I src -o /tmp/test_mdl tests/test_mdl.c src/mdl.c -lm
$ /tmp/test_mdl 2>&1 | tail -3
  PASS: bases_covered == 1200 (got 1200)
  PASS: compression_ratio in [0,1] (got 0.9932)
34 passed, 0 failed
```

### 4c. `tools/test_bed_pr.py` — 17/17 PASS

```
$ python3 tools/test_bed_pr.py 2>&1 | tail -3
  PASS: main() returns 0 on valid inputs
17 passed, 0 failed
```

### Standalone testD verification

```
$ bin/mdl-repeat -sequence tests/data/testD.fa ... 2>/dev/null
$ awk '/^>/{...}' /tmp/.../test.fa
>R=1 1000      # LINE  (truth: 1000bp)
>R=0  200      # SINE  (truth: 200bp)
```

Recovers the exact truth lengths.

## 5. chr4 wall-clock benchmark

Single-threaded discover dominates wall-clock; refine is the part the
banded DP touches.

| Run         | Wall clock | User time | Peak RSS |
|-------------|------------|-----------|----------|
| chr4_v2 (column-vote, baseline per HANDOFF) | ~8 min | n/a       | n/a       |
| chr4_K (banded DP)                          | 7:49.59 | 494.19 s  | 1523 MB |

Banded DP runtime is essentially the same as column-vote.  Per-column
DP is more expensive, but the workload is data-bound: extension stops
at the same boundary, and discover (single-thread l-mer counting per
HANDOFF.md sec 6) is still the dominant cost.

## 6. chr4 family count + length distribution

| Metric                  | chr4_v2 (baseline) | chr4_K (banded DP) |
|-------------------------|--------------------|--------------------|
| Families output         | 2932               | 2920               |
| Total consensus bp      | 746,631            | 750,728            |
| Mean consensus length   | 254.6              | 257.1              |
| Max consensus length    | **20014**          | **20014**          |
| ≥ 500 bp count          | 326                | 326                |
| ≥ 2000 bp count         | 25                 | 27                 |
| ≥ 5000 bp count         | 2                  | 2                  |
| ≥ 10 000 bp count       | 1                  | 1                  |

Top 10 longest in chr4_K: 3306, 3317, 3390, 3928, 4209, 4288, 4390,
4439, 7633, 20014.

Top 10 longest in chr4_v2: 3214, 3306, 3317, 3390, 4209, 4288, 4439,
4439, 7912, 20014.

## 7. Unexpected complications / observations

**A. The longest consensus did not exceed 20014 bp.**

The spec predicted at least one family to grow past 20014 bp.  In this
chr4 run, the maximum stayed at exactly 20014 bp — identical to the
column-vote baseline.

Investigation: 20014 = `2*L + l = 2*10000 + 14` is the discover.c
extension cap.  Multiple discover families saturated this cap (`Warning:
extended right all the way to 20014`).  In the merge step, the maximum
briefly grew to 20255 (`Consensus length: max=20255`), confirming refine
did extend SOMETHING past 20014 — but that family was subsequently
pruned by MDL or coalesced, leaving the final max at 20014.

The fact that the banded DP did not grow this particular family beyond
20014 is consistent with two possibilities:
1. The instances genuinely end at the chromosome boundary or non-repeat
   flank — both column-vote and banded DP would stop at the same place.
2. The DP score plateau / `WHEN_TO_STOP` triggers immediately past
   20014 because the flanking sequence is non-repetitive across copies.

Total bp went up by 4097 bp (750728 vs 746631) and the ≥2000 bp count
went up by 2 (27 vs 25), so the banded DP is producing slightly more
total content — just not via a single longer family.

**B. Score initialization choice.**

Unlike discover.c which initialises `bestbest[n] = l*MATCH` (l-mer seed
score), the refine version starts from `bestbest[n] = 0`.  Rationale:
refine extension begins at the consensus edge, not after a fresh seed,
so there's no "accumulated seed score" to carry in.  Keeping the
relative scale matters — if we started bestbest at a positive value,
the `CAPPENALTY` floor would also be raised, suppressing instances
whose flanking diverges, and refine would terminate too quickly on
high-divergence families.

**C. MAXOFFSET in refine vs discover.**

`g_align_maxoffset` (default 12) is larger than discover's MAXOFFSET=5.
Band width is 2*12+1 = 25 vs discover's 11.  The DP cost scales as
band_width² inside Case C of `ext_compute_score`, so refine's per-
column cost is ~5× discover's — but column count per refine call is
capped at ALIGN_MAX_EXTENSION=10000, so total refine extension cost is
still bounded.  No measurable wall-clock regression.

## 8. Conclusion

Banded DP extension is in place, builds clean, passes all 7 synthetic
tests, all 34 mdl unit tests, all 17 bed_pr tests, and the chr4 smoke
run completes in 7:49 (no regression vs ~8 min column-vote baseline).
The output distribution is essentially equivalent to baseline with
slightly more total bp (4 kb spread across long families) and the
nested-TE testD case (LINE + SINE) is recovered exactly.

The single open caveat is that the longest consensus stayed pinned at
20014 bp (the discover.c cap).  Whether banded DP can push that family
further would require dropping ALIGN_MAX_EXTENSION or examining the
instance-edge mask for that specific family — both are out of scope
for option K (which only swaps the algorithm, not the limits).
