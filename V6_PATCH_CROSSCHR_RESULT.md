# V6 Tier 1.5a Cross-Chromosome Correctness Patch — Result

Applied: 2026-04-28

All four bundled fixes landed in `src/refine.c` (only). No other source files
modified. Tests/data/multichr/ added; tests/run_tests.sh extended with Test E.

## Region 1 — B+Q6: `nested_containment_fraction` (refine.c:219-268)

* Added cross-chromosome guard as the first body statement of the inner
  loop (`if (shorter->instances[i].seq_index != longer->instances[j].seq_index) continue;`).
* Tightened the gate predicate via the **"solo evidence"** route (chosen
  per spec — simpler than re-aligning the non-overlap region). The function
  now also tracks `solo` = number of shorter instances that have NO overlap
  with any longer instance. If `solo < 2`, the function returns 1.0 so the
  caller's `ctf < 0.50f` veto does NOT fire — i.e. structural prefix/suffix
  relationships (longer = X + shorter, shorter never living on its own)
  no longer block merges.
* Lines changed: ~25 (function body partially rewritten with solo
  accounting + comment block).

## Region 2 — ENG-N11: `check_instance_overlap` (refine.c:283-345)

* Added cross-chromosome guard as the first statement inside both inner
  loops:
  - Sorted-path scan-forward loop (line 286).
  - Linear-fallback inner loop (line 314).
* Confirmed `Instance.seq_index` already exists at candidates.h:22 (no
  header change required).
* Lines changed: 2 inserts + 4 comment lines.

## Region 3 — ENG-N8: `refine_coalesce_tandem_instances` (refine.c:2538-2558)

* Inserted the cross-chromosome guard as the first body statement of the
  gap-merging walk loop. When `cur->seq_index != active->seq_index`,
  advance the active pointer (and `active_idx` for consistency) and
  `continue`. The gap test never fires across chromosomes.
* Per R3 directive: did NOT change the sort key. Runtime guard alone is
  sufficient.
* Lines changed: 9 inserts (guard + comment).

## Region 4 — ENG-N9: `refine_assemble_fragments` (refine.c:1864-2007)

* Added `int seq_index;` field to `InstanceEntry` struct (line 1900).
* Populated the field at line 2000 (`entries[ei].seq_index = inst->seq_index;`).
* Modified `cmp_instance_entry` to sort by `(seq_index, start)` —
  primary key seq_index ascending, secondary key start ascending
  (line 1916-1942).
* Per R3 directive: did NOT add an in-body `seq_index` check; the existing
  `break` at the head of the inner sweep loop continues to fire on the
  start-distance bound.
* Lines changed: 1 struct field + 1 populate + ~12 comparator (with comment).

## Test Status

```
=== MDL-Repeat Test Suite ===
  PASS: Build succeeds with zero warnings
  PASS: mdl.c unit tests (34/34)
  PASS: Test A: found families (>= 1)
  PASS: Test A: found instances (>= 50)
  PASS: Test C: found at least 1 family
  PASS: Test D: found families (>= 1)
  PASS: Test E (-threads 1): >= 3 families found (no cross-chr merge collapse)
  PASS: Test E (-threads 1): <= 8 families found (no spurious split runaway)
  PASS: Test E (-threads 1): no BED interval spans chr boundary (ENG-N8)
  PASS: Test E (-threads 1): every BED row has valid strand
  PASS: Test E (-threads 4): >= 3 families found (no cross-chr merge collapse)
  PASS: Test E (-threads 4): <= 8 families found (no spurious split runaway)
  PASS: Test E (-threads 4): no BED interval spans chr boundary (ENG-N8)
  PASS: Test E (-threads 4): every BED row has valid strand
  PASS: Test E: family count agrees within 2 across thread modes (4 vs 4)
  PASS: Human3M: found families (>= 10)

Results: 16 PASS, 0 FAIL
```

(Test B's tandem-array detection is a known MVP limitation noted in
CLAUDE.md / pre-patch baselines and is reported as INFO, not FAIL.)

```
$ python3 tools/test_bed_pr.py
17 passed, 0 failed
```

## Multi-chromosome Test Result

* Synthetic 2-chromosome genome: 2 × 100 kb chromosomes, each with 12
  repeat instances drawn from 2 distinct random consensus sequences (4
  truth families total: RepA/B on chr1, RepC/D on chr2).
* `-threads 1`: 4 families produced, 29 BED instances, no row crosses
  chr1/chr2 boundary, no spurious cross-chr family collapse.
* `-threads 4`: 4 families produced, 32 BED instances, identical
  invariants hold.
* Family count agrees exactly between thread modes.

## chr4 Smoke Test

* Input: /tmp/ath_bench/chr4.fa (18.6 Mb)
* Command: `bin/mdl-repeat -sequence chr4.fa -threads 4 -output ... -instances ... -stats ...`
* Wall time: **9m 8s** (within the 5-10 min spec).
* Result:
  - Discovered families: 4022 (post-compaction)
  - Post-merge: 3952 (66 merged)
  - Post-split: 3952 (0 split)
  - Post-fragment-assembly: 3933 (19 pairs assembled)
  - MDL-accepted: 2936 / 3933
  - Final library: 2914 families
* Baseline `chr4_full.fa`: 2725 families.
* Recent baselines `chr4_v2.fa` / `chr4_traced.fa`: 2932 / 2769 families.
* Patched output (2914) sits in the recent-baseline band; +6.9% vs. the
  oldest `chr4_full.fa` baseline. Sub-component counts (merge, split,
  assemble) are within unit-digit deltas of the historical numbers.
  No crash, no infinite loop, no warning newly emitted.

## Files Modified

* `src/refine.c` — all 4 fixes (no other source files touched).
* `tests/run_tests.sh` — added "Test E: Multi-chromosome cross-boundary
  correctness" block (9 assertions, both -threads 1 and -threads 4).

## Files Added

* `tests/data/multichr/gen_multichr.py` — deterministic generator
  (Python random with fixed seed) for the 2-chr synthetic FASTA.
* `tests/data/multichr/multichr.fa` — 2 × 100 kb chromosomes.
* `tests/data/multichr/multichr_truth.bed` — 24 truth intervals.

## Decisions and Trade-offs

* **B+Q6 gate variant**: chose "solo evidence" (≥ 2 shorter instances
  with NO overlap to any longer instance ⇒ block merge) over the
  alternative "non-overlap-region identity" check. Reason: simpler,
  no DP, no extra alignment, easy to reason about.
* **ENG-N8**: did not change the sort key (per R3); the tandem coalesce
  walk uses position-only sort and the runtime cross-chr guard simply
  resets the active pointer when crossing a sequence boundary.
* **ENG-N9**: per R3 directive, relied on the (seq_index, start)
  composite sort key + the existing `break` for boundary handling rather
  than adding a per-iteration seq_index check. Note that since positions
  are stored in a concatenated padded coordinate system without an
  inter-sequence padding gap, the `break` predicate alone may not catch
  every cross-boundary pair within distance D — but the new sort
  predicate keeps all chr-1 entries before any chr-2 entry, and the
  cross-chr leak window is narrow (<= D, ≈ 30 kb out of 100 kb in the
  test fixture, much smaller for real genomes) and was empirically
  benign in Test E and chr4. If a tighter guarantee is required later,
  a single-line `if (entries[j].seq_index != entries[i].seq_index) break;`
  inside the inner sweep would be the cleanest extension.
* **Tandem coalesce loop variable**: kept `active_idx` updated alongside
  `active` in the new cross-chr branch even though the spec only
  required `active = cur`. Reason: the existing else-branch updates
  both, so consistency is preserved and the variable stays meaningful
  throughout the loop.
