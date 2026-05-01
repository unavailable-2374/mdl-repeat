# V6 Patch Simple Result

Date: 2026-04-28

---

## Patch 1 — ENG-N3: Dynamic HASH_SIZE

### Files Changed

- `src/discover_internal.h`
  - Line 22: `#define HASH_SIZE 16000057` → `#define HASH_SIZE_MIN 16000057`
  - DiscoverContext struct: added `gpos_t hash_size` field (before the `int l` field)

- `src/discover.c`
  - `hash_function()` (lines 49, 52): `% HASH_SIZE` → `% (int)C->hash_size`
  - `build_headptr_from_freq()` (line 211): init loop uses `C->hash_size`
  - `build_headptr_internal()` (line 275): init loop uses `C->hash_size`
  - `build_headptr_parallel()` (line 345): init loop uses `C->hash_size`
  - `trim_headptr()` (line 468): loop uses `C->hash_size`
  - `build_all_pos()` (line 539): TANDEMDIST pass loop uses `C->hash_size`
  - `find_besttmp()` (lines 601, 615): both loops use `C->hash_size`
  - `free_headptr()`: signature changed to `free_headptr(struct llist **headptr, gpos_t hash_size)`; all 3 call sites updated to pass `C->hash_size`
  - `dump_freq_table()` (lines 1412, 1426): loops use `C->hash_size` (loop vars widened to `gpos_t`)
  - `discover_families()`: dynamic size computed as `max(HASH_SIZE_MIN, 4 * orig_length / l)`, rounded to odd; `headptr = calloc(C->hash_size, ...)` (was `malloc(HASH_SIZE, ...)`)

Confirmed: `discover_mask.c` uses only `SMALLHASH_SIZE` — not touched.

### Test Status: PASS
- `make`: clean build, zero new warnings
- `bash tests/run_tests.sh`: 7/7 PASS
- `test_mdl`: 34/34 PASS
- `test_bed_pr.py`: 17/17 PASS

---

## Patch 2 — ENG-N7: -coalesce-factor CLI

### Files Changed

- `src/main.c`
  - Added `float coalesce_factor = 20.0f;` variable declaration (~line 811)
  - Added CLI parsing block for `-coalesce-factor` (after `-max-dp-cells` block)
  - Replaced hardcoded `refine_coalesce_tandem_instances(candidates, 20.0f, verbose)` with conditional call using `coalesce_factor`; factor=0 skips coalescing entirely
  - Usage string: added `-coalesce-factor` entry under `Optional (refinement)`

- `README.md`
  - Added `-coalesce-factor` entry to the Refinement options table

- `TECHNICAL_DOC.md`
  - Added Step 8c (tandem-instance coalescing) to the pipeline step list

### Test Status: PASS
- `make`: clean build, zero new warnings
- `bash tests/run_tests.sh`: 7/7 PASS
- `test_mdl`: 34/34 PASS
- `test_bed_pr.py`: 17/17 PASS

---

## Patch 3 — J': Bump ALIGN_MAX_EXTENSION

### Files Changed

- `src/align.c`
  - Line 21: `#define ALIGN_MAX_EXTENSION 10000` → `#define ALIGN_MAX_EXTENSION 20000`
  - Comment updated: "max bases to extend per direction per iteration"

Total per-iteration cap: 2 × 20000 = 40 000 bp. Across 10 iterations: 400 000 bp max.

### Test Status: PASS
- `make`: clean build, zero new warnings
- `bash tests/run_tests.sh`: 7/7 PASS
- `test_mdl`: 34/34 PASS
- `test_bed_pr.py`: 17/17 PASS

---

## chr4 Smoke Test (All 3 Patches Combined)

Command:
```
bin/mdl-repeat -sequence /tmp/ath_bench/chr4.fa -output /tmp/v6_patch_chr4_smoke.fa \
  -instances /tmp/v6_patch_chr4_smoke.bed -threads 4
```

Result: **PASS — no crash, completed normally**

- Genome: 18,596,056 bp (chr4)
- hash_size: 16,000,057 (genome too small to trigger dynamic expansion; min applies)
- Families accepted: 2915 / 3940 candidates
- Bases covered: 4,532,600 / 18,596,056 (24.4%)
- Coalesced: 1558 tandem instance pairs (coalesce_factor=20.0, default)
- Wall time: 8m53s (threads=4)
- Exit code: 0
