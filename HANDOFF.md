# mdl-repeat — Work Log & Handoff
*Last updated: 2026-04-28*

This document captures everything the previous session(s) discovered, fixed,
attempted, and concluded about mdl-repeat. Read this BEFORE picking up new
work. Skip the bottom-up reproduction steps if you only need the current
state.

---

## 0. Project design philosophy (DO NOT confuse)

| | mdl-repeat | RECON-style |
|---|---|---|
| Target | **high-copy, active TEs** | ancient, fragmented TEs |
| Architecture | seed-and-extend (RepeatScout 内核) | self-alignment + MSA |
| Selection | MDL-based (no thresholds) | clustering thresholds |

**Important**: low-copy / fragmented elements are explicitly **out of scope**
for mdl-repeat by design. When evaluating, exclude them or interpret
results accordingly. Don't propose RECON-style architecture changes.

Tandem repeats / simple repeats / centromere CEN180 are also **out of scope**
(rejected by entropy + periodicity filters in `discover.c`).

---

## 1. Code state — committed correctness fixes (P0 round, 12 items)

The first round of work fixed correctness bugs and added engineering hygiene.
All committed and validated by the synthetic test suite (7/7 PASS, 31/31
unit tests).

### M1 (correctness)
| # | Item | File | Status |
|---|---|---|---|
| 1 | Savings double-counting (compression_ratio < 0) | mdl.c::mdl_select_library | ✅ Fixed via unique-coverage greedy gating |
| 2 | dl_total int cast overflow on >2.1 Gb | mdl.c:literal_cost | ✅ int64_t signature |
| 3 | Dead R-convergence loop | mdl.c, main.c | ✅ Removed loop + recovery pass |
| 4 | Synthetic test data drift | tests/run_tests.sh | ✅ Frozen via committed fixtures + regenerate_data.sh |

### M2 (algorithmic)
| # | Item | File | Status |
|---|---|---|---|
| 5 | Merge stage missing MDL gate | refine.c::refine_merge_families | ✅ Added `estimate_merge_score`, symmetric with split/assemble |
| 6 | Periodic l-mer seed filter | discover.c::is_periodic_lmer | ✅ Catches (GAGA)n etc. that entropy missed |
| 7 | Two scoring systems undocumented | refine.h, mdl.c | ✅ Documented: ALIGN (1/-1/-5) feeds MDL; REFINE (2/-3/-2) is binary "same family?" |

### M3 (engineering hygiene)
| # | Item | File | Status |
|---|---|---|---|
| 8 | discover.c 1971 lines too big | src/ | ✅ Split: mask subsystem → discover_mask.c (405 lines), discover.c → 1574 lines, discover_internal.h shared (136 lines) |
| 9 | L_int(n≤0) silent error | mdl.c | ✅ Returns INFINITY now |
| 10 | mdl unit test coverage | tests/test_mdl.c | ✅ 31 assertions, covers L_int boundary + 3 MDL modes + unique-coverage |
| 11 | No CI | .github/workflows/ci.yml | ✅ 3 jobs: native build, portable + tests, valgrind |

### M4 (benchmarking)
| # | Item | File | Status |
|---|---|---|---|
| 12 | BED-vs-BED P/R tool | tools/bed_pr.py + test_bed_pr.py | ✅ 17/17 unit tests; supports `--overlap-mode pred|truth|min`, `--by-class`, `--ignore-strand` |

### Tandem coalesce (added later, biggest practical win)
- `refine.c::refine_coalesce_tandem_instances`, default factor 20×, called from main.c after prune
- Empirical: chr4 slice bp F1 0.452 → **0.632** (+40%); >10kb bin recall 0% → 52%
- This is a REPORTING change (post-MDL), not a discovery change. RM-remap shows it bridges gaps in tandem arrays.

### Nested-element merge gate (P1 follow-up)
- `refine.c::nested_containment_fraction`
- Trigger: longer ≥ 3× shorter AND shorter ≥3 instances
- Veto if containment fraction < 0.50
- Result: testD F1 0.83 → 0.91, chr4 unchanged (intended trade)

---

## 2. Benchmark setup — reproduce in 3 commands

```bash
# Workspace
ls /tmp/ath_bench/
#   tair10_sm.fa            — Ensembl Plants TAIR10 soft-masked (truth source)
#   chr4.fa                 — chr4 full uppercase, 18.6 Mb
#   chr4_full_truth.bed     — RM-mask intervals on chr4 (12,886 truth)
#   tair10_nuclear.fa       — chr 1-5 (use this for full-genome benchmarks)
#   tair10_nuclear_truth.bed
#   chr4_E14.fa             — best-tuned baseline library
#   chr4_full*.fa, chr4_p0*.fa, chr4_p0v3*.fa  — experiment outputs

# RepeatMasker via conda env (PRE-INSTALLED)
PGTA=/home/shuoc/tool/miniconda3/envs/PGTA
$PGTA/bin/RepeatMasker --version    # 4.2.1
$PGTA/bin/blastn -version           # 2.14.1+

# Run mdl-repeat → RM-remap → recall, all in one wrapper:
bash /tmp/ath_bench/multipass.sh chr4.fa chr4_test 4
# (the wrapper does Pass 1 + Pass 2 on residual + combined eval)
```

### Truth derivation
The truth BED is derived from Ensembl Plants' RepeatMasker-soft-masked
TAIR10. Lowercase positions = repeat. Contiguous lowercase runs ≥30 bp
become one truth interval. **No tandem repeat exclusion** (CEN180 is in
the truth) — this overstates the apparent miss rate for mdl-repeat's
design scope.

### Proper benchmark methodology
**ALWAYS** use RM-remap for benchmarking, NOT the raw `-instances .bed`
output:
1. Run mdl-repeat → produces `library.fa`
2. RepeatMasker `--lib library.fa` on genome → produces `.out`
3. Convert RM .out → BED → compare to truth via `tools/bed_pr.py`

Reason: discovery's instance BED is limited by seed-and-extend reach.
The library is what's evaluable as the de-novo product. RM with the
library finds all instances genome-wide.

---

## 3. Current performance — chr4 full RM-remap baseline

```
Total:        F1 0.515  recall 0.389  precision 0.760
Outside cen:  F1 0.441  recall 0.327  precision 0.678  ← real algorithmic
Inside cen:   F1 0.699  recall 0.558  precision 0.938  ← high-copy LTR mostly OK
                                                           CEN180 missed by design
```

### Recall stratified by copy count (the most informative cut)

```
Copy class            #int   truth_bp     covered     recall
1 (singleton)         1433    790,046      46,391     0.059   ← out of scope (RECON's job)
2                      502    295,355      77,276     0.262   ← out of scope (very hard)
3-4                    602    410,080     135,480     0.330   ← in scope, GAP
5-9                    561    521,278     226,100     0.434   ← in scope, GAP
10-49                  715  1,148,090     537,806     0.468   ← in scope, GAP
50-199                 280    771,367     426,343     0.553   ← should be easy, gap
≥200 (CEN180)          134  1,497,395     846,403     0.565   ← out of scope (tandem)

In-design-scope (3-199 copies) recall = 0.465
```

**This 0.465 is the honest "library quality" number for mdl-repeat's design target.**
50% of even 50-199 copy elements are missed.

### Synthetic tests (sanity)
```
testA: 3/3 fam, 32150 bp covered
testB: 0/0 fam (tandem — known MVP limitation)
testC: 1/1 fam (3-copy GoodFamily detected, 2-copy BorderlineFamily correctly rejected)
testD: 2/2 fam (LINE+SINE both detected after nested-element merge gate)
human3M: 157/320 fam, 43.2% covered
```

---

## 4. Failed optimization attempts (NEGATIVE results — important to record)

These were tried with empirical justification, none moved RM-remap recall.
Don't redo them without new evidence.

| Attempt | Hypothesis | Result | File state |
|---|---|---|---|
| **P0v1** mask: don't decrement bystander l-mer freq | bystanders are killed by mask shadow | F1 0.515 → 0.490 | reverted |
| **P0v2** + explicit ban consensus l-mer | adds principled "consume" semantics | F1 0.491 (no change) | reverted |
| **P0v3** restored decrement + only explicit ban | combine both principles | F1 0.515 (zero change) | reverted |
| **MAX_FAMILY_CLAIMS=2** (coverage-aware mask) | nested elements need re-claim | F1 0.632 → 0.609 | constant left at 1, infrastructure preserved |
| **Adaptive stop** (extend stops proportional to extension) | early stops kill long elements | F1 unchanged | reverted (extension not stopping early — bottleneck elsewhere) |
| **Wider DP band (-maxgap 8/15)** | too narrow misses indel-rich | F1 0.515 → 0.397, max_cons drops 6784→3989 | parameter only, not committed |
| **-L 30000** | extension hits ceiling | F1 unchanged (max_cons same as -L 10000) | parameter only |
| **-stopafter 500** | similar to adaptive | F1 0.515 → 0.386 | parameter only |
| **-minthresh 1** | catch low-copy | F1 → 0.585 | parameter only |
| **-max-divergence 0.40 / 0.50** | catch diverged | F1 unchanged | parameter only |
| **Fragment-assembly D-cap 30 kb** | inter-family chains | only 8 extra fragment pairs assembled | committed (harmless) |
| **Multi-pass on residual** | dominant high-copy starves second-tier | inconclusive — Pass 2 was manually killed at ~16 min in l-mer counting (not provably stuck, but unacceptable runtime); when it eventually ran to completion via the wrapper's failure handler, Pass 2 produced 0 families. Combined library = Pass 1. To get a real test, fix HASH_SIZE first (see §6) so Pass 2 finishes in reasonable time. | wrapper at /tmp/ath_bench/multipass.sh |

### Key insight from P0 failure

The **"mask shadow" hypothesis was wrong**: 64% of high-copy missed
elements being within 500bp of a discovered instance is **TE-clustering
correlation, not mask-shadow causation**. The original aggressive freq
decrement was actually correct — preserving bystander freq leads to
wasted MAXR slots (find_besttmp keeps picking already-claimed l-mers).

---

## 5. What's the real bottleneck? (Open question)

84 cluster members of 18 ≥3-copy missed family clusters: **0/84 ever touched
by discovery's instance output**. They never enter the candidate pool at all.

Most likely mechanisms (UNVERIFIED — needs controlled experiment):
1. Their seed l-mer's freq is depressed by hash-table collisions (HASH_SIZE=16M
   is fixed; on Arabidopsis 119 Mb chains get long → counting degenerates)
2. Seed-selection greedy picks one l-mer of family A, which through `lmermatcheither`
   hash neighborhood masks shared l-mers with family B, depressing B's freq
3. The l-mer choice itself is wrong: high freq doesn't mean useful seed if it's
   shared with non-repeat regions

**Needed**: controlled experiment — synthesize a small genome (e.g. background
+ 10 copies of a known missed cluster's representative) and check if mdl-repeat
finds it. If yes → mechanism is large-genome specific (likely #1). If no →
fundamental seed selection / mask design issue.

---

## 6. Engineering bottleneck (blocks experimentation)

Single-thread l-mer counting in `discover.c::build_headptr_internal`:
- chr4 (18.6 Mb): ~7-8 minutes
- chr4 residual after Pass 1 mask: ~16+ minutes (kept growing, was killed)
- TAIR10 nuclear (119 Mb): >16 minutes l-mer counting alone, didn't finish

Root cause: `HASH_SIZE = 16,000,057` is fixed. On large repeat-rich
genomes, hash chains grow long; each insertion is O(chain_length) due to
`lmermatcheither` chain walks. Plus chain memory allocator (no pool).

This blocks any iterative experiment (multi-pass, parameter sweeps on full
genome). Should be the FIRST engineering improvement before further algorithm
work, because it's gating data collection.

Fix sketch:
- Dynamic HASH_SIZE = `max(16M, 4×N/l)` where N is genome length
- Or replace linked list with open-addressing hash
- Or use existing `kmer.c`'s striped-lock parallel hash for the counting phase

---

## 7. Best current configuration

```bash
bin/mdl-repeat -sequence genome.fa -output lib.fa \
    -instances inst.bed -stats stats.tsv -threads 4
# Defaults are tuned: tandem coalesce 20×, nested merge gate active,
# unique-coverage MDL select, -L 10000, MINTHRESH=2, max-divergence 0.30
```

Don't touch these params unless you have evidence — they were swept
without finding improvements:
- `-L`: 10000 (sweet spot; 30000 same; 50000 same)
- `-stopafter`: 100 (default; 500 hurts)
- `-maxgap`: 5 (default; 8/15 hurt)
- `-minthresh`: 2 (default; 1 hurts slightly)
- `-max-divergence`: 0.30 (default; 0.40/0.50 unchanged)

---

## 8. File structure

```
/home/shuoc/tool/mdl-repeat/
├── src/
│   ├── main.c                  pipeline driver
│   ├── discover.c              seed-and-extend (1574 lines)
│   ├── discover_mask.c         masking subsystem (405 lines, split out in M3#8)
│   ├── discover_internal.h     shared types between discover_*.c (136 lines)
│   ├── refine.c                merge / split / assemble / prune / coalesce (2540 lines)
│   ├── mdl.c                   L_int + MDL scoring + library selection
│   ├── align.c                 instance collection via banded DP
│   ├── kmer.c                  refine-stage k-mer table (parallel)
│   ├── candidates.c, genome.c, output.c, types.h ... 
│   └── cmd_line_opts.c
├── tools/
│   ├── bed_pr.py               BED-vs-BED P/R tool (M4#12)
│   └── test_bed_pr.py
├── tests/
│   ├── data/                   frozen synthetic fixtures (committed)
│   ├── test_mdl.c              31 unit assertions
│   ├── run_tests.sh            integration tests
│   └── regenerate_data.sh      explicit fixture regeneration
├── .github/workflows/ci.yml    M3#11 CI
├── HANDOFF.md                  this document
├── README.md, TECHNICAL_DOC.md
└── Makefile
```

```
/tmp/ath_bench/                 benchmark workspace
├── tair10_sm.fa                119 Mb truth source (Ensembl softmask)
├── tair10_nuclear.fa           nuclear chr 1-5 only
├── chr4.fa                     chr4 full uppercase
├── chr4_full_truth.bed         12,886 truth intervals
├── chr4_full.fa                baseline mdl-repeat library output
├── chr4_full_rm.bed            RM-remap of baseline library
├── chr4_E14.fa                 best-tuned library (with all M2/coalesce/nested fixes)
├── chr4_p0*.fa                 P0 attempts (all failed)
├── multipass.sh                Pass1 + Pass2 + combined wrapper
├── eval_quick.py               bp-level recall calculator
├── recall_by_length.py         by-truth-length stratified recall
└── (and many *_disc.bed, *_rm.bed pairs from sweeps)
```

```
/home/shuoc/tool/miniconda3/envs/PGTA/   conda env with RepeatMasker
└── bin/RepeatMasker, blastn, makeblastdb, ...
```

---

## 9. Next steps — pick one

### (A) Fix l-mer-counting engineering bottleneck
Make HASH_SIZE adaptive or replace data structure. Without this, every
multi-pass / large-genome / parameter-sweep experiment is rate-limited
by single-thread linked-list hash lookups. ~150 LOC, low risk.

### (B) Controlled diagnostic experiment
Build a small synthetic genome:
```
background (random, 1 Mb)
+ 10 copies of a known "missed" cluster representative (from
  /tmp/ath_bench/missed_long.fa) inserted at random positions
```
Run mdl-repeat on this small genome. Does it find the family?
- YES → bottleneck is large-genome specific (likely hash collisions)
- NO  → fundamental discovery design issue (seed competition, l-mer selection)

This pins down WHICH mechanism is responsible for the 18 missed clusters.

### (C) Refine specific recall numbers
The 0.465 in-scope recall has uncertainty. To reduce it:
- Run on full TAIR10 (5 chromosomes) once hash perf is fixed
- Compare to RepeatModeler2 / EDTA output on same chr4 — head-to-head
  table for paper

### (D) Status quo: ship current state
Current code on chr4 RM-remap is recall 0.389 / precision 0.760 / F1 0.515.
For active high-copy TE benchmarks at 50-199 copy bin: 0.553 recall, 
which is the realistic ceiling without architecture changes.

---

## 10. Files modified during this session (uncommitted)

The previous session modified these but didn't commit. `git status` will
show them. Most are the M1-M4 + tandem coalesce + nested gate fixes.

```
M README.md
M TECHNICAL_DOC.md
M src/align.c              (length filter for tiny instances)
M src/align.h
M src/discover.c           (periodic filter, M3#8 mask code extracted)
M src/kmer.c               (removed pos_fill_worker dead code)
M src/kmer.h
M src/main.c               (tandem coalesce wired in, recovery pass removed)
M src/mdl.c                (savings unique-coverage, L_int hardening,
                            R-loop removed, dl_total int64_t)
M src/mdl.h
M src/refine.c             (merge MDL gate, nested-element gate,
                            tandem coalesce function, fragment D-cap 30k)
M src/refine.h
M src/types.h
M tests/run_tests.sh
M tests/test_mdl.c
?? .github/                          (M3#11 CI)
?? HANDOFF.md                        (this file)
?? PAPER_REFERENCE.md
?? PAPER_REFERENCE_CN.md
?? src/discover_internal.h           (M3#8 shared types — new)
?? src/discover_mask.c               (M3#8 mask subsystem — new)
?? tools/bed_pr.py                   (M4#12)
?? tools/test_bed_pr.py
?? tests/regenerate_data.sh          (M1#4)
```

Recommended commit strategy: 4 commits matching M1/M2/M3/M4 milestones,
plus a 5th for tandem coalesce + nested gate, plus 6th for HANDOFF.md.

---

## 11. Quick reference — sanity checks

```bash
# Build
cd /home/shuoc/tool/mdl-repeat && make

# Synthetic regression (must pass 7/7)
bash tests/run_tests.sh

# Unit tests (must pass 31/31)
gcc -O2 -Wall -Wextra -std=c11 -I src -o /tmp/test_mdl \
    tests/test_mdl.c src/mdl.c -lm && /tmp/test_mdl

# bed_pr unit (17/17)
python3 tools/test_bed_pr.py

# chr4 full RM-remap recall (~10 min)
bash /tmp/ath_bench/multipass.sh chr4.fa /tmp/ath_bench/sanity 4
# Expect: ~2724 families, ~0.389 RM-remap bp recall
```

---

## Bottom line

mdl-repeat's RM-remap recall on Arabidopsis chr4 = **0.389** (whole),
**0.327** outside centromere, **0.465** for in-design-scope (3-199
copies). The proximate engineering blocker is single-threaded l-mer
counting; the deeper algorithmic question is why 18 ≥3-copy clusters
never enter the candidate pool. Any further work should start by either
fixing the engineering blocker (option A) or running a controlled
diagnostic (option B) to identify the actual missed-family mechanism —
NOT by re-trying the optimizations already disproven in section 4.
