# CLAUDE.md — mdl-repeat design & logic reference

A model-facing reference for the **current source** in `src/`. It describes what
the tool is, how the algorithm works, and the load-bearing invariants — grounded
in the code, not in the older prose docs. Where this document and
`README.md` / `TECHNICAL_DOC.md` disagree, **this document and the code win**;
those two are partially stale (see [§13 Known doc drift](#13-known-doc-drift)).
Session history, benchmark numbers, and the log of reverted experiments live in
`FINAL_REPORT.md` and are deliberately not repeated here.

---

## 1. What it is

mdl-repeat is a **de novo repeat / transposable-element (TE) family library
builder** for genomic DNA. Given a FASTA genome and no prior repeat library, it
reports a set of consensus sequences (one per repeat family) plus the genomic
instances of each.

The core idea is **RepeatScout-style seed-and-extend discovery** combined with
**Minimum Description Length (MDL) model selection**: a family is kept only if
encoding its copies as edits against a shared consensus costs fewer bits than
encoding those copies literally. MDL is what decides *which* families and *how
many* end up in the library — there is no manual recall/precision threshold to
tune.

**In scope:** high-copy, low-divergence interspersed repeats — LTR
retrotransposons, LINEs, SINEs, MITEs, DNA transposons. Families need roughly
**≥3 copies** and a consensus of meaningful length to clear the MDL gate.

**Out of scope (by design, not bugs):** 1–2 copy fragments (RECON's job),
tandem/satellite arrays and centromeric repeats (TRF's job, and actively
filtered out via tandem-distance and periodicity filters), and low-complexity
sequence (entropy-filtered).

Implementation: **~12,000 lines of C11**, no build-time external libraries
(`-lm -pthread` only). Single optional *runtime* tool: `seqkit` (opt-in QC only;
see [§10](#10-external-tool-reality)).

---

## 2. Build / run / test

```bash
make                 # bin/mdl-repeat  (gcc -O3 -march=native -flto)
make PORTABLE=1      # without -march=native (portable binary)
make clean

# Minimal run
bin/mdl-repeat -sequence genome.fa -output families.fa -v

# Full outputs + threads
bin/mdl-repeat -sequence genome.fa -output families.fa \
    -instances inst.bed -stats stats.tsv -threads 4 -vv

bin/mdl-repeat            # no args → prints the authoritative flag list
```

Tests:

```bash
bash tests/run_tests.sh                       # end-to-end synthetic genomes
gcc -O2 -std=c11 -I src -o /tmp/t tests/test_mdl.c       src/mdl.c    -lm && /tmp/t
gcc -O2 -std=c11 -I src -o /tmp/s tests/test_sweepline.c src/refine.c src/mdl.c -lm && /tmp/s
```

CI (`.github/workflows/ci.yml`) runs native build, portable build + test suite,
and valgrind on a small case on every push.

---

## 3. Theory: the MDL objective (`mdl.c`)

MDL minimizes the total cost of describing the genome as a library plus the
genome given that library:

```
DL(Genome) = DL(Library) + DL(Genome | Library)
```

All costs are in **bits**. The implemented pieces:

- **Rissanen universal integer code** `L_int(n)` — `mdl.c:18`. Encodes a positive
  integer in `log2*(n) + log2(c0)` bits, where `log2*` is the iterated-log sum
  and `LOG2_C0 = 1.5179605508986484` (= log2(2.865…)). `n ≤ 0` returns `INFINITY`
  on purpose, so malformed inputs reject a family rather than score as free.

- **Library cost per consensus** `mdl_model_cost = L_int(len) + 2·len`
  (`mdl.c:130`): a length header plus 2 bits/base for the literal consensus.

- **Per-instance cost** `mdl_instance_cost_full` (`mdl.c:73`):
  `L_int(a) + L_int(m+1) + m·log2(3) + [position term]`, where `a` = aligned
  length, `m` = edit count. The position term depends on `-mdl-mode`:
  - `none` — no position term.
  - `exact` (default) — `log2 C(a, m)` via `lgamma`, the exact number of ways
    to place `m` edits among `a` positions.
  - `upper` — `m·log2(a)`, a looser upper bound.

- **Per-family score** `mdl_score_family` (`mdl.c:143`):
  `total_savings = Σ_i (2·a_i − instance_cost_i)`, and
  `mdl_score = total_savings − model_cost`. A family is worth keeping when its
  copies compress (positive savings) by more than the cost of storing its
  consensus. Edit counts `m_i` come from the **discovery/alignment** DP
  (`1/−1/−5` scoring), not the refinement DP; a negative `num_edits` falls back
  to `divergence · a_i`.

> **Load-bearing caveat — costs are R-independent.** The current per-instance
> formula intentionally drops the family-count overhead (type bit, `log2(R)`
> family id, strand bit). `consensus_length` and `num_families` are passed in but
> `(void)`-ignored (`mdl.c:73`). Consequence: a family's score does **not** change
> with the number of accepted families R. This is why the old "iterate R to
> convergence" loop and the "post-prune recovery pass" were removed — they would
> be inert. Do not re-add R-dependent control flow without first re-introducing an
> R-dependent cost term.

---

## 4. Pipeline (exact current order, `main.c`)

```
genome_load                                              (1)
  └─ [if raw_length > sample-size] genome_sample_windows (1b, very large genomes)
default_k  → l-mer length
discover_families  OR  discover_chunked                  (2)  branch on raw_length > chunk-size
  └─ [if -recall-rescue] recall_rescue_run               (2b, optional rescue pass)
kmer_count → kmer_trim → kmer_build_positions            (3)  k-mer table + position index
compact (drop families with <2 instances or short cons)  (4)
refine_merge_families        (80-80-80 + union-find)     (5)
refine_split_families        (Otsu bimodality)           (6)
refine_assemble_fragments    (spatial co-occurrence)     (6b)
  └─ [if -recruit-short] align_blast_recruit_short_fam…   (6c, optional rmblastn recruit)
mdl_select_library           (single-pass MDL selection) (7)
refine_prune_families        (exclusive-coverage rule)   (8)
  └─ [if coalesce-factor > 0] refine_coalesce_tandem…    (8b, reporting-only)
output_fasta / output_bed / output_stats                 (9)  remap coords if sampled
  └─ [if -external-qc <file>] external_qc_run_seqkit…    (9b, opt-in QC)
```

Notes on ordering that matter for reasoning:

- **MDL selection runs *after* merge/split/assembly** — so it scores the refined
  family set, not the raw discovery output.
- **Short-family recruitment (`-recruit-short`, opt-in) sits after assembly and
  before selection** — it only *adds* instances to existing short families, so the
  boosted counts help the MDL/standalone gate without ever feeding the assembly
  sweep. Off by default; needs `rmblastn` ([§10](#10-external-tool-reality)).
- **Prune runs *after* MDL selection** and uses its own exclusive-coverage rule
  (not the MDL score) to drop marginal families.
- **Tandem coalescing and external QC are post-decision** — they change reported
  instances / emit a QC file but do not feed back into selection.
- Empty-output early exit if 0 families survive compaction.
- `-trace-dir` dumps the family set after stages 01–08 for debugging.

---

## 5. Discovery engine (`discover.c`, `align.c`, `kmer.c`, `discover_mask.c`)

Seed-and-extend over the genome, which is loaded **front-padded by `PADLENGTH =
11000`** so extensions never run off the array.

1. **L-mer frequency table.** Count every l-mer with **symmetric/canonical
   hashing** (a forward l-mer and its reverse complement share a bucket;
   `hash_function`, `discover.c:40`). Two filters at insert time:
   - **Tandem-distance filter** — a hit increments frequency only if it is
     ≥ `TANDEMDIST` (500 bp) from the previous same-strand hit, so tandem arrays
     don't inflate counts.
   - **Entropy + periodicity filters** — Shannon-entropy reject (`MAXENTROPY
     −0.70`) and a period scan (`PERIODIC_MATCH_PCT 85`) drop low-complexity and
     short-period seeds.
   - **Auto l-mer length** `l = ceil(1 + log4(N))` (`discover.c:1496`). Hash table
     size is **dynamic**: `max(16000057, 4N/l)` rounded odd (16M is only a floor).
     A 31-base packing cap exists only on the parallel counting path (`kmer.c`).

2. **Seed selection** (`find_besttmp`, `discover.c:594`): greedily take the
   highest-frequency l-mer not yet masked, with a locality shortcut that reuses
   the previous best hash before falling back to a full scan.

3. **N-sequence simultaneous banded extension** (`extend_right`/`extend_left`):
   align *all* N occurrences of the seed at once and grow the consensus one
   column at a time in both directions. For each candidate base ∈ {A,C,G,T} the
   engine re-runs the banded DP with that base imposed and sums the best score
   across all occurrences; the base maximizing that sum becomes the consensus
   column. Band half-width `MAXOFFSET = 5`; scoring `match 1 / mismatch −1 /
   gap −5`, with a per-occurrence floor at `CAPPENALTY −20`. Stop when no column
   improves total score by ≥ `MINIMPROVEMENT` (3) within an **adaptive** quiet
   window `max(WHEN_TO_STOP=100, extended/10)`. Keep the consensus only if it is
   ≥ `GOODLENGTH` (30 bp).

4. **Masking** (`mask_headptr`, `discover_mask.c`): re-find the new family's
   occurrences via 1-vs-1 banded DP against the consensus and mark those genome
   positions claimed (`MAX_FAMILY_CLAIMS = 1`; tandem-pruned positions get the
   `CLAIM_PERMANENT = 255` sentinel) so the next seed iteration can't re-discover
   the same family. Then go back to step 2.

**Alignment-time instance recruiter** (`align.c`): separate from discovery, the
position-indexed k-mer table is used to recruit/realign instances —
`seed_genome_scan` collects consensus k-mer hits, `cluster_seed_hits` groups them
into per-locus anchors (cluster window `1.5·cons_len`), then `align_banded`
refines each. Strand is resolved from `cons_is_rc == genome_is_rc`.

**Parallelism (deterministic).** L-mer counting — both the recruiter k-mer table
(`kmer.c`) and discovery's own l-mer table (`build_headptr_parallel`,
`discover.c`) — uses **bucket-ownership**: every worker scans the whole genome in
coordinate order but only processes hash buckets it exclusively owns
(`h % num_threads == tid`). Because the hash is RC-symmetric, an l-mer and its
reverse complement share a bucket and thus one owner, so each k-mer's occurrences
— and the order-dependent TANDEMDIST frequency filter — are evolved by a single
thread in genome order, byte-identically to the serial path and independent of
thread count. No locks. (This replaced an earlier striped-lock counter whose
TANDEMDIST filter raced across genome chunks, and unified discovery's formerly
separate serial/parallel counters which used different strand conventions.) The
position index is counted in parallel but filled in genome order under the
`KMER_MAX_POSITIONS = 50000` truncation cap. The discovery main loop itself is
serial; large genomes parallelize across *chunks* instead (see [§8](#8-large-genome-scaling)).

---

## 6. Refinement (`refine.c`)

Four transformations between discovery and MDL selection, plus a post-selection
prune and a reporting-only coalesce. Each transformation that changes family
content is **MDL-gated** so it can only help compression.

- **Merge — 80-80-80** (`refine_merge_families`). k-mer Jaccard pre-screen
  (`REFINE_SCREEN_K 8`, `MIN_JACCARD 0.15`) then semi-global DP on both strands.
  Merge when identity ≥ **0.80**, coverage ≥ **0.80**, aligned ≥ **80 bp**
  (a relaxed 0.70/0.70 tier is allowed only with instance-overlap confirmation).
  Guards: **length-ratio ≥ 0.7**, a **nested-element veto** (when one consensus
  is ≥3× the other and the short one has ≥3 copies, require ≥50% containment),
  and an **MDL veto** (`estimate_merge_score ≤ 0 → skip`). Transitive merges
  resolved by union-find. The consensus-vs-consensus DP scoring (`2/−3/−2`)
  produces identity/coverage only and never feeds MDL. **Deterministic
  resolution:** the per-component representative is the max-instance family with
  ties broken by lowest family index, and every representative that absorbed
  others is re-refined — both root-independent so the merged consensus is
  byte-identical regardless of thread scheduling (see [§14](#14-invariants--gotchas-read-before-editing)).

- **Split — Otsu bimodality** (`refine_split_families`). Build a 100-bin
  divergence histogram per family, find the Otsu threshold, and split when
  bimodality (inter-/total-variance) ≥ **0.20** (a valley-depth check handles the
  borderline 0.20–0.40 band). Needs ≥3 instances, ≥3 per cluster, min divergence
  gap 0.03. Accepted only if the **MDL** score of the two sub-families clears a
  relaxed gate (≥ 0 when the original scored positive).

- **Fragment assembly** (`refine_assemble_fragments`). A **spatial
  co-occurrence sweep-line** (not sequence similarity) finds families whose
  instances repeatedly sit within proximity `D = median_cons·4` (clamped
  500–30000 bp) of each other in the genome. Guards: ≥3 co-occurrences,
  same-direction ≥0.80, size-ratio ≥0.10, a **nesting guard** (skip if ≥50%
  containment), and gap-geometry sanity (median gap, MAD). Build the joined
  consensus and accept only if its MDL score exceeds the **sum** of the parts.

- **Prune — exclusive coverage** (`refine_prune_families`, sweep-line). Sort
  accepted families weakest-first; for each instance compute exclusive bases
  (positions covered by no other family). An instance "counts" only if ≥ **25%**
  of its length is exclusive; a family is pruned (`CAND_ACCEPT_PRUNED`) if **no**
  instance clears 25%. The operative criterion is this 25%-exclusive rule, not a
  direct score-vs-cost comparison.

- **Tandem coalesce** (`refine_coalesce_tandem_instances`, post-MDL,
  reporting-only). Merge consecutive same-strand instances separated by
  `−10 ≤ gap ≤ coalesce_factor·consensus_length` (floor 50). Default
  `coalesce-factor = 20.0`. Recomputes divergence; does not change which families
  are in the library.

**Cross-record safety (correctness-critical).** Every position comparison across
instances must guard on `Instance.seq_index` so positions from different FASTA
records are never treated as adjacent. The merge containment/overlap helpers, the
assembly sweep, and coalesce all carry these guards. **Any new code comparing
`.position` across families MUST check `seq_index` first** — this was a systemic
multi-record bug class.

---

## 7. MDL selection (`mdl_select_library`, `mdl.c:281`)

**Single pass, no R-convergence loop** (`R_estimate = num_families`, used once).
Families are sorted by score descending and admitted greedily under a
**two-branch gate**; only **exclusive** (non-overlapping) savings accumulate into
the reported total, preserving the two-part bound.

A family is admitted if **either**:

1. **Exclusive branch** — it has uncovered ("exclusive") instance bases and
   `exclusive_savings − model_cost > 0`; **or**
2. **Standalone-fallback branch** — `standalone_score > 0` **and** consensus
   length ≥ 50 **and** ≥ 3 instances. This rescues real families that a pure
   unique-coverage greedy would destroy because their copies overlap families
   already admitted. Admitted-by-fallback families are flagged
   `CAND_ACCEPT_STANDALONE` / `CAND_QF_STANDALONE_FALLBACK`.

Pre-screen reject: `mdl_score ≤ 0` or `< 2` instances.

**Coverage tracking is an interval sweep-line** (`MdlInterval`, binary search +
two-pointer merge), **not** a per-base bitmap — memory is O(#instances), which is
what makes multi-gigabase genomes feasible.

**Compression ratio** (`mdl.c:713`): `dl_total / (2N)` where
`dl_total = 2N − total_savings + dl_library`, clamped to ≥ 0. The clamp is the
only hard invariant enforced here.

---

## 8. Large-genome scaling (`main.c`)

Two mechanisms activate automatically on `raw_length` (unpadded bases):

- **Genome sampling** (> `sample-size`, default 1 Gb). Reduce to a representative
  sample of 1 Mb tiles selected by Fisher–Yates over `window-size` tiles
  (`-seed` for reproducibility). Discovery runs on the sample; instance
  coordinates are **remapped to the original genome** before output.

- **Chunked discovery** (> `chunk-size`, default 200 Mb). Segment long sequences
  with boundary overlap, LPT-balance the segments into per-thread bins, discover
  in parallel, and remap each chunk's coordinates. Each chunk computes its own
  l-mer length for sensitivity. Recursive halving kicks in above `chunk-size·1.8`.

All positions/lengths are 64-bit (`gpos_t`/`glen_t = int64_t`) to avoid overflow.

---

## 9. Recall rescue (`rescue.c`, opt-in `-recall-rescue`)

An optional second discovery pass (pipeline step 2b) that trades runtime for
sensitivity. It re-runs the **in-process** discovery engine with a **shorter
l-mer** (`l − rescue-l-delta`, floored at 8) to catch families the primary pass
missed, then appends non-duplicate families back into the candidate list flagged
`CAND_QF_RESCUE_DISCOVERY`.

- **Targeted mode** (default): build rescue segments from **uncovered genome
  gaps** (existing instance intervals merged; gaps ≥ `rescue-min-gap` flanked by
  ±L), re-discover only there, remap back.
- **Full-genome mode** (`-rescue-full-genome`): re-discover over the whole genome
  (chunked if large).
- **Duplicate gate**: a rescue family is dropped as a duplicate when
  length-ratio ≥ 0.80, identity ≥ 0.80, and containment ≥ 0.80 versus an existing
  family (forward/RC, sliding-window identity). `-rescue-audit <file>` records
  per-target / per-candidate decisions. No external tools involved.

---

## 10. External-tool reality

**By default the binary forks no subprocess** — every external-tool path below is
opt-in and off by default. This matters because older prose implies BLAST is
always in the loop; it is not.

- **rmblastn / BLAST short-element recruitment is WIRED but opt-in
  (`-recruit-short`).** `align_blast_recruit_short_families` (`align.c`) has a
  live caller at `main.c:1341` (pipeline **Step 6c**, see [§4](#4-pipeline-exact-current-order-mainc)),
  gated on the `-recruit-short` flag (default off). When enabled it batch-BLASTs
  every short (< `RMBLAST_SHORT_THRESHOLD` = 500 bp) family consensus against the
  genome and adds non-overlapping instances (div ≤ `max-divergence`, capped per
  family) so weakly-supported short families (SINEs/MITEs) clear the
  MDL/standalone admission gates. It runs **after assembly / before selection**,
  so the boosted instance counts never feed the assembly sweep (the 12×-explosion
  that caused the 2026-05-01 revert lived in the old bundled `align_refine_all`
  pre-pass, not here). Determinism is preserved: parsed hits are sorted into a
  canonical `(family,position,strand,score)` order before the order-sensitive
  overlap de-dup, so the recruited set is byte-identical at `rmblastn
  -num_threads 1` vs 16. `find_rmblastn` honors `$RMBLASTN_BIN`, else resolves
  `which rmblastn` from PATH — no hardcoded path; a missing tool is a soft
  warning (the flag is skipped). RepeatMasker is never invoked anywhere.
- **Two external tools the pipeline can run, both opt-in:**
  - **`rmblastn` + `makeblastdb`** — only under `-recruit-short` (above). The
    recruit path invokes them via raw **`system()`/`popen()`** command strings
    (`align.c:1447/1487/1615`), **not** through `tool_runner.c`, so it carries no
    timeout and builds/removes a scratch BLAST db next to the query FASTA.
  - **`seqkit stats`** — only when `-external-qc <file>` is given (which
    auto-promotes `-external-tools` to `auto`). It runs **non-mutating QC** on the
    already-written FASTA and writes a TSV; it never feeds back into
    discovery/refinement. `seqkit` is located via PATH or `-seqkit <path>`.
    Hard-fails only under `-external-tools require`; otherwise a missing tool is a
    soft warning.
- **`tool_runner.c`** is a generic, dependency-free `fork`/`execv` launcher with
  timeout (default 300 s); currently only `external_qc.c` (seqkit) calls it — the
  `-recruit-short` path deliberately does not (it uses `system()`), an
  inconsistency worth unifying if the recruit path grows.

---

## 11. Core data structures (`candidates.h`, `types.h`, `genome.h`)

- **`gpos_t` / `glen_t` = `int64_t`** — every genome position and length. 64-bit
  is mandatory for genomes > 2 Gb.
- **DNA encoding**: A/C/G/T = 0/1/2/3, N = 99; complement = `3 − c` (N-safe).
  `PADLENGTH = 11000` front-pad on every loaded genome.
- **`Genome`**: numeric `sequence` (padded), `length` (incl. pad),
  `raw_length` (real bases — drives sampling/chunk thresholds), `boundaries`
  (per-record), `sequence_ids`.
- **`Instance`**: `position` (padded coords), `aligned_length`, `cons_start/end`,
  `num_edits`, `divergence` (0–1), `score`, `strand` (±1), **`seq_index`** (FASTA
  record — correctness-critical for any cross-instance comparison).
- **`CandidateFamily`**: `consensus` (numeric), `consensus_length`, `topology`
  (linear/complex/cyclic), `estimated_copies`, `instances[]`, `discovery_flags`,
  and a `CandidateMdlState mdl` (model cost, standalone/exclusive savings &
  scores, exclusive bases/instances, accept state, quality tier). Legacy
  `mdl_score`/`model_cost` aliases are kept in sync during migration.
- **`CandidateAcceptState`**: `UNSCORED=0, REJECTED, EXCLUSIVE, STANDALONE,
  PRUNED` (0 = not-yet-selected so a `memset` init is valid).
- **`CAND_QF_*`** quality-flag bitset — notably `CAND_QF_RESCUE_DISCOVERY` and
  `CAND_QF_STANDALONE_FALLBACK`.

---

## 12. Parameters & defaults (current code)

| Param | Default | CLI | Stage |
|---|---|---|---|
| l (seed length) | auto `ceil(1+log4(N))` | `-l` | discovery |
| L (extension half-cap) | 10000 | `-L` | discovery |
| MINTHRESH | 2 | `-minthresh` | discovery |
| TANDEMDIST | 500 | `-tandemdist` | discovery |
| MAXOFFSET (band) | 5 | `-maxgap` | discovery |
| match / mismatch / gap / cap | 1 / −1 / −5 / −20 | `-match` … `-cappenalty` | discovery |
| MINIMPROVEMENT | 3 | `-minimprovement` | discovery |
| WHEN_TO_STOP | 100 (adaptive ·/10) | `-stopafter` | discovery |
| MAXENTROPY | −0.70 | `-maxentropy` | discovery |
| GOODLENGTH | 30 | `-goodlength` | discovery |
| hash size | `max(16000057, 4N/l)` odd | (auto) | discovery |
| max-divergence | 0.30 | `-max-divergence` | refine |
| refine-gap | −5 | `-refine-gap` | refine |
| refine-maxoffset | 12 (max 32) | `-refine-maxoffset` | refine |
| merge identity/coverage/aligned | 0.80 / 0.80 / 80 | (const) | refine |
| split bimodality | 0.20 | (const) | refine |
| max-dp-cells | 10,000,000 | `-max-dp-cells` | refine |
| coalesce-factor | 20.0 (0 = off) | `-coalesce-factor` | refine |
| MDL mode | exact | `-mdl-mode` none\|exact\|upper | select |
| standalone gate | len ≥ 50, ≥ 3 inst | (const) | select |
| chunk-size | 200 Mb | `-chunk-size` | large genome |
| sample-size | 1000 Mb | `-sample-size` | large genome |
| window-size | 1000 kb | `-window-size` | large genome |
| seed | 42 | `-seed` | large genome |
| recall-rescue | off | `-recall-rescue` | rescue |
| rescue-l-delta / maxrepeats / min-gap | 1 / 2000 / 200 | `-rescue-*` | rescue |
| recruit-short (rmblastn, < 500 bp fams) | off | `-recruit-short` (`$RMBLASTN_BIN`) | recruit |
| external QC | off | `-external-qc`, `-seqkit`, `-external-tools` | QC |

`bin/mdl-repeat` with no args prints the authoritative, code-accurate list.

---

## 13. Outputs

- **FASTA** (`-output`): one consensus per accepted family,
  `>R=<id> length=<L> copies=<n> mdl=<score>`.
- **BED6** (`-instances`): one row per instance; column 5 = `1000·(1−divergence)`,
  strand in column 6. Coordinates are remapped to the original genome if sampling
  was used.
- **TSV** (`-stats`): per-family consensus length, copy count, mean divergence,
  MDL score, model cost, topology.

---

## 14. Invariants & gotchas (read before editing)

- **`seq_index` guards** on every cross-instance position comparison. Omitting
  one silently fuses unrelated loci across FASTA records.
- **64-bit positions everywhere** — never store a genome offset in `int`.
- **MDL costs are R-independent** in the current formula. Do not add control flow
  that assumes scores change with the family count without first restoring an
  R-dependent cost term ([§3](#3-theory-the-mdl-objective-mdlc)).
- **`PADLENGTH` (11000) must stay ≥ the maximum single-side extension distance**,
  or boundary math in extension/masking/rescue breaks.
- **MDL invariants**: compression ratio ∈ [0,1]; `dl_total` clamped ≥ 0;
  only exclusive savings accumulate into the reported total. Preserve the
  exclusive-OR-standalone admit semantics if you touch selection.
- **Refinement transforms are MDL-gated.** Merge/split/assembly each have an MDL
  acceptance check so they cannot reduce compression — keep that property.
- **Fragment-assembly co-occurrence hash must stay growable, and the sweep must
  stay budgeted** (`refine_assemble_fragments`, `refine.c`). The pair-counter is
  an open-addressing hash over distinct co-occurring family pairs; it is resized
  (`coocc_ht_grow`, load factor ≤ ~0.7) so probe chains stay O(1). A **fixed**
  table saturates whenever distinct pairs exceed its slots (trivial on
  repeat-dense or over-split inputs): every lookup degrades to a full-table scan
  *and* excess pairs are silently dropped — a multi-hour spin plus corrupt
  co-occurrence counts. The `n > 50000` early-skip guards the *family* count,
  **not** the pair count, so it does not protect against this. The sweep is
  O(entries × window-density) and `entries` can reach millions once
  `-max-instances` is high and split is iterated; `SWEEP_BUDGET` (comparison cap)
  converts a pathological scale into a logged, safe skip (assembly is an optional
  MDL-gated enhancement) instead of an unbounded run. Do not reintroduce a
  fixed-size pair hash or remove the work budget.
- **Headline metric is family-level recall** (cd-hit + BLAST, RepeatModeler2/EDTA
  style). Per-instance bp recall is a noisy proxy and understates library quality;
  don't use it as the primary number.
- **Parallel stages must be result-deterministic — never key a decision on a
  union-find root, rank, scheduling order, or pointer.** `refine_merge_families`
  applies merge pairs to union-find in worker-stolen order, so a component's uf
  *root* is nondeterministic; only its *partition* (the set of members) is stable.
  Two bugs once leaked that into output and made the merged consensus vary
  run-to-run at genome scale (both fixed):
  1. **Representative tie-break.** The representative was seeded at
     `representative[root] = root` with a strict `>`, so an instance-count tie
     kept the order-dependent root. It now breaks ties by **lowest family index**
     (iterate `i` ascending from a `-1` sentinel) — root-independent.
  2. **Re-refine gate.** Re-refinement of a merged family was gated on
     `i == root`, so whether a family's consensus got rebuilt depended on whether
     its representative happened to be the root. It now re-refines **every
     representative that absorbed others**, independent of the root.
  Rule for any new parallel code: the result must depend only on the work
  *content*, never on the order it was scheduled or on any root/rank/pointer that
  order produces. The whole pipeline is now byte-identical run-to-run and
  thread-count-invariant; verify with `-trace-dir` stage dumps diffed across two
  runs (`01_discover` … `08_coalesce` must all match).

---

## 15. Known doc drift (this doc / code are correct)

`TECHNICAL_DOC.md` and `README.md` predate several changes. The current code
differs from them in these load-bearing ways:

1. **No iterative R-convergence loop** in selection — it is single-pass.
2. **No post-prune recovery pass** — removed (would be inert under R-independent
   costs).
3. **Per-instance cost is R-independent** — the family-count overhead term was
   dropped.
4. **Coverage is a sweep-line, not a bitmap** (in both selection and prune).
5. **Discovery hash size is dynamic** `max(16000057, 4N/l)` (odd), not a fixed
   prime.
6. **Selection gate** is exclusive-OR-standalone (with the 50 bp / 3-instance
   standalone branch), not simply `score > 0 & ≥2 instances`.
7. **Split** uses a relaxed MDL gate (≥ 0 for originally-positive families), not
   strict `split > original`.
8. **Prune** decides on the per-instance 25%-exclusive-coverage rule.
9. **coalesce-factor default is 20.0** (some inline comments still say 1.5).
10. **New stages** not in the README architecture: `-recall-rescue` (`rescue.c`)
    and opt-in `seqkit` QC (`external_qc.c` + `tool_runner.c`); `discover_mask.c`
    is the masking subsystem split out of `discover.c`.

---

## 16. Source map

```
src/
  main.c              Pipeline driver, CLI parse, sampling, chunked discovery, usage()
  types.h             gpos_t/glen_t (int64_t), DNA encoding, PADLENGTH, constants
  genome.c/.h         FASTA load, front-padding, per-record boundaries, sampling
  kmer.c/.h           Canonical k-mer count (bucket-ownership parallel), position index
  discover.c/.h       Seed-and-extend engine; l-mer table; N-seq banded extension
  discover_internal.h Shared types between discover.c and discover_mask.c
  discover_mask.c     Post-family masking subsystem (1-vs-1 banded DP)
  align.c/.h          Multi-k-mer seeding, instance recruit/realign; opt-in rmblastn short-family recruit (-recruit-short)
  candidates.c/.h     CandidateFamily / Instance / CandidateList + MDL state
  refine.c/.h         Merge / split / fragment-assembly / prune / tandem-coalesce
  mdl.c/.h            L_int (Rissanen), per-family scoring, single-pass selection
  rescue.c/.h         Optional recall-rescue second discovery pass
  external_qc.c/.h    Opt-in seqkit-stats QC policy layer
  tool_runner.c/.h    Generic fork/execv subprocess launcher with timeout
  output.c/.h         FASTA / BED6 / TSV writers
  cmd_line_opts.c/.h  Generic CLI value parsers
```
