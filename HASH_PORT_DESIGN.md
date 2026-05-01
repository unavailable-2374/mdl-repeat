# HASH_PORT_DESIGN.md

## Port Design: Replace discover.c's Linked-List Hash with kmer.c's Striped-Lock Hash

**Scope:** Drop-in replacement for the l-mer counting and seed-lookup
structures in `src/discover.c` and `src/discover_mask.c`.
**Status of this document:** Code audit only.  No source files were
modified.  Line numbers reference the codebase as of the audit session.

---

## 1. Discover-Side Function Inventory

Every function that touches l-mer counting, hash insertion/lookup, or
the mask subsystem is listed below.

### 1.1 Hash structure definitions (`src/discover_internal.h`)

| Symbol | Lines | Purpose |
|--------|-------|---------|
| `HASH_SIZE` | 22 | Prime `16,000,057` — fixed bucket array size |
| `SMALLHASH_SIZE` | 23 | 5003-bucket mini-hash used by mask |
| `SMALLL` | 24 | 6 — short l-mer prefix for the mini-hash |
| `struct posllist` | 45–48 | Singly-linked list node: one genome position |
| `struct llist` | 50–57 | Hash chain node: freq counter, last-plus/minus-occ, position list, next pointer |
| `struct repeatllist` | 59–63 | Mini-hash node used only inside `mask_headptr` |

`headptr` is declared at call sites as `struct llist **headptr` with
`HASH_SIZE` buckets.  The bucket array is a plain `malloc` of
`HASH_SIZE * sizeof(*headptr)`.

---

### 1.2 Functions in `src/discover.c`

#### `hash_function` — lines 39–54
**Purpose:** Map an l-mer (as a `char[]` of 0-3 values) to a bucket
index in `[0, HASH_SIZE)`.  
**Hash operation:** Two polynomial hashes (forward and RC), returns the
larger; returns -1 if any base is `DNA_N`.  
**Called by:** `build_headptr_internal`, `build_headptr_from_freq`,
`build_all_pos`, `find_besttmp`, `trim_headptr`, `dump_freq_table`, and
`mask_headptr` (in `discover_mask.c`).

**Key property:** The hash is *symmetric* — it returns the same bucket
for an l-mer and its reverse complement.  This is used pervasively; all
lookups rely on `lmermatcheither` after bucket lookup.

---

#### `smallhash_function` — lines 56–67
**Purpose:** 5003-bucket polynomial hash over the first `SMALLL=6`
bases of an l-mer.  Used only in `mask_headptr` to build a small
consensus hash.  
**Called by:** `mask_headptr` only.

---

#### `lmermatch` / `lmermatchrc` / `lmermatcheither` — lines 73–92
**Purpose:** Three byte-by-byte string comparators for forward, RC, and
either-strand matching.  
**Called by:** Every hash lookup site (build_headptr_internal,
build_all_pos, mask_headptr, build_pos, etc.).

---

#### `build_headptr_from_freq` — lines 196–262
**Purpose:** Read a RepeatScout `.freq` file and populate `headptr`.  
**Hash operation:** `hash_function` → chain insert; skips duplicates
via `lmermatcheither`; applies entropy and periodicity filters.  
**Called by:** `discover_families` when `params->freq_file != NULL`.

---

#### `build_headptr_internal` — lines 268–324
**Purpose:** Single-threaded O(N) genome scan to count l-mer
frequencies.  **This is the bottleneck.**  
**Hash operation:** For each genome position, computes `hash_function`,
walks chain looking for forward (`lmermatch`) or RC (`lmermatchrc`)
match; increments `freq` with TANDEMDIST guard on `lastplusocc` /
`lastminusocc`; if no match, allocates new `struct llist` node via
`malloc`.  Applies entropy + periodicity filters only at *insertion
time* (new entries only).  
**Called by:** `discover_families` when no freq file.

RISK: The loop bounds are `x = PADLENGTH .. C->orig_length - C->l`
(line 282–283).  It scans the *raw* genome only (not the doubled copy).
The doubled genome is only used by `build_all_pos` (which fills
positions for extend/mask).  Any replacement must maintain this same
scan range.

---

#### `trim_headptr` — lines 330–356
**Purpose:** Remove entries with `freq < MINTHRESH`; reset `freq = 0`
for survivors (prep for `build_all_pos`).  
**Hash operation:** Iterates all `HASH_SIZE` buckets, unlinks and frees
below-threshold nodes.  
**Called by:** `discover_families` after `build_headptr_internal` or
`build_headptr_from_freq`.

---

#### `build_all_pos` — lines 363–454
**Purpose:** Second genome scan (on the doubled genome) to populate
`struct posllist *pos` chains.  Also applies TANDEMDIST filtering by
marking within-distance same-strand pairs as negative and pruning them
from the chain; sets `C->removed[pos] = CLAIM_PERMANENT` for pruned
positions.  
**Hash operation:**
- Pass 1 (lines 374–401): `hash_function` + `lmermatcheither` → append
  `struct posllist` node to matching `llist`.
- Pass 2 (lines 406–451): Mark/remove same-strand within-TANDEMDIST
  pairs.  For removed positions sets `C->removed[-pos] =
  CLAIM_PERMANENT`.  
**Called by:** `discover_families` after `trim_headptr`.

RISK: The scan range on line 374 is `x = C->l - 1 .. C->length - C->l + 1`
i.e. it covers the *doubled* genome (length = 2*orig + PADLENGTH).
This is intentional: having RC occurrences in the position list doubles
the number of seeds available for extension.  Any replacement must
also populate positions from the doubled copy.

---

#### `find_besttmp` — lines 461–495
**Purpose:** O(HASH_SIZE) linear scan to find the entry with the
highest `freq`.  Caches `prevbesthash` / `prevbestfreq` to avoid full
rescans when the best does not change.  
**Hash operation:** Sequential walk of all buckets and chains.  
**Called by:** Main loop in `discover_families` on every iteration.

RISK: The caching optimisation (`prevbesthash`, `prevbestfreq`) provides
O(1) amortised best-pick when the top frequency does not drop.  With a
parallel hash, the cached hash index is meaningless unless we preserve
a sorted top-K heap.

---

#### `build_pos` — lines 502–538
**Purpose:** Given `besttmp`, populate `C->pos[]` / `C->rev[]` / `C->upperBoundI[]`
from the stored `struct posllist *pos` chain.  Checks `C->removed` for
claim count.  
**Hash operation:** None — reads from `llist.pos` which was built by
`build_all_pos`.

---

#### `dump_freq_table` — lines 1269–1321
**Purpose:** Write `.freq` file (optional, for reuse).  
**Hash operation:** Iterates all buckets; collects (freq, occ) pairs;
sorts by frequency; writes.  
**Called by:** `discover_families` if `params->freq_output != NULL`.

---

#### `free_headptr` — lines 1070–1094
**Purpose:** Free `headptr` and all chained `llist` + `posllist` nodes.  
**Called by:** End of `discover_families`.

---

### 1.3 Functions in `src/discover_mask.c`

#### `mask_headptr` — lines 263–405 (entry point; see internal helpers below)
**Purpose:** Mask the most-recently-found family (masters[R-1]) out of
the genome and decrement `llist.freq` for every covered l-mer.  
**Hash operations performed (complex — see below):**

1. **Build consensus mini-hash** (lines 280–310): `smallhash_function`
   on each l-mer window within the consensus.  Allocates
   `struct repeatllist` nodes.

2. **For each consensus l-mer, look it up in main headptr** (lines
   313–326): `hash_function` + `lmermatch` or `lmermatchrc` —
   O(chain) lookup.

3. **Walk each matching `llist.pos` chain** (lines 330–389): For each
   unclaimed position, run banded-DP extension to get exact mask
   boundaries, then:
   - For every genome position `y` in `[startmask, endmask]`:
     - If `C->removed[y] == 0` (first claim): call `hash_function(C,
       C->sequence + y)` and walk chain with `lmermatcheither` to find
       the l-mer entry; **decrement its `freq` by 1**.
     - Set `C->removed[y] = prev_claim + 1`.

**Called by:** `discover_families` main loop after each successful
extension.

RISK: The per-position decrement inside `mask_headptr` is the
**most performance-critical write path** in the mask subsystem.  For a
100 kb element with 1000 copies, it performs ~100M hash lookups
(100k positions × 1000 occurrences).  Replacing this with a parallel
hash that requires a lock per position would be extremely expensive.

---

## 2. kmer.c Striped-Lock Hash: Internals

### 2.1 Data structures

```
KmerTable {
    KmerEntry **buckets;    // array of size table_size (prime)
    int          k;
    size_t       table_size; // scaled to genome length: max(genome/4, 16M)
    int64_t      num_entries;
    int64_t      num_frequent;
}

KmerEntry {
    uint64_t  kmer;           // canonical (min of fwd, RC), 2-bit packed
    freq_t    frequency;      // int32_t
    gpos_t    last_plus_occ;  // last fwd occurrence (TANDEMDIST guard)
    gpos_t    last_minus_occ; // last RC occurrence (TANDEMDIST guard)
    gpos_t   *positions;      // array, allocated after kmer_build_positions()
    int32_t   num_positions;
    int32_t   cap_positions;  // capped at KMER_MAX_POSITIONS = 50,000
    KmerEntry *next;          // chain within bucket
}
```

KmerEntry nodes are allocated from `KmerPool` (pool of 4096-node
blocks).  There is a **global** pool (`g_kmer_pool`) and per-thread
pools that are merged into global after parallel counting.

### 2.2 Concurrency model

**Counting phase (kmer_count, lines 311–391):**
- `NUM_STRIPES = 4096` mutexes.
- Each thread is assigned a contiguous genome chunk.
- Before modifying a bucket chain, the thread acquires
  `stripe_locks[h % 4096]`.
- New entries are allocated from a **per-thread** `KmerPool`; after all
  threads finish, pools are merged into the global pool.

**Position-building phase (kmer_build_positions):**
- Phase 1 (parallel): atomic increment of `entry->num_positions`
  (`__atomic_fetch_add`, RELAXED).
- Phase 2 (sequential): fill positions in genome-coordinate order for
  determinism.

### 2.3 Available APIs

| Function | Purpose |
|----------|---------|
| `kmer_count(g, k, tandemdist, nthreads)` | Count; returns `KmerTable*` |
| `kmer_trim(kt, min_freq)` | Remove entries below threshold |
| `kmer_build_positions(kt, g, nthreads)` | Populate positions[] arrays |
| `kmer_lookup(kt, canonical_kmer)` | O(chain) lookup by packed uint64 |
| `kmer_free(kt)` | Teardown |
| `kmer_print_stats(kt)` | Diagnostic |

**What is NOT available (gaps vs. discover's needs):**

| Missing capability | Needed by | Notes |
|--------------------|-----------|-------|
| `freq` decrement by 1 | `mask_headptr` (line 379) | No `kmer_decrement` API |
| Iterate all entries in bucket order for `find_besttmp` | `discover_families` main loop | Only `kmer_lookup` by exact key |
| Access raw `char[]` sequence representation of a k-mer | All lmermatch* callers | kmer.c stores packed `uint64_t`; discover uses `char[]` (0–3) |
| Strand-agnostic lookup returning the stored `lastocc` position (as `char*` into genome) | mask, build_pos | kmer.c has `last_plus_occ` / `last_minus_occ` but no pointer into genome char array |
| Per-position strand tag stored inside positions[] | `build_pos` (C->rev[N]) | kmer.c positions[] stores negative = RC, positive = fwd — this DOES match, see §6.1 |
| Top-K priority extraction for `find_besttmp` | `discover_families` main loop | kmer.c has no such structure |
| Preserve entropy / periodicity-filtered entries as absent | discover filters at insert time | kmer.c has no entropy or periodicity filter |
| Support l > 31 | discover: `l` auto = ceil(1+log4(N)), on large genomes can reach 14–16 (well under 31) | Not a concern in practice |
| Support the "doubled genome" coordinate system | `build_all_pos`, all position lookups | kmer.c only ever sees the flat genome |

---

## 3. Strand Handling and Encoding Compatibility

### discover.c strand model
- Genome is doubled: `[fwd | PADLENGTH | RC_copy]`.
- `build_headptr_internal` scans only the forward copy
  (lines 282–283: `x <= C->orig_length - C->l`).
- Each `llist` entry stores `lastplusocc` (last fwd hit) and
  `lastminusocc` (last RC hit).  TANDEMDIST is applied separately
  to each strand.
- The `hash_function` returns `max(hash(fwd), hash(RC))` — same bucket
  for a sequence and its RC.  `lmermatcheither` resolves exact strand
  at lookup time.
- `C->rev[n]` in `build_pos` is set by comparing the occurrence against
  `besttmp->lastocc` with `lmermatch`; 0 = forward, 1 = RC.

### kmer.c strand model
- No genome doubling.  Scans original genome only.
- `kmer_canonical(fwd, rc) = min(fwd, rc)` — same canonical form as
  discover's `max(hash(fwd), hash(rc))` bucket rule in spirit, but the
  actual canonical key is different (kmer.c: min packed uint64; 
  discover: max polynomial hash).
- `positions[]` stores `-(gpos_t)i` for RC occurrences, positive for
  forward.  This matches what `build_pos` needs for `C->rev[]`.
- No PADLENGTH prefix in kmer.c's Genome (kmer.c starts at
  `PADLENGTH` per `kmer.c:248` — it skips the pad).

RISK: The doubled-genome architecture is deeply embedded in discover's
position list (`build_all_pos` line 374: scans to `C->length` = 2*N +
PADLENGTH, including the RC copy).  kmer.c never builds a doubled
genome.  Any adapter that feeds kmer.c positions to the discover loop
must synthesise the RC occurrence set, either by running a second pass
or by storing both positions explicitly.

---

## 4. Port Design Options

---

### Option α — Thin Adapter Layer (thinnest possible)

**Concept:** Keep all of discover.c's internal structures (`struct
llist`, `struct posllist`, `headptr`, `HASH_SIZE`) unchanged.
Introduce a new function `build_headptr_parallel(C, headptr)` that
replaces only `build_headptr_internal`.  It uses `kmer_count`
(parallel) to get frequencies, then translates `KmerEntry` records back
into `struct llist` nodes and inserts them into `headptr`.

**Mechanism:**
1. Call `kmer_count(doubled_genome, l, TANDEMDIST, nthreads)` where
   `doubled_genome` is a synthetic `Genome` wrapping `C->sequence[0..orig_length-1]`
   (not the doubled copy, matching existing scan bounds).
2. Iterate all `KmerEntry` in `kt` with `freq >= MINTHRESH`:
   - Decode packed `uint64_t` back to `char[l]` (0–3 encoding).
   - Apply entropy filter and periodicity filter.
   - Compute `hash_function(C, decoded_lmer)`.
   - Allocate `struct llist` node; copy `frequency`, `lastocc` (from
     `kmer_entry->last_plus_occ` or `last_minus_occ`).
   - Insert into `headptr[h]`.
3. Call `kmer_free(kt)`.  From this point, all existing code is
   unchanged.

**LOC delta:** ~100 lines added (new function).  Zero lines modified in
existing code paths.

**Tests at risk:** None — the rest of the pipeline is untouched.
Synthetic test results (7/7) should be stable; border risk of
frequency value differences due to TANDEMDIST strand logic (see §6.3).

**Invariants at risk:**
- TANDEMDIST applied in kmer.c is per-canonical-strand (RC treated as
  separate).  discover's TANDEMDIST is also per-strand.  However,
  kmer.c counts from position 0 of the flat genome (after PADLENGTH),
  whereas discover's internal counter starts from `PADLENGTH`.  These
  are the same if the input `Genome` struct carries `PADLENGTH`
  properly — verify `genome->sequence[0..PADLENGTH-1]` are DNA_N.
- kmer.c's TANDEMDIST guard uses `last_plus_occ` initialised to -1000000
  (line 197), identical to discover's `lastminusocc = -1000000` initial
  value — consistent.
- **Entropy and periodicity filters** must be re-applied by the adapter
  because kmer.c does not know about them.  This requires decoding each
  packed k-mer back to `char[]` for `compute_entropy` and
  `is_periodic_lmer`.

**Parallelism gain:** Only the counting phase becomes parallel.
`trim_headptr`, `build_all_pos`, `find_besttmp`, `mask_headptr` remain
serial.

**Estimated speedup on chr4 (45 Mb):** 4–8× on counting phase (which
is currently ~7–8 min); total wall-clock improvement depends on the
fraction of time spent in counting vs. extension (likely >50%).

---

### Option β — Rewrite discover to use kmer.c's API directly

**Concept:** Replace `struct llist **headptr` throughout discover.c and
discover_mask.c with a `KmerTable*`.  Adapt all lookup, frequency
decrement, position iteration, and best-seed selection to use
`kmer_lookup` and direct `KmerEntry*` access.

**Mechanism:**
- Remove `HASH_SIZE`, `struct llist`, `headptr` malloc/free, all
  hash_function / lmermatch* calls from counting/lookup paths.
- Add `kmer_entry_decrement_freq(KmerEntry *e)` to kmer.c (1 line).
- Replace `find_besttmp` scan with a max-heap (priority queue) built
  over KmerEntry* after trimming.
- Replace `build_all_pos`'s position chain with kmer.c's
  `positions[]` array.
- Handle the doubled-genome RC occurrences by running
  `kmer_build_positions` on the doubled genome struct (a new struct
  wrapping the full `C->sequence[0..length-1]`).
- Replace `mask_headptr`'s per-position `hash_function` + chain walk
  with `kmer_lookup(kt, kmer_pack(seq+y, l))` + direct
  `entry->frequency--`.

**LOC delta:** ~500–800 lines modified/deleted, ~200 added.

**Tests at risk:** All 7 synthetic tests and the smoke test baseline
are at risk because the algorithm's exact integer behaviour would
change:
- Heap-based best-seed selection changes tie-breaking order vs.
  the current linear scan.  This changes *which* seed is chosen at
  each round, potentially changing the number and content of families.
- Determinism is broken unless the heap has a stable secondary sort key.

**Invariants at risk (RISK):** The entire `mask_headptr` subsystem
relies on `llist.freq` being decremented in situ.  Replacing this with
`KmerEntry.frequency` decrement is semantically equivalent but requires
ensuring that the `mask_headptr` main loop (which scans `posllist`
chains, not `KmerEntry.positions[]`) is rebuilt to use the kmer.c
positions array.  The doubled-genome position set and the
`CLAIM_PERMANENT` semantics in `C->removed[]` must be preserved
exactly.

---

### Option γ — Parallel Counting Only, Merge Back, Keep Rest Serial

**Concept:** Same as Option α, but even thinner: spawn `kmer_count`
with nthreads, then copy the counted frequencies directly into the
existing `headptr` structure.  Do not attempt to reuse KmerEntry's
`positions[]`.  This is a strict subset of α — the key difference is
that the `kmer_count` call must wrap the flat (non-doubled) genome
exactly matching `build_headptr_internal`'s scan range.

This is the "frequency boost only" variant: addresses the bottleneck
(counting is the slow O(N) serial scan), leaves the rest of the
pipeline identical.

**LOC delta:** ~80 lines added; zero lines modified.

**Tests at risk:** Same as α (low).

**Invariants at risk:** Same as α but smaller surface area.  The main
risk is the TANDEMDIST / strand frequency counting discrepancy
described in §6.3.

---

### Option Comparison Table

| Criterion | α (adapter) | β (full rewrite) | γ (count-only) |
|-----------|-------------|-----------------|----------------|
| LOC delta | +100 | ±700 | +80 |
| Existing logic changed | No | Yes (heavily) | No |
| Mask subsystem risk | None | High | None |
| Determinism risk | Low | High | Low |
| Tests at risk | 7/7 (but low P) | 7/7 (high P) | 7/7 (but low P) |
| Parallelism gain | Counting only | Counting + lookup | Counting only |
| Implementation effort | 1–2 days | 5–10 days | 0.5 day |
| Correctness confidence | High | Medium | High |

---

## 5. Recommendation

**Recommended: Option γ (count-only parallel) as Phase 1; Option α as
Phase 2 if γ is insufficient.**

Rationale:

1. The HANDOFF.md engineering blocker states that the l-mer counting
   phase (the single-thread O(N) scan in `build_headptr_internal`) is
   the observed bottleneck — chr4 takes ~8 min, TAIR10 119 Mb did not
   finish in 16 min.  Option γ directly and conservatively attacks that
   specific bottleneck.

2. The mask subsystem (`mask_headptr` in `discover_mask.c`) is the
   highest-complexity invariant in the system and has been empirically
   tuned (the comment at lines 359–365 documents that a seemingly minor
   change to the decrement behaviour dropped recall from 0.389 to 0.362).
   Options β and even α are unnecessary risks to this subsystem.

3. Option α's decode-back-to-char[] step (unpacking uint64 → char[l])
   adds no algorithmic value for Phase 1; γ avoids it by calling the
   existing `lmermatch*` + `hash_function` infrastructure.

4. Option β should be deferred until there is a stable regression
   baseline from the Milestone 4 BED-PR test suite (17/17 passing) and
   until the discover subsystem has end-to-end integration tests at the
   frequency-count level.

---

## 6. Precise Change List for the Implementer (Option γ)

### File: `src/discover.c`

**Change 1 — add helper include at top**  
After line 22 (`#include "discover_internal.h"`), add:
```c
#include "kmer.h"
```

**Change 2 — add a new function `build_headptr_parallel`**  
Insert after `build_headptr_internal` (currently ending at line 324),
before `trim_headptr` (line 330).

The function signature must be:
```c
static void build_headptr_parallel(DiscoverContext *C,
                                   struct llist **headptr,
                                   int num_threads);
```

Internal logic:
1. Build a temporary `Genome` struct (stack-allocated or heap) that
   points to `C->sequence` (which at this point is the raw genome,
   not yet doubled — this matches `build_headptr_internal`'s scan range
   of `PADLENGTH .. C->orig_length - C->l`).
   Set `g_tmp.sequence = C->sequence_owned` (the original, pre-doubled
   allocation), `g_tmp.length = C->orig_length`,
   `g_tmp.raw_length = C->orig_length - PADLENGTH`,
   `g_tmp.num_sequences = C->disc_num_sequences`,
   `g_tmp.boundaries = C->disc_boundaries`.

   RISK: `C->sequence_owned` at call time points to the doubled genome
   (see `discover_families` line 1375: `C->sequence = C->sequence_owned`
   where `C->sequence_owned` is the result of `build_doubled_genome`).
   The scan bound in `build_headptr_internal` is
   `x <= count_end = C->orig_length - C->l` (line 282).  The kmer.c
   call must therefore use a genome of length `C->orig_length` (not
   `C->length`).  Construct `g_tmp.length = C->orig_length` and pass
   that as the Genome to `kmer_count`.

2. Call `kmer_count(&g_tmp, C->l, C->TANDEMDIST, num_threads)`.

3. Iterate `kt->buckets[0..kt->table_size-1]`.  For each `KmerEntry *e`:
   - Skip if `e->frequency < C->MINTHRESH` (will be handled by
     `trim_headptr`, but skip early to avoid decoding).
   - Decode `e->kmer` (packed uint64, canonical = min(fwd, RC)) back
     to `char[l]` in 0–3 encoding.  Use the formula:
     ```c
     for (int i = l-1; i >= 0; i--) {
         decoded[i] = e->kmer & 3;
         e->kmer >>= 2;
     }
     ```
     (copy `e->kmer` to a local variable first; do not modify the entry.)
   - Apply `compute_entropy(C, decoded)` — skip if > `C->MAXENTROPY`.
   - Apply `is_periodic_lmer(decoded, C->l)` — skip if periodic.
   - Compute `h = hash_function(C, decoded)`.  If `h < 0` skip.
   - Check if already in `headptr[h]` via `lmermatcheither` chain
     walk.  (Should not occur if kmer.c deduplicates; defensive check.)
   - Allocate `struct llist *node = malloc(sizeof(*node))`.
   - Set `node->freq = e->frequency`.
   - Determine `node->lastocc`:  prefer `e->last_plus_occ` if >= 0,
     else `e->last_minus_occ`.
   - Set `node->lastplusocc = e->last_plus_occ`.
   - Set `node->lastminusocc = e->last_minus_occ`.
   - Set `node->pos = NULL`, `node->next = headptr[h]`.
   - `headptr[h] = node`.

4. `kmer_free(kt)`.  (Note: `kmer_free` calls `pool_free_all()` which
   resets the **global** pool.  If `kmer_count` is subsequently called
   again — e.g. in `main.c` Step 3 — the global pool is reset cleanly.
   There is no double-free risk because KmerEntry nodes are pool-managed
   and `kmer_free` frees the pool blocks, not individual entries.)

   RISK: `kmer_free` resets `g_kmer_pool` (a module-level global in
   kmer.c).  If `discover_families` is called more than once per
   process (e.g. in chunked discovery), each call will reset the pool.
   The pool is only used by kmer.c functions, so as long as no `KmerTable`
   is live when `kmer_free` is called, this is safe.  Verify in
   `discover_chunked` (main.c line 356) that each chunk's
   `discover_families` call returns and frees its `kt` before the next
   chunk begins.

**Change 3 — add `num_threads` parameter to `discover_families`**  
`src/discover.h` line 57–58: change the signature to:
```c
CandidateList *discover_families(const Genome *genome,
                                 const DiscoverParams *params,
                                 int num_threads);
```

Add `int num_threads` to `DiscoverParams` (or pass it as a separate
argument — the separate argument is preferred to avoid changing
`DiscoverParams` serialisation or the `.freq` file reading path).

**Change 4 — wire `build_headptr_parallel` into `discover_families`**  
At `discover_families` lines 1418–1424:
```c
    if (params->freq_file) {
        ...
        build_headptr_from_freq(C, headptr, params->freq_file);
    } else {
        ...
        build_headptr_internal(C, headptr);   // <-- replace or branch
    }
```
Change to:
```c
    if (params->freq_file) {
        build_headptr_from_freq(C, headptr, params->freq_file);
    } else if (num_threads > 1) {
        build_headptr_parallel(C, headptr, num_threads);
    } else {
        build_headptr_internal(C, headptr);
    }
```
This preserves the single-thread path exactly, minimising regression
risk.

### File: `src/main.c`

**Change 5 — pass `num_threads` to `discover_families`**  
Line 857: `candidates = discover_families(genome, &dparams);`  
Change to: `candidates = discover_families(genome, &dparams, num_threads);`  
Line 345 (in `chunk_discover_worker`): similarly update the call.  The
`ChunkWorkerArgs` struct needs a `num_threads` field, or pass 1
explicitly for chunk workers (since chunks already run in parallel at
the chunk level, intra-chunk parallelism would over-subscribe CPUs).

Recommended: pass `1` inside `chunk_discover_worker` — the outer loop
in `discover_chunked` is already parallel across chunks.  Only pass
`num_threads` for the non-chunked path.

### File: `src/discover_internal.h`

No changes required for Option γ.

### File: `src/discover_mask.c`

No changes required.

### Build system (`Makefile`)

No changes required: `kmer.c` is already in the build (used by the
refine pipeline in main.c).

---

## 7. Gotchas and Risk Register

### 7.1 Strand handling: doubled genome vs. kmer.c flat genome

**Risk level: MEDIUM**

`build_headptr_internal` scans only the forward half of the doubled
genome (lines 282–283).  However, `build_all_pos` scans the *full*
doubled genome (lines 374–454).  The counting phase and the position
phase therefore have different scan domains.

Option γ replaces only the *counting* phase.  The position phase
(`build_all_pos`) is unchanged and still uses the doubled genome.
The critical constraint is: `kmer_count` in the adapter must be given
a genome of length `C->orig_length` (not `C->length`) to match the
existing scan bound.  PADLENGTH is still present at the start of
`C->sequence_owned` — so `g_tmp.sequence = C->sequence_owned` and
`g_tmp.length = C->orig_length` is correct.

RISK: `C->sequence_owned` points to a buffer of length `C->length = 2*orig + PADLENGTH`
allocated by `build_doubled_genome`.  The kmer.c worker will compute
`end = g_tmp.length - l + 1 = orig_length - l + 1` and scan
`[PADLENGTH, end)`.  The doubled genome's forward copy occupies
`[PADLENGTH, orig_length)` and the RC copy is beyond.  The scan bound
is correct.

### 7.2 Memory pool conflict (kmer.c global pool)

**Risk level: MEDIUM**

`kmer.c` uses a file-static `KmerPool g_kmer_pool`.  `kmer_free(kt)`
resets this pool.  If `kmer_count` is called inside `build_headptr_parallel`
and `kmer_free` is then called, the pool is reset before Step 3 of
main.c calls `kmer_count` again for the refine pipeline.

This is **safe** because the `discover_families` call completes and
frees `kt` before returning to `main.c`, which then starts Step 3.
There is no overlap.

RISK for chunked discovery: `discover_chunked` runs multiple
`discover_families` calls in parallel threads (up to `max_parallel`
simultaneous, see main.c lines 512–548).  Each thread calls `kmer_count`
→ `kmer_free` independently, all sharing `g_kmer_pool`.  This is a
**data race** on `g_kmer_pool`.

Resolution: For Option γ, when called from a chunk worker, pass
`num_threads = 1` (sequential kmer_count path, no per-thread pool).
Serial kmer_count allocates from `g_kmer_pool` but each
`chunk_discover_worker` runs serially within its own thread.  The
race is on `g_kmer_pool.head` — `pool_alloc` is called without a lock.
This remains a race if two chunks run simultaneously.

**Concrete fix for the chunked case:** Either (a) add a mutex around
`pool_alloc` / `pool_free_all`, or (b) use separate per-call pools.
The simplest fix is (b): in `build_headptr_parallel`, allocate a local
`KmerPool lpool = {NULL, POOL_BLOCK_SIZE}`, use `pool_alloc_local(&lpool)`
for all counting, and free it with a local `pool_free_local` instead of
calling `kmer_free(kt)` (which hits the global pool).

However, kmer.c's public API does not expose local pool allocation.
The cleanest fix is to add a `kmer_count_local(g, k, td, nthreads,
out_pool)` variant, or to serialize the chunked-discover path at the
counting phase.  Since chunk workers already run in parallel at the
chunk level, passing `num_threads = 1` to `build_headptr_parallel`
inside a chunk worker is correct and avoids the race.

### 7.3 TANDEMDIST semantics: frequency values may differ slightly

**Risk level: LOW-MEDIUM**

`build_headptr_internal` applies TANDEMDIST per-strand: forward uses
`lastplusocc` and RC uses `lastminusocc`.  kmer.c does the same
(lines 177–184 of kmer.c).

However, discover's canonical lookup uses `hash_function` (returns
`max(hash(fwd), hash(rc))`), which can place the same canonical k-mer in
a *different* bucket depending on which strand happened to be the first
occurrence.  kmer.c uses `min(fwd, rc)` as the canonical form.  The
resulting `freq` count should be numerically identical because the
TANDEMDIST logic only depends on position gaps, not on canonical form.

RISK: The very first occurrence's strand determines whether `lastplusocc`
or `lastminusocc` is initialised.  If the first occurrence arrives as RC
in discover but as forward in kmer.c (due to different initial sentinel
semantics or scan order), there could be a single off-by-one in frequency
for a small number of l-mers.  This would change which seeds are chosen
at high frequency, potentially altering the discovered families slightly.

This is an acceptable risk: the existing entropy/periodicity filters and
the MDL post-filter absorb minor count noise.  However, it should be
verified by running the 7/7 synthetic test suite and confirming output
is identical or statistically equivalent before committing.

### 7.4 Entropy and periodicity filters: must be re-applied

**Risk level: LOW**

`build_headptr_internal` applies `compute_entropy` and `is_periodic_lmer`
only when creating a *new* node (first occurrence).  kmer.c has no
knowledge of these filters.  The adapter in Option γ must apply them
during the translation loop.  Both functions are `static` in
`discover.c`; they need to remain in scope at the `build_headptr_parallel`
call site (which is also in `discover.c`), so no visibility changes are
needed.

The `compute_entropy` function (lines 98–112) takes a `DiscoverContext*`
and a pointer into `C->sequence`.  The adapter needs a representative
position for the decoded k-mer.  The correct approach is to use
`e->last_plus_occ` as the position (or `e->last_minus_occ` if
`last_plus_occ < PADLENGTH`), pointing into `C->sequence` at the
decoded k-mer's last forward occurrence.  This is consistent with how
`build_headptr_internal` applies the filter.

RISK: For a k-mer whose only occurrences are on the RC strand,
`last_plus_occ` is `-1000000` (sentinel).  The adapter must detect this
and use `last_minus_occ` instead.

### 7.5 Determinism: find_besttmp linear scan order

**Risk level: LOW for Option γ**

`find_besttmp` (lines 461–495) scans `headptr[0..HASH_SIZE-1]` in
order and returns the first entry with maximum frequency.  The parallel
counting in Option γ fills the same `headptr` with the same frequencies
(modulo §7.3 above), and `find_besttmp` still runs serially over the
same bucket structure.  Determinism of seed selection is preserved
provided the frequencies are identical.

RISK for Option β: replacing with a heap breaks tie-breaking order
entirely.

### 7.6 `kmer_trim` vs `trim_headptr` interaction

**Risk level: LOW**

In Option γ, the adapter does NOT call `kmer_trim` before translating
into `headptr`.  Instead it lets `trim_headptr` handle below-threshold
entries (as before), or it can skip them during translation (recommended,
to avoid unnecessary `struct llist` allocation).  Either is correct.

Do NOT call `kmer_trim` and `trim_headptr` both — `kmer_trim` only
affects the `KmerTable`; `trim_headptr` only affects `headptr`.  They
are independent structures.

### 7.7 Pool reset and kmer_free call sequence

**Risk level: LOW if serialized**

As noted in §7.2, `kmer_free` resets the global pool.  In the
non-chunked path, the sequence is:
1. `build_headptr_parallel` → `kmer_count` → `kmer_free` (pool reset).
2. `trim_headptr` / `build_all_pos` (no kmer.c involvement).
3. main.c Step 3: `kmer_count` (fresh pool allocation).

This is safe: pool is reset before Step 3 allocates into it.

### 7.8 MAXR (C->MAXR) and memory safety

**Risk level: LOW**

Unrelated to the hash port.  Noted for completeness: the adapter
allocates `struct llist` nodes with `malloc`, not from a pool.  If
the frequency table contains millions of distinct l-mers (possible for
highly repetitive genomes), this could produce many allocations.  The
existing `build_headptr_internal` has the same property — no regression.

### 7.9 `build_headptr_from_freq` path is not parallelised

**Risk level: NOT APPLICABLE**

When a `.freq` file is provided, `build_headptr_from_freq` is used
instead.  Option γ does not touch this path.  It remains single-threaded.
If the `.freq` file path is also slow, a separate optimisation would be
needed.  Currently it is not identified as a bottleneck.

---

## 8. Summary Table

| What to change | File | Lines affected | Nature |
|----------------|------|----------------|--------|
| Add `#include "kmer.h"` | `src/discover.c` | 22 | Add 1 line |
| Add `build_headptr_parallel` function | `src/discover.c` | After line 324 | New ~90 LOC |
| Branch on `num_threads` in `discover_families` | `src/discover.c` | 1418–1424 | Modify 3 lines |
| Add `num_threads` param to `discover_families` signature | `src/discover.c` + `src/discover.h` | discover.c:1327; discover.h:57 | Modify 2 lines |
| Update `discover_families` call site (non-chunked) | `src/main.c` | 857 | Modify 1 line |
| Update `discover_families` call site (chunk worker) | `src/main.c` | 345 + `ChunkWorkerArgs` | Modify 2 lines + 1 field |
| Pass `num_threads=1` inside `chunk_discover_worker` | `src/main.c` | 345 | Modify 1 line (for safety, §7.2) |

**Total: ~100 net new LOC, ~10 modified lines, 0 lines deleted.**

All other files (`discover_mask.c`, `discover_internal.h`, `kmer.c`,
`kmer.h`, `types.h`) require **no changes**.

---

*End of HASH_PORT_DESIGN.md*
