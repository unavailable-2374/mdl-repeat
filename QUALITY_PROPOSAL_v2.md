# mdl-repeat Output Quality Enhancement Proposal — v2

**Changes from v1**: integrated R1 reviews (BIO-N1..N6 from biology, ENG-N1..N7 from engineering). Re-prioritized based on TE-class classification gating downstream work. Removed two items (H, I); rescoped four (D, E, F, J).

**Goal**: Improve mdl-repeat library output quality beyond Stage A + Stage B + Option K baseline.

**Locked baseline (after this session's work)**:
- Option K shipped (banded DP in refine, replaces column-vote)
- chr4 family-level recall: 80×50=0.823, 80×80=0.769, 90×80=0.599
- TAIR10 nuclear: 80×50=0.811, 80×80=0.767, 90×80=0.601
- K + bumped -L=30000 ablation: equivalent to K alone within noise → K alone confirmed final

**User constraint**: Quality first. Runtime, memory, implementation cost are subordinate.

---

## Tier 1: Diagnostic + classification gate (must run first)

### A+BIO-N1: Diagnose AND classify the 24 missed in-scope families

**Combined item** (R1 bio collapsed A and BIO-N1).

For each of the 24 missed in-scope (≥3 copy) families on chr4:
1. Extract sequence via `bedtools getfasta` from chr4_full_truth.bed
2. BLAST against Dfam / Repbase Arabidopsis library
3. Assign Wicker class/order/superfamily
4. Group by class. Possible categories (per Heitkam & Schmidt 2020):
   - **Solo LTRs** (300–750 bp, matching LTR sequence of a known Gypsy/Copia family)
   - **MITEs** (100–500 bp, TIR-flanked, non-autonomous; Arabidopsis: ATREP3, HELITRONY3, ATREP10D, ATREP15 each >1000 copies)
   - **TSI / pericentric tandem variants** (170–360 bp, possibly out-of-scope as tandem)
   - **TIR/Helitron** (200–1500 bp, non-autonomous DNA TEs)
   - **Full-length elements with low chr4 copy but higher genome-wide copy** (chr4-specific evaluation underestimates)

**Output**: `MISSED_FAMILIES_CLASSIFICATION.md` with the 24 families tagged by class.

This gates everything downstream — different TE classes need different fixes. Effort: 1 day.

### A.aux: Sanity stratifications also performed at this stage
- BIO-N4: Verify any of the 24 are CEN180/TSI fragments that should be removed from the recall denominator
- BIO-N6: Verify ATHILA intervals are correctly inside the recall denominator (not excluded with CEN180)
- ENG-Q2 (R1): Check `awk '/^>/{next}{print length($0)}' chr4_v2.fa | sort -n | tail -20` to see if library has consensus ≥10 kb (informs whether ALIGN_MAX_EXTENSION is binding)

### ENG-N1: Run small-genome controlled diagnostic

**Independent of A but cheap**: build a 2 Mb synthetic genome (1 Mb random Markov background + 10 copies of one known-missed cluster representative). Run mdl-repeat. If the family is found → bottleneck is genome-scale (hash collisions, seed competition at large N). If not found → bottleneck is core seed-and-extend on this TE class.

**Output**: yes/no result per representative; published in same `MISSED_FAMILIES_CLASSIFICATION.md`.

Effort: 2-3 hours after A picks targets.

---

## Tier 2: Confirmed bugs + biology-anchored fixes

### B: Fix nested_containment_fraction over-merge bug

**Confirmed by both reviewers**. Genome C: 2194 bp LINE merged into 2897 bp LINE+LTR composite (perfect prefix relationship). Both reviewers note:
- Eng (refine.c:557, 748): both parallel `merge_worker_fn` and sequential path call same `nested_containment_fraction`; fix must be applied to both.
- Eng concern (open Q3): `nested_containment_fraction` uses raw `position` without sequence_index check — multi-chromosome correctness needs verification first.
- Bio (BIO-N2 follow-up): biologically critical for plants where LTR-LINE / LTR-DNA-TE nested topologies are routine; RepeatModeler2 explicitly handles this via reciprocal filtering.

**Tightened spec**:
1. First, audit `nested_containment_fraction` (refine.c:219-243) for cross-chromosome bug — fix if present
2. Tighten the gate predicate: require longer/shorter NOT be a strict prefix/suffix (compute identity over the non-overlapping region; if = 0%, it's a structural prefix → don't merge)
3. Require ≥2 instances of the inner element to occur OUTSIDE any longer-element instance (= "solo evidence")
4. Apply the change to BOTH `merge_worker_fn` and the sequential path; refactor into a shared helper if not already

**New test required (ENG test coverage)**: synthetic with nested element 300 bp inside 900 bp parent, plus 5 solo copies of the 300 bp; under both `-threads 1` and `-threads 4`, assert two distinct families remain.

Effort: 1 day (was 0.5).

### BIO-N2: Solo-LTR / full-element separation regression test

After fixing B, build a synthetic genome containing:
- 8 copies of a full-length ATHILA-like element (~11 kb internal + 600 bp LTR each side = 12 kb total)
- 12 copies of the same 600 bp LTR sequence as solo elements

Expected output: TWO library families:
- One ~11 kb internal consensus
- One ~600 bp LTR consensus

If only one family is produced (mixed), the merge stage is over-aggressive between solo-LTR and full-element variants and needs separate fix beyond B.

**Go/no-go gate** for any work in the 500 bp – 12 kb length-bin recovery.

Effort: 0.5 day.

### J': Increase ALIGN_MAX_EXTENSION conditionally

**Bio-R1 elevated this from Tier 3 to Tier 1-aligned.**

Justification: ATHILA internal regions are ~11 kb (PMC7385966). `ALIGN_MAX_EXTENSION = 10000` (align.c:21) caps refine-stage extension per direction per iteration. With 10 iterations × 10000, total cap is 100 kb — should not bind. But the per-iteration cap might bind during the FIRST iteration when a family's instances are at 10 kb and the refine extension wants to grow past 10 kb in one shot.

**Action**:
1. First check (open question Q2): does any chr4 library consensus reach 10 kb? If max is < 10 kb across many families, ALIGN_MAX_EXTENSION is not binding and skip this item.
2. If yes (some at 10 kb), bump ALIGN_MAX_EXTENSION to 20000.
3. Re-run chr4 + TAIR10 family recall; check if 10-99 copy bin (high-confidence long elements) improves.

Note: Eng-R1 cautions HANDOFF §4 already swept `-L` to 50000 with zero F1 change. But that was per-instance F1 (now-deprecated metric). At family-level metric this MAY differ.

Effort: 0.5 day (1 line + benchmark).

---

## Tier 3: Architectural improvements (data-conditioned)

### C: Improve subfamily splitting (instrument first)

**Both reviewers strongly endorsed instrumentation before any threshold change.**

1. **Step 1 (instrument, ~10 LOC, both threading paths)**:
   - Log how many families are skipped by `REFINE_MIN_SPLIT_INSTANCES = 10` floor (refine.h:51) — Eng-R1 open Q2
   - For families that pass the floor, log why each split was rejected: bimodality < 0.20, valley not deep enough, n_lo/n_hi < 3, MDL gain ≤ 0
   - Run on chr4 v2 + chr4_K libraries

2. **Step 2 (analyze)**:
   - If most are skipped by the 10-instance floor: lower it to 5 (most plant TE subfamilies are still detectable at 5 copies)
   - If most fail bimodality: this is biologically expected for young Copia-class families (ATCOPIA78/93). Skip threshold change.
   - If most fail MDL gain: subtle issue with split-vs-orig MDL accounting (refine.c:1303); needs separate investigation
   - **CAUTION (Eng-R1)**: do NOT change `REFINE_BIMODALITY_THRESH = 0.20` in isolation — there's a paired check at refine.c:1056 that skips valley check when bimodality ≥ 0.40; lowering the lower threshold without revisiting the upper inverts the safeguard

3. **Step 3 (fix only if instrumentation reveals a clear cause)**:
   - Otsu's reliance on `inst->divergence` (substitution-only, no gap signal — Eng-R1) is a known limitation. Splits that differ by indel rate appear unimodal here.
   - If indel-mode bimodality is the issue, replace Otsu input with `(inst->num_edits / inst->aligned_length)` which includes gaps.

**New test required**: synthetic family with two subfamilies at clearly-separated divergence (5% vs 20%); assert n_splits=1, both halves accepted; under both threading modes.

Effort: 1-2 days for instrumentation + analysis; 1-3 more days if MDL-gate fix needed.

### D' (rescoped per BIO-R1): Multi-pass with FULL-genome seeding

**Reframed**: Original v1 said "seed on residual genome" — Bio-R1 correctly noted residual-masked genome creates l-mer frequency artifacts at masked-unmasked boundaries. RepeatModeler2 does NOT do residual seeding; it re-seeds on full genome each pass and uses growing library only as acceptance filter.

**Revised design**:
1. Pass 1: standard discovery → Library 1
2. Pass 2: discovery on FULL genome (NOT masked), but at the merge stage, exclude any candidate that is ≥80-80-80 contained in Library 1 (BLAST screen)
3. Combined library: Library 1 ∪ (Pass 2 candidates surviving the screen)

This relies on the standard merge logic to deduplicate, not on mask state. No mask-export API needed (Eng-R1 ENG-N1 concern partially addressed — still need the controlled diagnostic to verify multi-pass is the right hypothesis).

Effort: 2-3 days (was 1).

**Pre-condition**: ENG-N1 controlled diagnostic must come back showing seed-competition is a real bottleneck. Otherwise this work has unknown payoff.

### F' (rescoped per BIO-R1): BLAST-based instance recruitment for SHORT elements

**Reframed**: Bio-R1 noted the design scope is "low divergence", so BLAST's divergence handling adds little. Instead, BLAST helps where banded DP fails on SHORT elements (200-500 bp): fixed bandwidth handles indels poorly in short contexts.

**Revised target**: replace `align_collect_instances` ONLY for families with consensus_length < 500. Use dc-megablast for instance discovery. Banded DP remains for ≥500 bp families.

Effort: 3-4 days (was 2-3, accounting for cons_start/cons_end mapping correctness — Eng-R1).

---

## Tier 4: Defer until after Tier 1-3 (high uncertainty)

### E'' (heavily rescoped): Structural-feature-aware seed selection

**Both reviewers warned**: profile-from-Pass-1 is circular (Bio); discover.c is fragile (Eng, HANDOFF §4 had 12 failed attempts there).

**Revised concept**: instead of k-mer profile bootstrap, use STRUCTURAL features:
- TIR pair detection (palindromic short k-mer windows)
- LTR pair detection (long-range identity within ~15 kb)
- Boost l-mer score if l-mer participates in a TIR or LTR pair

This is structurally grounded (not circular), but is a multi-week project. Defer until Tier 1-3 results are in.

### BIO-N3: TIR-seeding failure mode investigation for MITEs

If A+BIO-N1 reveals MITEs are a major missed-family class, investigate whether seed l-mers land in TIRs (palindromic short repeats at element ends). Symptom: a family's top-frequency l-mer's reverse-complement also has high frequency at similar genomic positions. Diagnostic ~0.5 day; fix (deprioritize palindromic-self-RC seeds) ~1 day if needed.

### G: Loosen merge length-ratio guard
**Both reviewers caution**: do NOT loosen below 0.70.
- Bio-R1: a 500 bp + 1000 bp merge (ratio 0.50) likely represents solo LTR + full element — biologically wrong to merge
- Eng-R1: estimate_merge_score's centroid-distance approximation is unreliable below 0.70 length ratio

**Decision**: keep at 0.70. Skip this item.

### H: Loosen Stage B fallback to 2-instance
**Both reviewers reject**:
- Bio-R1: 2-copy in pericentric region is statistically indistinguishable from spurious near-identity; FP would land in highest-uncertainty truth region (ATHILA + CEN180)
- Eng-R1: HANDOFF §0 explicitly marks 2-copy as out of design scope

**Decision**: do NOT implement.

### I: Loosen max-divergence
**Both reviewers (and HANDOFF §4) reject the global change**.
- Conditional version (only for families that already passed Stage B's standalone fallback) is acceptable per Bio-R1; would add diverged copies without global noise
- Eng-R1: HANDOFF §7 explicitly forbids without new evidence

**Decision**: defer until ENG-N1 controlled diagnostic provides new evidence.

---

## NEW Tier 1.5: Engineering robustness for big genomes (added per Eng-R1)

### ENG-N2: Fix `refine_prune_families` O(genome_len) memory allocation

`refine.c:1769` does `calloc((size_t)genome_len, 1)`. On Arabidopsis (119 MB) this is fine; on maize (2.3 GB) it's 2.3 GB single allocation; on **wheat (17 GB) it's fatal**.

**Fix**: replace per-base coverage byte array with sorted instance intervals + sweep line — O(num_instances) memory.

This is a **correctness blocker for the universal-tool stated goal**. Effort: 0.5 day.

### ENG-N3: Make `HASH_SIZE` dynamic

`HASH_SIZE = 16000057` is a `#define` in `discover_internal.h:22`. For 119 Mb genome with l=14, expected distinct l-mers up to 4^14 = 268M, 17x the table size. Forces avg chain length ~17, O(N × chain) counting — the HANDOFF §6 bottleneck.

Stage A's parallel hash addressed this in `kmer.c` for the refine stage but `discover.c::build_headptr_internal` (NOT replaced — only `build_headptr_parallel` was added) still uses the fixed table for some code paths.

**Action** (per HANDOFF §9 sketch): set `HASH_SIZE` dynamically as `max(16M, 4*N/l)`. Audit every `for (h = 0; h < HASH_SIZE; h++)` loop in discover.c and discover_mask.c to use the dynamic value.

Effort: 0.5 day.

### ENG-N4: Document/quantify parallelism non-determinism

`merge_worker_fn` uses `__atomic_fetch_add` (refine.c:445) for work distribution; thread interleaving causes pair-collection ordering variance. Two runs at `-threads 4` may differ by up to 0.5-1 pp in family-level recall.

**Action**: run baseline benchmark TWICE with same seed, document variance. If variance > 0.5 pp, single-pp claims (G/I/J magnitude changes) cannot be evaluated reliably.

Effort: 1 hour. Run before/after every other change.

### ENG-N7: Expose `-coalesce-factor` as CLI option

Currently hard-coded at `coalesce_factor = 20.0` in main.c:1170. Compact genomes (fungal, protozoan) where adjacent same-family copies are NOT in tandem arrays would have spurious mega-intervals.

**Action**: add `-coalesce-factor <float>` CLI option (default 20.0; 0 = disabled). Document in usage string.

Effort: 1 hour.

---

## NEW: Library-utility correctness (added per BIO-N5)

### BIO-N5: Verify canonical strand orientation for LTR retrotransposon consensuses

For Gypsy/Copia retrotransposon families, the GAG-POL-INT domain order should run 5' → 3' in the consensus output. If consensus is built in arbitrary orientation (whichever strand the seed l-mer happened to be on), downstream TE age estimation, K-divergence dating, and CpG decay analysis produce incorrect results.

**Action**: 
1. For each library consensus ≥ 5 kb, BLAST against UniProt for retrotransposon protein domains (RT, RH, IN)
2. If domain order is INT-RH-RT (5'→3'), the consensus is reverse-complemented; flip to canonical
3. Add this as a post-MDL polishing step

This does not affect family-level recall but materially affects library utility for downstream analyses.

Effort: 2-3 days (requires UniProt domain DB integration).

---

## Recommended execution order (revised)

**Phase 1 — Classification & diagnosis (1-2 days, MUST be first)**:
1. **A + BIO-N1 + BIO-N4 + BIO-N6 + Q2** (1 day): classify 24 missed families, sanity-check truth purity, check library max consensus length, determine ATHILA scope
2. **ENG-N1** (3 hours): controlled small-genome diagnostic for whichever class A reveals as dominant
3. **ENG-N4** (1 hour): quantify parallelism noise floor

Output: classification report + diagnostic verdict. **All subsequent decisions branch on these results.**

**Phase 2 — Confirmed bugs (1.5 days)**:
4. **B + ENG-Q3 audit** (1 day): fix nested_containment_fraction over-merge + multi-chr correctness
5. **BIO-N2** (0.5 day): synthetic ATHILA solo-LTR test (gates downstream)

**Phase 3 — Big-genome robustness (1 day, can parallel with Phase 2)**:
6. **ENG-N2** (0.5 day): fix refine_prune_families memory
7. **ENG-N3** (0.5 day): dynamic HASH_SIZE
8. **ENG-N7** (1 hour): expose -coalesce-factor

**Phase 4 — Conditional fixes (3-5 days, gated on Phase 1 results)**:
9. **J'** (0.5 day): conditional ALIGN_MAX_EXTENSION bump if Q2 shows binding
10. **C** (1-2 days instrumentation; 1-3 days fix if warranted)
11. **BIO-N3** (1.5 days): MITE TIR-seeding fix if A reveals MITEs as major class
12. **D'** (2-3 days): full-genome multi-pass with library-screen, only if ENG-N1 supports it

**Phase 5 — Polishing (3-4 days)**:
13. **F'** (3-4 days): BLAST-based instance recruitment for short elements
14. **BIO-N5** (2-3 days): canonical strand orientation polish

**Skipped per R1 review**:
- E (defer to E'' structural-feature variant, deferred indefinitely)
- G (don't loosen below 0.70)
- H (don't implement)
- I (don't redo blanket; conditional only with new evidence)

---

## Open questions persisting after R1

1. **What TE class are the 24 missed families?** (drives entire roadmap; answered by A+BIO-N1)
2. **Does any chr4 library consensus reach 10 kb?** (drives J'; answered by 1-line awk)
3. **Are ATHILA intervals counted in recall denominator or excluded with CEN180?** (affects baseline numbers)
4. **What does multipass.sh actually do** — residual-masked seeding or full-genome with screen? (affects D' implementation reference)
5. **Run-to-run variance with `-threads 4`?** (sets noise floor for evaluating any fix; ENG-N4)
6. **Is `nested_containment_fraction` cross-chromosome correct?** (audit needed before B's fix)

---

## Test coverage requirements summary (per Eng-R1)

| Item | New test required |
|---|---|
| A | Classification report file (data analysis) |
| B | Nested element synthetic; both threading modes |
| BIO-N2 | ATHILA solo-LTR synthetic; assert 2 families |
| C | Bimodal subfamily synthetic; both threading modes |
| D' | Pass 2 recovery test (low-copy starved by Pass 1 high-copy) |
| F' | Short-element instance BED coordinate test (0-based, PADLENGTH offset) |
| ENG-N2 | malloc cap regression test |
| ENG-N3 | full 7/7 + chr4 smoke after dynamic HASH_SIZE |
| ENG-N7 | -coalesce-factor 0 vs default behavior diff |
