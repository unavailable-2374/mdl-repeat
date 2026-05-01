# mdl-repeat Output Quality Enhancement Proposal — v4

**Changes from v3**: integrated R3 reviews. Bio R3 = 0 new items (1st converged round). Eng R3 = 1 new item ENG-N11 (`check_instance_overlap` 5 unsafe sites) + 2 spec clarifications (ENG-N8 runtime guard alone is sufficient; ENG-N9 sort change alone is sufficient). v3's ENG-N9.aux was incomplete; v4 enumerates all 10 `.position` sites in refine.c and disposes each inline.

**Changes from v2 (carried forward)**: integrated R2 reviews. 4 bio items (BIO-N7..N10), 3 eng items (ENG-N8..N10). Q4/Q6 audits. Tier 1.5a multi-chromosome correctness theme. Stage B Fix 1 baseline clarified. BIO-N2 redesigned.

**Goal**: Improve mdl-repeat library output quality beyond Stage A + B + K baseline.

**Locked baseline (clarification per BIO-N8)**:
- Stage A (parallel l-mer hash) ✅ shipped
- **Stage B Fix 1 (MDL standalone-fallback admit gate) ✅ shipped** — accepts families with `standalone_savings > 0 AND consensus_length ≥ 50 AND num_instances ≥ 3`. The −0.241 bp_rec drop documented in REFINE_TRACE_REPORT was BEFORE this fix; current baseline already includes it.
- Stage B Fix 2 (merge length-ratio guard at 0.7) ✅ shipped
- Option K (banded DP in refine extension) ✅ shipped
- chr4 family-level recall (current baseline): 80×50=0.823, 80×80=0.769, 90×80=0.599
- TAIR10 nuclear: 80×50=0.811, 80×80=0.767, 90×80=0.601
- K + L=30000 ablation: equivalent to K alone within noise → K alone confirmed final

---

## Tier 1: Diagnostic + classification gate (must run first)

### A+BIO-N1: Diagnose and classify the 24 missed in-scope families

**As v2.** Sequence extraction → BLAST against Dfam/Repbase → Wicker class assignment for each of the 24 chr4 missed families. Branch entire downstream roadmap on the result.

Output: `MISSED_FAMILIES_CLASSIFICATION.md`. Effort: 1 day.

### A.aux: Sanity stratifications + answers to open questions
- **Q2** (single awk, ~30 sec): does any chr4 / TAIR10 library consensus reach 10 kb? Answers J' urgency.
- **BIO-N4**: are any of the 24 missed CEN180/TSI fragments? If yes, remove from denominator → real recall is higher.
- **BIO-N6**: are ATHILA intervals correctly inside recall denominator? Per HANDOFF the inside-cen recall=0.699 implies yes, but verify positionally.

### ENG-N1: Controlled small-genome diagnostic

2 Mb synthetic (1 Mb random Markov + 10 copies of one missed-class representative). Tests whether bottleneck is genome-scale or core algorithm. Effort: 3 hours after A picks targets.

### ENG-N4: Quantify parallelism noise floor

**Run before any recall-delta claims are trusted.** Two baseline runs at `-threads 4`, same genome. Measure family-level recall variance. If > 0.5pp, single-pp claims (J', G, I-conditional) cannot be evaluated reliably. Effort: 1 hour.

---

## Tier 1.5a: Multi-chromosome correctness (NEW THEME, R2)

**Systemic finding from Eng-R2**: `grep -n seq_index refine.c` returns empty. Every positional comparison in refine.c is potentially subject to cross-chromosome bugs because raw padded `.position` values are compared without sequence boundary checks. Three confirmed code paths (Q6, ENG-N8, ENG-N9); audit needed for the rest.

### B+Q6: Fix `nested_containment_fraction` cross-chromosome bug

Q6 confirms the bug: refine.c:219-243 compares raw `.position` values without `.seq_index` check. On multi-chromosome inputs, two instances from different chromosomes with proximate raw position values trigger spurious overlap.

**Fix step 1**: at refine.c:229, add `if (shorter->instances[i].seq_index != longer->instances[j].seq_index) continue;`

**Fix step 2** (the original B): tighten the gate predicate to require non-prefix/suffix relationship + ≥2 solo instances of the inner element.

**Apply to BOTH `merge_worker_fn` and sequential path** (both call `nested_containment_fraction`).

Effort: 1 day (audit + fix both paths + new synthetic test from BIO-N2 + dual-thread test).

### ENG-N8 (NEW R2): Fix `refine_coalesce_tandem_instances` cross-chromosome bug

refine.c:2487-2496 sorts instances by raw `.position` only and applies `coalesce_factor * consensus_length` (default 20× = 200kb for a 10kb element) gap threshold. **Cross-chromosome instances within PADLENGTH of each other in raw coordinates can be silently coalesced into one BED interval spanning a chromosome boundary.** Output BED has invalid local coordinates (overflow chromosome length or negative).

**Severity: BLOCKER for multi-chromosome inputs** (TAIR10 has 5 chromosomes; affects every published number).

**Fix**: add `seq_index` to the sort key; insert `if (cur->seq_index != active->seq_index) { active = cur; continue; }` before gap test at refine.c:2494.

Effort: 0.5 day (includes new dual-chromosome test).

### ENG-N9 (NEW R2): Fix `refine_assemble_fragments` cross-chromosome sweep-line

refine.c:2000-2007. Sweep-line break condition `entries[j].start - entries[i].start > D` (D=30kb default) doesn't filter by chromosome. Cross-chromosome instances within 30kb in raw coordinates are counted as co-occurring → spurious fragment assembly. MDL gate provides partial protection but not guaranteed.

**Severity: WARNING** (MDL partial protection).

**Fix**: add `seq_index` to `InstanceEntry` (refine.c:1864), populate at refine.c:1946, sort by `(seq_index, start)`, guard inner loop.

Can be folded into the same patch as ENG-N8 (shared `seq_index` infrastructure).

Effort: 0.5 day if combined with ENG-N8; 1 day standalone.

### ENG-N11 (NEW R3): Fix `check_instance_overlap` (refine.c:252-334) — 5 unsafe sites

**R3 audit pinpointed exactly which sites need fixing**. Eng-R3 ran `grep -n '\.position' src/refine.c` (10 sites total). 5 are unsafe and not yet covered by B/N8/N9: ALL inside `check_instance_overlap`.

| refine.c line | Path | Status |
|---|---|---|
| 268 | sorted, outer loop | UNSAFE |
| 275 | sorted, binary-search midpoint | UNSAFE |
| 284 | sorted, forward-scan | UNSAFE |
| 307 | linear-fallback, outer loop | UNSAFE |
| 311 | linear-fallback, inner loop | UNSAFE |

The function is called at the merge gate (refine.c:557 + 748) to decide whether to block a merge due to instance overlap. Cross-chromosome false positives block desirable merges; false negatives admit incorrect ones.

**Fix**:
- Linear-fallback inner loop (line 310): add `if (a->instances[i].seq_index != b->instances[j].seq_index) continue;` as first body statement
- Sorted path inner loop (line 283): same guard before the `b_start >= a_end` break check
- `seq_index` already exists in `Instance` struct (candidates.h:22) — no struct change needed

Effort: 15 minutes if folded into the ENG-N8 + N9 + B+Q6 combined patch.

### ENG-N9 spec clarification (NEW R3)

Eng-R3 noted v3's "guard inner loop" was ambiguous. **Resolution**: sorting `cmp_instance_entry` by `(seq_index, start)` alone is sufficient — the existing `break` at refine.c:2005 fires naturally at the chromosome boundary because `entries[j].start - entries[i].start` becomes the full chromosome offset (>> D). No in-body `seq_index` check needed. Update implementation spec accordingly.

### ENG-N8 spec clarification (NEW R3)

Eng-R3 noted that the runtime guard alone (without changing the sort key) is sufficient for ENG-N8 — adjacent same-`seq_index` instances would only be coalesced when their gap is ≤ threshold AND on the same chromosome. **Decision**: ship the runtime guard only (`if (cur->seq_index != active->seq_index)`) without changing the sort key — minimal change, no need to also rewrite `cmp_inst_by_start`.

### ENG-N9.aux (revised): Position-comparison audit COMPLETED inline

R3 enumerated all 10 `.position` sites in refine.c. Audit result:
- 5 unsafe: B+Q6 (lines 225, 230), ENG-N11 (lines 268, 275, 284, 307, 311)
- 1 covered by ENG-N8 (line 2474)
- 2 safe by design (refine.c:1776, 1834 — `refine_prune_families` uses padded global coordinates intentionally for the cov[] array)

**No further audit needed.** Effort line removed (folded into ENG-N11).

---

## Tier 1.5b: Big-genome scalability (per Eng-R1)

### ENG-N2: Fix `refine_prune_families` O(genome_len) memory

refine.c:1769 allocates `calloc(genome_len, 1)`. Wheat (17 GB) → fatal. Replace with sorted instance intervals + sweep line, O(num_instances).

**Bundle with ENG-N10 (NEW R2)**: replace the O(n²) selection sort at refine.c:1757-1766 with `qsort`. Both touch `refine_prune_families` so combined patch.

Effort: 1 day (was 0.5; +0.5 for the qsort + multi-genome scalability test).

### ENG-N3: Make `HASH_SIZE` dynamic

`#define HASH_SIZE 16000057` → `max(16M, 4*N/l)`. Audit all `for (h = 0; h < HASH_SIZE; h++)` in discover.c + discover_mask.c. SMALLHASH_SIZE confirmed out of scope.

Effort: 0.5 day.

### ENG-N7: Expose `-coalesce-factor` CLI

Currently hard-coded `coalesce_factor = 20.0` at main.c:1170. Add CLI option (default 20.0; 0=disabled). Document.

**BIO-N9 follow-up**: for nested-stack-rich genomes (maize, wheat), default 20× may be too aggressive (spurious mega-intervals across stacked LTRs). Document recommended downward adjustment (~5×) for those use cases.

Effort: 1 hour.

---

## Tier 2: Confirmed bugs + biology-anchored fixes

### BIO-N2 (REVISED per BIO-N7 + N10): Solo-LTR / full-element regression test

**v2 specified 8 full-length + 12 solo. R2 bio noted this is biologically inverted (real solo:full ratio is 3:1 to 6:1) AND too close to MDL noise floor; test may pass for wrong reason.**

**Revised spec**:
- **20 full-length** ATHILA-like elements (12 kb each = 11 kb internal + 600 bp LTR each side)
- **80 solo LTR** copies (600 bp)
- **Critical**: full-length element must be seeded by an INTERNAL-sequence l-mer (distinct from LTR sequence), otherwise pipeline produces only solo-LTR consensus and test cannot distinguish "correct two-family output" from "incomplete one-family output"
- **Assertion (per BIO-N10)**: BOTH families must be in the FINAL post-MDL library, not just pre-MDL candidates. This simultaneously validates B's gate fix AND Stage B Fix 1 (standalone-fallback admit). If only one family appears, distinguish whether it's bug B (merge collapse) or MDL select (fallback didn't fire).

**Pre-condition**: B (cross-chr fix + nested gate tighten) must land first.

Effort: 1 day (synthetic build + run + dual-threading verification).

### J': Conditional ALIGN_MAX_EXTENSION bump

**Conditional on Q2 result**:
- If chr4 max consensus < 10 kb: cap not binding, skip
- If chr4 max consensus ≥ 10 kb: bump ALIGN_MAX_EXTENSION 10000 → 20000

**BIO-N9 amendment**: even if Arabidopsis doesn't need it, **ship the bump unconditionally** for non-Arabidopsis (maize, wheat have routine 10-19 kb retrotransposons that WILL hit the cap). The "skip if Arabidopsis doesn't need it" gate from v2 is wrong when judging universal-tool requirements.

**Decision**: ship the bump unconditionally. Verify chr4 + TAIR10 don't regress.

Effort: 0.5 day (1 line + dual benchmark).

---

## Tier 3: Architectural improvements (data-conditioned)

### C: Improve subfamily splitting (instrument first)

**As v2.** Step 1 instrument both threading paths (`refine_split_families` is duplicated, ~10 LOC each). Step 2 analyze. Step 3 fix only if instrumentation reveals clear cause.

**R1 Eng caution**: do NOT lower `REFINE_BIMODALITY_THRESH` in isolation; paired check at refine.c:1056 skips valley test when bimodality ≥ 0.40. Lowering threshold inverts the safeguard.

**R1 Eng note**: Otsu uses `inst->divergence` which is substitution-only (no gaps). Indel-mode bimodality is invisible.

Test required: bimodal subfamily synthetic; both threading modes (per ENG-N6).

Effort: 1-2 days instrumentation; 1-3 more days fix.

### D' (rescoped per BIO-R1 + ENG-R2): Multi-pass with full-genome seeding

**v2 design**: Pass 2 on full genome + BLAST screen at merge. **R2 Eng new concern (Q4 audit + new edge case)**:

- **Q4 confirmed**: existing multipass.sh uses RESIDUAL masking (the OLD v1 design that bio rejected). HANDOFF's "Pass 2 = 0 families" result is from this WRONG design under the hash bottleneck. D' as specified is a genuine rewrite, not a script tweak.
- **R2 representative-direction concern**: when Pass 2 candidates merge with Library 1 families, union-find rep selection at refine.c:886-892 picks the family with most instances. If a Pass 2 candidate happens to have more instances than its Library 1 match, Library 1 family is absorbed → loses its already-refined consensus.

**v3 spec addition**: explicitly mark Library 1 families as "protected" (preferred representative) during the D' merge call. OR enforce the absorption direction by family-id ordering (Library 1 IDs < Pass 2 IDs).

**Pre-condition**: ENG-N1 controlled diagnostic must show seed-competition is real bottleneck. Otherwise no payoff.

Effort: 3 days (was 2-3; +0.5 for representative-direction guard, +0.5 for the wholly-new wiring).

### F' (rescoped per BIO-R1 + ENG-R2): BLAST-based instance recruitment for short elements

**v2 design**: replace `align_collect_instances` only for consensus_length < 500 with dc-megablast.

**R2 Eng new concerns**:
- (1) Minus-strand BLAST hits: `Instance.cons_start`/`cons_end` must be in CONSENSUS coordinates even for minus-strand. dc-megablast returns `strand=-1` for RC hits. Test required to assert correct mapping.
- (2) PADLENGTH offset: `Instance.position` is padded coordinate. BLAST returns genomic. Must add PADLENGTH or silent off-by-N bug.
- **R2 Bio note**: for MITEs, BLAST HSP chaining can produce alignments spanning both TIRs as one hit, misrepresenting boundaries. Future concern; not blocking.

Effort: 3-4 days (unchanged from v2 — concerns are spec clarifications, not new work).

---

## Tier 4: Defer (high uncertainty / explicit reject)

### E'': Structural-feature seed selection
Defer indefinitely. **R2 Bio note**: 15 kb LTR-pair detection window is right for Arabidopsis but should be **configurable** (25-30 kb for maize/wheat) when implemented.

### BIO-N3: TIR-seeding investigation for MITEs
Conditional on A revealing MITEs as major missed class.

### G: Loosen merge length-ratio guard
**REJECT**: do not loosen below 0.70. Confirmed by both reviewers (R1 + R2).

### H: Loosen Stage B fallback to 2-instance
**REJECT**: out of design scope, FP risk in pericentric region. Confirmed by both reviewers.

### I: Loosen max-divergence
**REJECT (blanket)**: HANDOFF §4 + §7 forbid. Conditional version (only families that already passed standalone fallback) acceptable per Bio-R1 but defer until ENG-N1 provides new evidence.

---

## Tier 5: Library-utility correctness

### BIO-N5: Verify canonical strand orientation for retrotransposon consensuses

For ≥5 kb library consensuses, BLAST against UniProt for RT/RH/IN domains. If domain order is INT-RH-RT (5'→3'), flip consensus to canonical. Add as post-MDL polishing step.

Doesn't affect family-level recall but materially affects library utility for downstream TE age estimation, K-divergence dating, CpG decay analysis.

Effort: 2-3 days (UniProt domain DB integration).

---

## NEW Tier 1.7: Cross-genome calibration documentation (per BIO-N9)

mdl-repeat must work on all genomes. v3 documents calibration adjustments needed for non-Arabidopsis:

- **MINTHRESH**: 3 may be too low for very large genomes (maize 2.3 Gb, wheat 17 Gb) where sequence-shared spurious families bloat library. Future investigation: MINTHRESH should scale with genome size (e.g., `max(3, log10(N) - 5)`).
- **ALIGN_MAX_EXTENSION**: bump unconditionally (per J' amendment) — maize/wheat have 10-19 kb routine retrotransposons.
- **coalesce-factor default**: 20× too aggressive for nested-stack genomes. Document recommended 5× via the new -coalesce-factor flag (ENG-N7).
- **L_pair detection window** (when E'' is implemented): configurable 15-30 kb based on genome.

These are advisory for future work; not current blockers.

---

## Recommended execution order (revised v3)

**Phase 1 — Classification, diagnosis, robustness (2-3 days, all parallel-able)**:
1. **A + BIO-N1 + BIO-N4 + BIO-N6 + Q2** (1 day): classify 24 missed families, sanity-check truth, library max length, ATHILA scope
2. **ENG-N1** (3 hours): controlled small-genome diagnostic
3. **ENG-N4** (1 hour): noise floor measurement
4. **ENG-N9.aux** (0.5 day): grep + audit all positional-comparison sites in refine.c
5. **ENG-N3** (0.5 day): dynamic HASH_SIZE
6. **ENG-N7** (1 hour): -coalesce-factor CLI

**Phase 2 — Multi-chr + scalability bug fixes (2 days, can parallel)**:
7. **B + Q6** (1 day): fix nested_containment cross-chr + tighten gate predicate; both threading paths
8. **ENG-N8 + ENG-N9** (0.5 day): coalesce + assemble cross-chr fix; bundled patch
9. **ENG-N2 + ENG-N10** (1 day): refine_prune memory + qsort

**Phase 3 — Biology-validated regression tests (1 day)**:
10. **BIO-N2** (1 day): synthetic ATHILA solo-LTR test, post-MDL assertion (per BIO-N7 + N10)

**Phase 4 — Conditional improvements (gated on Phase 1 results)**:
11. **J'** (0.5 day): bump ALIGN_MAX_EXTENSION unconditionally (per BIO-N9)
12. **C** (1-5 days): subfamily splitting instrument + fix
13. **BIO-N3** (1.5 days): MITE TIR-seeding fix if A reveals MITEs
14. **D'** (3 days): multi-pass with library-screen + rep-direction guard

**Phase 5 — Polishing (3-4 days)**:
15. **F'** (3-4 days): BLAST short-element recruitment + minus-strand + PADLENGTH spec
16. **BIO-N5** (2-3 days): canonical strand orientation polish

**Skipped (R1+R2 confirmed reject)**:
- E (deferred indefinitely; if ever, use E'' structural design)
- G (don't loosen below 0.70)
- H (don't implement)
- I (blanket reject; conditional only with new evidence)

---

## Open questions remaining after R2

1. (Unchanged) **What TE class are the 24 missed families?** — gates entire roadmap, A+BIO-N1
2. (Unchanged) **Does any chr4 library consensus reach 10 kb?** — Q2, single awk
3. (Unchanged) **Are ATHILA intervals in recall denominator or excluded?** — BIO-N6
4. (Unchanged) **Run-to-run variance at -threads 4?** — ENG-N4, before any recall-delta claim
5. (NEW R2) **In D', when Library 1 + Pass 2 merged, which side becomes union-find representative?** — needs explicit guard (v3 spec amends this)
6. (NEW R2) **What is the solo:full-length ratio for the 24 missed families if any are solo LTRs?** — drives BIO-N2 final test design (BIO-N7)
7. (NEW R2) **Are there other cross-chromosome unsafe sites in refine.c beyond N8/N9/Q6?** — ENG-N9.aux audit answers

---

## Test coverage requirements summary (R1 + R2)

| Item | New test required |
|---|---|
| A | Classification report; data analysis only |
| B+Q6 | Multi-chromosome synthetic + nested element synthetic; both threading modes |
| BIO-N2 | ATHILA 20+80 with internal-l-mer seeding; post-MDL assertion; both threading modes |
| C | Bimodal subfamily synthetic; both threading modes |
| D' | Pass 2 recovery test + rep-direction guard test (Library 1 family preserved) |
| F' | Short-element BLAST: 0-based BED, PADLENGTH offset, minus-strand cons_start/cons_end |
| ENG-N2 | malloc cap regression test (large genome simulation) |
| ENG-N3 | full 7/7 + chr4 smoke after dynamic HASH_SIZE |
| ENG-N7 | -coalesce-factor 0 vs default behavior diff |
| ENG-N8 | Multi-chromosome synthetic; assert no BED interval spans chr boundary |
| ENG-N9 | Multi-chromosome synthetic; assert no spurious cross-chr fragment assembly |
| ENG-N9.aux | Document audit findings for any other unsafe sites |
| ENG-N10 | n_accepted scaling regression (ensure prune doesn't hang on 10k+ families) |
