## Round 1 Engineering Review

_Reviewer perspective: C/algorithm engineering, scalability, correctness risk_
_Files read: QUALITY_PROPOSAL_v1.md, HANDOFF.md, src/main.c, src/refine.c, src/discover.c, src/discover_internal.h, src/align.c, src/mdl.c, src/refine.h, src/align.h_

---

### Confirmed (items I agree are well-scoped)

- **A (diagnostic first)**: Extracting the 24 missed families and characterizing them is the right first move. The tool produces trace dumps via `-trace-dir` already (main.c:51-168), so stage-by-stage FASTA can be diffed against truth without any code change. Only risk is interpretation time; no regression possible from a read-only experiment.

- **B (nested containment gate tightening)**: The gate is isolated in `refine.c:219-243` (`nested_containment_fraction`). It is triggered by two guards at refine.c:557-561 (parallel path) and refine.c:748-751 (sequential path). Both paths are structurally identical and call the same helper, so a fix only needs to touch `nested_containment_fraction` and possibly the trigger predicate. Low regression risk because the gate only fires when `longer >= 3 × shorter AND n_short_instances >= 3`; the trigger condition itself is not frequently satisfied on chr4 (the HANDOFF says testD was the main beneficiary). The proposal's fix direction (requiring non-prefix/suffix relationship and evidence of solo instances) is feasible without touching `discover.c` or mask state.

- **C (instrument split before relaxing thresholds)**: Step 1 of C (instrumentation) is strongly endorsed. `refine_split_families` currently emits a log line only when a split is *accepted*; adding a `-v` line at each rejection point (bimodality too low, valley not deep, n_lo/n_hi too small, MDL gain negative) is ~10 LOC and zero risk. Do this before touching any threshold. The sequential and parallel paths in `refine_split_families` are structurally duplicated (refine.c:1374-1551 vs 1553-1727); any instrumentation change must be applied to both copies.

---

### Concerns / risks

**Item B — nested gate trigger asymmetry**

The trigger requires `longer >= 3 × shorter` (refine.c:557), meaning a 500 bp element inside a 1497 bp element fires the gate, but a 500 bp element inside a 1499 bp element does NOT. This 3× hard cutoff is arbitrary and the proposal to tighten the gate will inherit the same asymmetry. The bigger concern: both the parallel merge worker (`merge_worker_fn`) and the sequential path have nearly identical logic that must be kept in sync. Any change to the gate predicate in one must be duplicated in the other, or a shared helper introduced. A divergence already happened once (the parallel worker lacks the relaxed-identity branch's nested check that the sequential path has—compare refine.c:740-762 sequential vs the parallel path around refine.c:539-573).

**Item C — Otsu on discover-stage divergence values has a circularity**

The divergence values fed to `otsu_threshold` come from `inst->divergence`, which is set by `align_banded` (align.c:444). The banded DP scores substitutions only, not gaps (align.c:427 note). This means the divergence distribution is systematically narrower than the true edit-distance distribution. A bimodal family where both modes differ primarily by indel rate (as opposed to substitution rate) will appear unimodal here and Otsu will never fire. This is a known limitation already documented (align.c:424-431) but the proposal does not account for it when claiming a threshold relaxation will fix the 0-split-on-chr4 observation.

**Item C — threshold relaxation before instrumentation inverts the scientific order**

The proposal lists instrumentation as step 1 and relaxation as step 2, which is correct. The risk is that pressure to see improvement leads to relaxing REFINE_BIMODALITY_THRESH (refine.h:54) before the logs are understood. REFINE_BIMODALITY_THRESH=0.20 is a hard-coded `#define`. If it is changed, both the parallel and sequential paths in `refine_split_families` use it, but there is also a second threshold at refine.c:1056: the valley depth check is skipped when bimodality >= 0.40. Lowering REFINE_BIMODALITY_THRESH to 0.10 without revisiting the 0.40 skip bound would mean families with 0.10-0.20 bimodality bypass the valley check, admitting splits that may not be real subfamilies. Flag this interaction explicitly before any tuning.

**Item D — multi-pass depends on the engineering bottleneck (HANDOFF §6)**

The HANDOFF records that Pass 2 "produced 0 families" when it eventually completed via the failure handler. The proposal attributes this to the timeout, but this is not established. If Pass 2 produced 0 families because the residual genome genuinely contains no new seeds above MINTHRESH (all surviving high-freq l-mers have been masked), then fixing HASH_SIZE will not change the outcome. Multi-pass wiring in main.c would require: (1) a mask-export API from `discover_families` (currently mask state lives inside DiscoverContext on the heap and is freed at function exit — discover.c has no "export current mask" interface), (2) applying that mask to a second discover call on the same genome, (3) deduplicating the combined library via the merge stage. None of these interfaces exist. The effort estimate of 1 day is too low; 2-3 days is more realistic for a safe implementation that does not touch the mask internals beyond what already exists.

The HANDOFF §4 explicitly records the multi-pass result as "inconclusive — not provably stuck but unacceptable runtime; when it eventually ran to completion Pass 2 produced 0 families." Before building infrastructure, run HANDOFF §5's controlled diagnostic (option B in §9): a 1 Mb background + 10 copies of a known-missed cluster representative. If that small genome also misses the family, the bottleneck is seed-competition, not hash collisions, and multi-pass will not help regardless of HASH_SIZE.

**Item D — mask state is fragile (HANDOFF warning)**

HANDOFF §4 documents a 0.027 F1 drop from mask-adjacent changes (P0v1 attempt). Any multi-pass that touches `C->removed[]` — even by reading it — risks introducing a subtle interaction. The HANDOFF specifically warns "don't touch these params unless you have evidence." Multi-pass requires at minimum reading the mask state post-Pass-1, which is currently destroyed when `discover_families` returns and frees the DiscoverContext.

**Item E — repeat-discriminative seed scoring is high-risk, high-effort**

The "bootstrap from Pass 1 regions" approach changes `build_headptr_internal` or `find_besttmp`, both of which are in the fragile core of discover.c. The HANDOFF documents 12 failed optimization attempts in this file area. Minimizer/syncmer-based seeding would require a parallel implementation of the l-mer counting phase that doesn't break `build_all_pos`'s reliance on the `struct posllist` structure built from the headptr. Effort of 4-7 days is plausible for a prototype, but regression risk is high. This should be last.

**Item F — BLAST-based instance recruitment changes the output format contract**

`align_collect_instances` (align.c:464) currently produces `Instance` structs with fields `position`, `aligned_length`, `cons_start`, `cons_end`, `num_edits`, `divergence`, `score`, `strand`. BLAST's HSP output would need to populate all of these, particularly `cons_start`/`cons_end` which are used in `build_subset_consensus` (refine.c:1096-1133) to index back into the consensus during split. A BLAST wrapper that maps incorrectly to `cons_start`/`cons_end` would silently produce wrong split consensi. The BED coordinate output in `main.c:104-123` (the trace dump) and `output.c` both derive final coordinates from `inst->position - PADLENGTH`, so any BLAST integration must correctly handle PADLENGTH offset or the BED will be wrong by `PADLENGTH` bases.

**Item G — loosening length-ratio guard from 0.7 to 0.5**

The 0.7 guard was introduced as Stage B fix 2 specifically because merging vastly-different-length families destroyed 2.36 Mb of coverage (cited in refine.c:484). Loosening to 0.5 re-opens that risk. This is safe to test only after verifying that the merge MDL gate (estimate_merge_score) provides sufficient protection at lower length ratios, which depends on how accurate the centroid-distance approximation in `estimate_merge_score` is for 50%-length pairs. The approximation uses `half_extra = (1 - identity) / 2` (refine.c:397), which is only accurate near 80% identity; at 0.5 length ratio the shorter family's instances will be clipped much more aggressively, making the approximation unreliable. The MDL gate may underestimate the damage.

**Item H — loosening standalone fallback from 3 to 2 instances**

The threshold is at mdl.c:272 (`STANDALONE_MIN_INST = 3`). Lowering to 2 will admit families that currently fail the `num_instances < 2` pre-screen at mdl.c:281 — those are different. Families that pass the 2-instance pre-screen but have `num_instances == 2` are currently rejected if they also fail the standalone fallback because 2 < 3. Admitting them means accepting 2-copy elements, which the design philosophy explicitly marks as out of scope (HANDOFF §0: "2-copy — out of scope"). This change conflicts with the stated design boundary and should not be done without a design-philosophy review.

**Item I — loosening max-divergence from 0.30 to 0.40**

This parameter was swept and found unchanged at 0.40/0.50 (HANDOFF §4 and §7). The result is documented: "F1 unchanged." Do not redo this without new evidence from the controlled diagnostic experiment. It is in the HANDOFF's explicit "don't touch without evidence" list.

**Item J — increasing ALIGN_MAX_EXTENSION from 10000 to 20000**

Also swept: "-L 30000 → F1 unchanged (max_cons same as -L 10000)" and "-L 50000 → F1 unchanged" (HANDOFF §7). This is a dead end. The ALIGN_MAX_EXTENSION constant is at align.c:21 and is distinct from `-L` (the discovery-stage extension parameter), but the HANDOFF's -L sweeps were testing the same bottleneck (whether extension length is the limiting factor). The controlled diagnostic experiment is needed to determine if extension distance ever reaches the cap on known-missed families before repeating this test.

---

### NEW action items not in proposal

**ENG-N1: Run HANDOFF §5 controlled diagnostic before any code change**

The proposal's entire Tier 2 and Tier 3 rationale rests on an unverified assumption about the bottleneck mechanism. Build a 2 Mb synthetic genome (1 Mb background + 10 copies of the longest missed-cluster representative from the chr4 analysis in item A). If mdl-repeat finds the family on this small genome, the problem is genome-scale specific (hash collisions, seed competition on large N). If it does not find it, every multi-pass / parameter-tuning proposal is addressing the wrong layer. This is a 30-minute experiment once item A identifies the target sequence. Do this before committing to B, C, or D.

**ENG-N2: Fix `refine_prune_families` memory allocation on large genomes**

`refine_prune_families` calls `calloc((size_t)genome_len, 1)` at refine.c:1769 to allocate a per-base coverage byte array. On a 17 Gb wheat genome this is 17 GB of heap in a single allocation, which will fail or thrash. The array is used only to count coverage up to 255; a sparse representation (sorted instance intervals + sweep line, or even a hash of occupied ranges) would reduce this to O(instances) instead of O(genome_length). For Arabidopsis (119 Mb) this is 119 MB which is acceptable; for maize (2.3 Gb) it is 2.3 GB which is at the edge; for wheat it is fatal. The proposal's stated goal of being a "UNIVERSAL tool" makes this a correctness-class issue for large genomes even if chr4 is unaffected.

**ENG-N3: Fix `HASH_SIZE` being compile-time fixed before any large-genome work**

`HASH_SIZE = 16000057` is a `#define` in `discover_internal.h:22`. For a 119 Mb genome with an l-mer of ~14 (ceil(1 + log4(119e6)) ≈ 14), the expected number of distinct l-mers is up to 4^14 = 268 million, which is 17× the hash table size. This forces average chain lengths of ~17 on a genome that is already repeat-rich, making counting O(N × average_chain_length) rather than O(N). This is the bottleneck HANDOFF §6 identifies. Any item D (multi-pass), E (new seeding), or large-genome benchmark is gated on this fix. The fix sketch in HANDOFF §9 option A (dynamic `HASH_SIZE = max(16M, 4 × N / l)`) is the right direction and has low risk since it only changes initialization of the headptr array. However, `struct llist **headptr` is allocated as `malloc(HASH_SIZE * sizeof(struct llist*))` in `discover_families`; making HASH_SIZE dynamic requires passing it as a parameter or computing it at context-init time and allocating accordingly. Estimated effort: 2-3 hours for the counting phase; `build_all_pos` and `mask_headptr` also iterate `for (h = 0; h < HASH_SIZE; h++)` and must use the same dynamic value.

**ENG-N4: Parallelism non-determinism in merge stage**

`merge_worker_fn` (refine.c:445) uses `__atomic_fetch_add` to distribute rows across threads. The collected pairs are applied sequentially (refine.c:848-864), but the order in which workers add to their `pairs[]` arrays depends on scheduling. Two runs with the same seed can produce different merge sets when thread ordering differs. This means the baseline recall number (0.389) may have run-to-run variance. Before attributing a 0.5 pp change to any code modification, confirm the variance by running the baseline twice. If variance exceeds 0.5 pp, changes of items G/H/I cannot be evaluated reliably.

**ENG-N5: `split_analysis_worker` reads `cl->families` while main thread has snapshot pointer; parallel split applies in Phase 2a sequentially, but fam pointer at refine.c:1650 is re-acquired via `fam = &cl->families[fi]` after realloc — this is correct. No race here, but the `SplitAnalysisArgs.families` pointer (refine.c:1167) captures `cl->families` before the parallel phase; if Phase 2a causes a `realloc`, the snapshot pointer becomes stale. The split analysis phase is complete before Phase 2a begins (refine.c:1594), so the snapshot pointer is only read during Phase 1. This is safe as written, but is a latent hazard if the two phases are ever interleaved.**

**ENG-N6: Add `-trace-dir` output counting to the regression tests**

`run_tests.sh` currently checks for family detection (PASS/FAIL). For item C (split threshold tuning), the tests should also verify that `n_splits` on the synthetic bimodal test case (if one existed) equals the expected count. No such test currently exists. Before relaxing REFINE_BIMODALITY_THRESH, add a synthetic test with a known two-subfamily family and assert it is split exactly once. Without this, threshold relaxation has no regression anchor.

**ENG-N7: `discover_chunked` assigns sequential IDs to combined families (main.c:715), but the chunk-level discover calls have their own internal ID numbering. After concatenation, IDs are reassigned sequentially. The merge stage's union-find uses array indices, not IDs, so this is safe. The BED and FASTA outputs use `f->id` for the label (main.c:81 trace, output.c). IDs will be dense 0..N-1 after chunking but may not correspond to the per-chunk discovery order. This is a cosmetic issue only (no biological correctness impact), but worth noting if the paper describes family IDs as meaningful.**

---

### Should-not-do

- **Item I (max-divergence 0.30 → 0.40)**: HANDOFF §4 and §7 explicitly document that this parameter was swept to 0.40 and 0.50 with "F1 unchanged." Repeating it without new evidence wastes time and risks introducing subtle regressions. Do not redo.

- **Item J (ALIGN_MAX_EXTENSION 10000 → 20000)**: The equivalent parameter (-L) was swept up to 50000 with zero recall change (HANDOFF §7). The extension is not the bottleneck. Do not redo.

- **Item E before controlled diagnostic**: The "minimizer/syncmer-based seeding" path requires rewriting the core of `build_headptr_internal` (discover.c:269-325), which is the most fragile part of the codebase (HANDOFF §4 records 12 failed attempts in this region). Starting E without first knowing whether the problem is hash-collision-related (ENG-N1/N3) is building the wrong thing at high cost.

---

### Test coverage requirements per item

- **A**: No code test needed; this is a data analysis. Document findings in a report file.
- **B**: Add a synthetic test case where a nested element (SINE-like, ~300 bp) appears both solo (5+ copies) and embedded in a longer element (LINE-like, ~900 bp, 5+ copies). Verify that the two families are NOT merged. Currently testD covers the embedded case but may not assert the solo-preservation behavior explicitly.
- **C**: Before any threshold change, add a synthetic bimodal test: one family with instances at two clearly-separated divergence levels (e.g., 0.05 and 0.20). Verify `n_splits = 1` and that both halves survive MDL. Use `-vv` output counting to confirm the split fires. Apply this test to both the sequential (`-threads 1`) and parallel (`-threads 4`) paths, as both are separately implemented.
- **D**: If multi-pass is implemented, add a test where Pass 1 dominates on one family (high copy) and Pass 2 recovers a different family (low copy, but with instances masked in Pass 1). Verify both appear in the combined output.
- **F**: If BLAST is integrated, add a test verifying that instance BED coordinates are 0-based and that `cons_start`/`cons_end` in the Instance struct are correctly populated from the BLAST HSP alignment.
- **G/H**: Parameter-only changes. Add a test that a 2-copy family (item H) is correctly rejected (not admitted) when `STANDALONE_MIN_INST` is 3, and correctly admitted when lowered to 2. Keep both as a flag-controlled knob test.
- **ENG-N2**: Add a test that `refine_prune_families` does not allocate more than 500 MB when called with a synthetic large `genome_len` value. This can be a compile-time assertion or a runtime check with a configurable cap.
- **ENG-N3**: After making HASH_SIZE dynamic, run the full synthetic suite (7/7) and the chr4 smoke test to confirm no regression.

---

### Effort estimate corrections

- **Item A**: 1-2 hours estimate is accurate for data characterization, assuming the truth BED and missed-family BED are already available (they are, per HANDOFF §2). Add 0.5-1 hour if the `-trace-dir` stage-by-stage dumps need to be cross-referenced against truth to determine at which pipeline stage families are lost.
- **Item B**: 0.5 day estimate is realistic for the gate tightening alone. Add 0.5 day for the new synthetic test (ENG-N6 for B). Total: 1 day.
- **Item C**: 1-2 day estimate assumes threshold tuning is straightforward after instrumentation. If instrumentation reveals the split is blocked at the MDL gate (split_score <= orig_score, refine.c:1303) rather than the bimodality threshold, the fix is different (scoring, not thresholds) and would take 2-3 days more.
- **Item D**: 1 day is too low. Exporting mask state from discover.c requires a new API (discover.c has no export interface for `C->removed`). Building a residual-genome pass, a deduplication step, and the wiring in main.c is 2-3 days minimum. Add 1 day for the controlled diagnostic to first verify the mechanism. Realistic: 3-4 days if the diagnostic shows multi-pass is worth pursuing.
- **Item E**: 4-7 day estimate might be optimistic given that 12 prior discover.c attempts all failed. Estimated realistic effort for a working implementation that passes 7/7 synthetic tests: 7-14 days with high uncertainty.
- **Item F**: 2-3 day estimate assumes BLAST is a simple subprocess call. The difficulty is mapping BLAST tabular output to the `Instance` struct (especially `cons_start`/`cons_end`) correctly. If BLAST's coordinate convention differs from mdl-repeat's internal PADLENGTH-offset convention, silent coordinate bugs will result. Estimate 3-4 days including testing.

---

### Open questions

1. At which pipeline stage do the 24 missed chr4 families drop out? The `-trace-dir` dumps can answer this without any code change. If they are present post-merge but rejected by MDL, the fix is MDL sensitivity. If they never appear post-discover, the fix is seeding. The entire priority ordering of B through E depends on this answer.

2. Is item C's 0-split observation on chr4 explained by genuinely unimodal families (all copies of each TE class are at similar divergence — expected for recently active TEs), or by the `REFINE_MIN_SPLIT_INSTANCES = 10` floor (refine.h:51) filtering out families before Otsu even runs? The instrumentation in C step 1 should log how many families are skipped by the 10-instance floor vs how many fail bimodality.

3. Does the `nested_containment_fraction` function correctly handle multi-chromosome genomes? Instances from different chromosomes could spuriously appear to "not overlap" with the longer family's instances, leading the gate to veto a legitimate merge of fragments from different chromosomes. The function (refine.c:219-243) uses raw `position` values without sequence index checks. If two families each have instances on chr1 and chr4 but the inner family's chr4 instances don't overlap chr4 instances of the outer family, the containment fraction could be falsely low. This is a correctness question for item B's fix.

4. For item D: does the current `multipass.sh` wrapper at `/tmp/ath_bench/multipass.sh` implement residual masking by excluding known-instance positions from the Pass 2 seed scan, or does it simply re-run discovery on the same genome and deduplicate? The HANDOFF description is ambiguous. This determines whether multi-pass as currently scripted is testing the right hypothesis.

5. Is the `coalesce_factor = 20.0` currently in main.c:1170 documented anywhere as a tunable parameter? It is not exposed via CLI and not mentioned in the usage string (main.c:752-793). If a user needs to disable tandem coalescing for a non-tandem-array-heavy genome (e.g., a compact fungal genome), there is no mechanism to do so. Consider exposing `-coalesce-factor` as a CLI option with 20.0 as default and 0 as "disabled."
