## Round 2 Engineering Review

_Reviewer: independent bioinformatics code review_
_Files read: REVIEW_ENG_R1.md, QUALITY_PROPOSAL_v1.md, QUALITY_PROPOSAL_v2.md,
src/refine.c (all 2554 lines), src/main.c (all 1206 lines), src/discover.c (head),
src/discover_internal.h, src/refine.h, src/candidates.h, /tmp/ath_bench/multipass.sh_

---

### Integration check

- **ENG-N1 in v2: correctly integrated.** Appears verbatim as its own numbered item (§ENG-N1, Tier 1). The controlled diagnostic is made independent of item A ("independent of A but cheap") and is pre-conditioned for Tier 3 items D' and I. Priority is correct.

- **ENG-N2 in v2: correctly integrated.** Appears as §ENG-N2 in Tier 1.5. The fix direction (sorted interval sweep-line, O(num_instances)) is correct. The test requirement (malloc cap regression) is carried over to the test-coverage table. One gap: see NEW item ENG-N8 below — ENG-N2 addresses `refine_prune_families` at refine.c:1769 but two other functions (`refine_coalesce_tandem_instances` and `refine_assemble_fragments`) contain the same O(genome_len) position iteration pattern using raw `position` without `seq_index` guards — though those are positional comparisons rather than allocations, they share the multi-chromosome correctness class of bugs. The v2 text scopes ENG-N2 narrowly to the `calloc`; the broader class is not addressed.

- **ENG-N3 in v2: correctly integrated** as §ENG-N3, Tier 1.5. The audit scope note ("audit every `for (h = 0; h < HASH_SIZE; h++)` loop in discover.c and discover_mask.c") is accurate and sufficient. The one addition worth flagging: `SMALLHASH_SIZE = 5003` (discover_internal.h:24) and `SMALLL = 6` (discover_internal.h:25) are separate constants used for mask consensus lookup; ENG-N3 scoping to HASH_SIZE only is correct since SMALLHASH_SIZE does not scale with genome length.

- **ENG-N4 in v2: correctly integrated** as §ENG-N4, Tier 1.5. The 1-hour effort estimate is accurate. The description of the non-determinism source (`__atomic_fetch_add` at refine.c:454 in merge_worker_fn) is correct.

- **ENG-N5 in v2: absent.** R1 ENG-N5 described a latent hazard: `SplitAnalysisArgs.families` captures `cl->families` before the parallel phase, and a subsequent `realloc` would invalidate it. R1 concluded this is safe as written (Phase 1 is complete before Phase 2a). v2 does not mention ENG-N5 at all, not even as "reviewed safe." Since R1's own conclusion was "safe as written, latent hazard," omitting it is acceptable, but if anyone restructures the split phases the hazard reactivates. Consider adding a one-line comment in refine.c:1566 noting the ordering invariant.

- **ENG-N6 in v2: partially integrated.** The test-coverage table correctly adds "Bimodal subfamily synthetic; both threading modes" for item C. The specific assertion (`n_splits = 1`) and the requirement to test both `-threads 1` and `-threads 4` paths match R1. However, the analogous gap R1 noted for item B (the existing testD may not assert solo-preservation) was already present in v2's B spec as "assert two distinct families remain" under both threading modes. Integration is correct.

- **ENG-N7 in v2: correctly integrated** as §ENG-N7, Tier 1.5. The 1-hour estimate is accurate; the default value (20.0) and disabled value (0) are specified. The test requirement is present in the test-coverage table. Confirmed: `refine_coalesce_tandem_instances` is called at main.c:1170 with the literal `20.0f` and `verbose` — no CLI hook exists. The ENG-N7 fix is unambiguous and low-risk.

---

### Open question audits (Q4, Q6)

**Q4 (multipass.sh): RESIDUAL-MASKED SEEDING — not full-genome with library screen.**

Read /tmp/ath_bench/multipass.sh in full. The actual implementation:

1. Pass 1: standard `mdl-repeat` on full genome → `library_pass1.fa`
2. RepeatMasker remaps Pass 1 library on the original genome → `pass1_rm.bed`
3. A Python snippet reads `pass1_rm.bed` and writes `N` over all RM-hit positions → `residual.fa`
4. Pass 2: `mdl-repeat` on `residual.fa` (N-masked genome) → `library_pass2.fa`
5. Combine: simple concatenation (`cat pass1.fa pass2.fa`) with renumbered IDs

This is the OPPOSITE of what D' now proposes. The existing script does residual-masked seeding (the v1 D design), which Bio-R1 and Eng-R1 both warned against: l-mer frequency artifacts at masked-unmasked boundaries, and the HANDOFF §4 "Pass 2 produced 0 families" failure mode. The script does NOT implement the D' revised design (full-genome seeding with growing-library screen).

Consequence for v2: D' is correctly specified as a rewrite relative to the current script. The controlled diagnostic (ENG-N1) is essential before committing to D' because the existing script's failure ("Pass 2 produced 0 families") may have been caused either by the masking artifact or by the hash-size bottleneck — it is genuinely ambiguous which one. The v2 text acknowledges this ambiguity ("still need the controlled diagnostic") — that note is accurate.

**Q6 (nested_containment_fraction multi-chr): BUG CONFIRMED, AFFECTS ITEM B.**

Code at refine.c:219-243 (full function):

```c
static float nested_containment_fraction(const CandidateFamily *shorter,
                                         const CandidateFamily *longer)
{
    ...
    for (int i = 0; i < shorter->num_instances; i++) {
        gpos_t s_start = shorter->instances[i].position;   // raw padded genome pos
        gpos_t s_end   = s_start + shorter->instances[i].aligned_length;
        ...
        for (int j = 0; j < longer->num_instances; j++) {
            gpos_t l_start = longer->instances[j].position; // raw padded genome pos
            ...
            gpos_t ov_e = (s_end < l_end) ? s_end : l_end;
            if (ov_e <= ov_s) continue;
            ...
        }
    }
```

There is no `seq_index` check anywhere in this function. `position` is the padded genome-wide coordinate (PADLENGTH + offset into the concatenated sequence buffer). In a multi-chromosome FASTA the sequences are concatenated with padding between them. For a genome with chr1 (100 Mb) + chr4 (20 Mb), a position value of 100,001,000 could refer to a chr1 instance of the longer family and a chr4 instance at a similar raw offset. The overlap check `ov_e <= ov_s` would correctly return no-overlap in most cases because the raw coordinate gap between chromosomes is the full chromosome length — but this relies on an implicit assumption that chromosomes are never laid out so that a chr4 position value falls within the range of a chr1 position + instance_length. In practice this is generally safe for large chromosomes, but it is a latent correctness issue for short contigs or scaffolds where the concatenated layout could bring coordinates from different sequences close together.

More concretely: the containment fraction computes whether shorter's instances are inside longer's instances. If shorter family has instances on chr1 and chr4, and longer family has instances only on chr1, the chr4 instances of shorter will correctly find no overlap (gap too large). But if a scaffold is 300 bp and the padding is 200 bp, a shorter instance at the end of scaffold A and a longer instance at the start of scaffold B could spuriously "overlap" in the raw coordinate space.

**Verdict**: Bug is present but only triggers on short-contig assemblies. For chr4 (18.6 Mb) and TAIR10 (5 chromosomes each >10 Mb), the practical risk is low. However, the proposal correctly requires fixing this before B's gate tightening, since the fix for B (adding `seq_index` check) is the standard correction. The v2 spec for B correctly lists "First, audit `nested_containment_fraction` for cross-chromosome bug" as step 1.

---

### NEW action items (Round 2 only)

**ENG-N8: `refine_coalesce_tandem_instances` does not check `seq_index` — can coalesce instances across chromosome boundaries.**

refine.c:2487-2496:
```c
for (int k = 1; k < f->num_instances; k++) {
    int    cur_idx = order[k].idx;
    Instance *cur  = &f->instances[cur_idx];
    gpos_t active_end = active->position + active->aligned_length;
    int    gap        = (int)(cur->position - active_end);
    if (cur->strand == active->strand &&
        gap >= -10 &&
        gap <= gap_threshold) {
```

Instances are sorted by raw `.position` (refine.c:2474-2477), which is the padded genome-wide coordinate. There is no `seq_index` check before the gap test. When two chromosomes are adjacent in the concatenated buffer and the boundary padding (PADLENGTH N-bases) happens to be smaller than `gap_threshold` (which is `20.0 * consensus_length`, potentially 20,000 × 20,000 = 400,000 bp for a long element), an instance at the end of chromosome A and an instance at the start of chromosome B would be coalesced into one BED interval spanning the chromosome boundary. This interval's coordinates would then be written to the BED file as if it were a single instance, with a `local_start`/`local_end` that is computed by subtracting the wrong `chr_offset` (whichever chromosome the position-lookup resolves to). The resulting BED line would have coordinates that exceed the chromosome length.

Fix: add `if (cur->seq_index != active->seq_index) { active_idx = cur_idx; active = cur; continue; }` before the gap test (refine.c:2494).

This is a correctness bug that can produce invalid BED coordinates in multi-chromosome inputs whenever `PADLENGTH < coalesce_factor * consensus_length`. For TAIR10 (5 chromosomes) with a 20 kb consensus and the default `coalesce_factor=20`, `gap_threshold` is 400 kb, which is far larger than PADLENGTH (which is typically a few hundred N-bases as a boundary delimiter). In a single-chromosome input (e.g., chr4 test) this cannot fire.

Severity: produces wrong BED coordinates silently. Classified 🔴 for multi-chromosome inputs.

**ENG-N9: `refine_assemble_fragments` sweep-line does not filter by `seq_index` — can detect spurious cross-chromosome co-occurrences.**

refine.c:2000-2007:
```c
for (int i = 0; i < entries_count; i++) {
    int fa = entries[i].family_idx;
    for (int j = i + 1; j < entries_count; j++) {
        if (entries[j].start - entries[i].start > D) break;
        int fb = entries[j].family_idx;
        if (fa == fb) continue;
```

`entries` is sorted by raw `.position` (refine.c:2956). When two chromosomes are concatenated, an instance at the end of chromosome A and an instance at the start of chromosome B are adjacent in the sorted `entries` array. If their raw position difference is ≤ D (up to 30,000 bp), they will be counted as co-occurring and their pair's `count` and `same_dir_count` will increment. If the pair's co-occurrence count meets the assembly threshold, the two families will be assembled as if they were fragments of the same TE — when they are actually on different chromosomes and unrelated.

This is a less critical bug than ENG-N8 because the assembly gate also checks MDL improvement (the assembled family must have better MDL than the sum of its parts), which provides some protection: a spurious cross-chromosome assembly will generally not have genome-wide instances supporting the combined consensus and will fail MDL. But it is not impossible for two families that happen to have real instances near the chromosome boundary to be misassembled. Fix: add `if (entries[j].seq_index != entries[i].seq_index) break;` (breaking is only safe if entries are sorted by `(seq_index, position)` rather than raw position — alternatively `continue` instead of `break`).

Note: the `InstanceEntry` struct at refine.c:1864-1870 does NOT have a `seq_index` field; it would need to be added:
```c
typedef struct {
    gpos_t start; gpos_t end; int family_idx; int instance_idx; int8_t strand;
} InstanceEntry;
```
To fix this, `seq_index` must be added to `InstanceEntry` and populated at refine.c:1946-1950. The sort comparator must be updated to sort by `(seq_index, start)` to make the break-on-D logic correct.

Severity: silently produces wrong family assemblies on multi-chromosome inputs. Classified 🟡: the MDL gate provides partial protection.

**ENG-N10: `refine_prune_families` sort is O(n^2) bubble sort — acceptable now but breaks at scale.**

refine.c:1757-1766:
```c
for (int i = 0; i < n_accepted - 1; i++) {
    for (int j = i + 1; j < n_accepted; j++) {
        if (order[j].score < order[i].score) {
            IdxScore tmp = order[i];
            order[i] = order[j];
            order[j] = tmp;
        }
    }
}
```

This is a selection-sort (O(n^2)) over `n_accepted` families. On chr4 with ~100 accepted families this is negligible. On a full wheat genome, if 10,000+ families are accepted, this is 10^8 operations in the prune sort alone. This is not a correctness issue, but combined with the O(genome_len) calloc from ENG-N2, `refine_prune_families` becomes the double bottleneck on large genomes. Fix is one line: replace with `qsort`. Classified 🔵 (nit, not bio-correctness), but noted because it pairs with ENG-N2.

---

### Verification of rescoped items

**D' (multi-pass full-genome seeding with library screen): accept with reservation.**

The rescoped design (full genome re-seeding, exclude candidates ≥80-80-80 contained in Library 1 at the merge stage, no mask-export API) correctly addresses R1's mask-state concern. The "no mask-export API needed" claim is accurate — D' as described does not require exporting `C->removed[]`.

One new technical concern: the "exclude candidates ≥80-80-80 contained in Library 1" step requires running the semi-global alignment of EVERY Pass 2 candidate against the ENTIRE Library 1. If Library 1 has 200 families and Pass 2 discovers 500 candidates, this is 100,000 pairwise alignments before the Jaccard pre-filter. The Jaccard pre-filter (REFINE_MIN_JACCARD = 0.15) would eliminate most pairs quickly, so this is probably fast in practice, but the proposal does not address this computational path. Worth confirming before implementation that the existing `refine_merge_families` call handles the combined CandidateList correctly — specifically, if Library 1 families are marked with `mdl_score > 0` before the merge call, the merge stage should not absorb them into Pass 2 candidates (it would absorb the Pass 2 candidate into Library 1 instead, since Library 1 has more instances). This directionality is correct but depends on which side becomes the "representative" in the union-find (refine.c:886-892: representative = family with most instances). If Library 1 families have more instances than Pass 2 candidates, absorption direction is correct. If a Pass 2 candidate somehow has more instances than a Library 1 family, the Library 1 family would be absorbed — which is the wrong direction. This edge case needs explicit protection.

**E'' (structural feature seeding): accept deferral.**

No further concern. The "defer to Tier 4" decision is correct given discover.c's fragility.

**F' (BLAST for short elements): accept scope, flag one unresolved concern.**

The `cons_start`/`cons_end` mapping concern from R1 is acknowledged but not resolved in v2. It is listed as a reason for the effort increase (3-4 days vs 2-3), which is the right framing. However, the proposal does not specify HOW the mapping will be done: BLAST tabular output gives query start/end (into consensus) and subject start/end (into genome), which maps directly to `cons_start`/`cons_end` and `position`/`aligned_length`. The concern is whether the BLAST query orientation (always 5'→3' on the + strand query) matches the convention in `Instance.strand`. If dc-megablast returns a reverse-complement hit, the strand field and the `cons_start`/`cons_end` interpretation must both be set correctly. The proposal should add an explicit test for a minus-strand BLAST hit to confirm the `strand` field is set correctly and `cons_start`/`cons_end` are set in consensus coordinates (not genomic coordinates).

The `PADLENGTH` offset concern from R1 (BED coordinates off by PADLENGTH bases) is also not explicitly resolved in v2. The proposal should state where in the integration PADLENGTH is added/subtracted. Confirmed from main.c:104: `gpos_t raw_start = inst->position - PADLENGTH` — so `inst->position` is stored as padded, and PADLENGTH is subtracted at output time. Any BLAST integration that sets `inst->position` must add PADLENGTH to the genomic coordinate.

**J' (ALIGN_MAX_EXTENSION conditional): accept.**

The conditional structure (check Q2 first, bump only if any consensus reaches 10 kb) is correct. The R1 concern about the prior `-L` sweep being on a now-deprecated metric is correctly noted. The 0.5-day effort estimate is accurate.

**Tier 1.5 (ENG-N2/N3/N7): priority correctly elevated.**

ENG-N2 and ENG-N3 are correctly classified as blocking for the "universal tool" stated goal. ENG-N7 is cosmetic but cheap. The ordering (Phase 3, parallel with Phase 2) is correct — these do not depend on diagnostic results.

---

### Effort re-estimates if v2's are wrong

- **ENG-N8 (new): 0.5 day.** Requires adding `seq_index` field to `InstanceEntry` struct, populating it at the fill loop, updating the sort comparator to sort by `(seq_index, start)`, and adding the cross-chromosome guard in the gap test. Ripple to multi-chr test.
- **ENG-N9 (new): 0.5 day** if ENG-N8 is done first (shares the `seq_index` addition to `InstanceEntry`). Could be folded into the same patch. If done separately: 1 day.
- **D' pre-condition on Library 1 representative protection: 0 extra days** — this is a design decision that can be made at implementation time, but the author should document it explicitly before coding to avoid re-discovering it as a bug post-integration.

---

### Open questions remaining

1. (Unchanged from R1 Q1) At which pipeline stage do the 24 missed chr4 families drop out? The `-trace-dir` dumps can answer this without code change. Entire priority ordering of B through E depends on this.

2. (Unchanged from R1 Q2) Is the 0-split observation on chr4 explained by the `REFINE_MIN_SPLIT_INSTANCES = 10` floor (refine.h:51) or by genuine unimodality? The C instrumentation resolves this.

3. (Unchanged from R1 Q3 / now verified as Bug confirmed) `nested_containment_fraction` multi-chr: bug confirmed at refine.c:219-243. Fix is step 1 of item B — adding `if (shorter->instances[i].seq_index != longer->instances[j].seq_index) continue;` inside the inner loop. Safe to close this open question with "bug present, fix specified."

5. (Unchanged from R1 Q5) Run-to-run variance with `-threads 4`? Sets noise floor for evaluating any fix. ENG-N4 resolves this — has it been run yet?

6. (Now answered) Q4 is answered: multipass.sh uses residual-masked seeding (old D design), NOT the new D' design. This confirms D' is a genuine rewrite.

7. (NEW) Does the `refine_assemble_fragments` sweep-line sort by raw `position` or by `(seq_index, position)`? Confirmed: sort is by raw `position` only (refine.c:2956 via `cmp_instance_entry` which compares `.start` only). This is the source of ENG-N9.

8. (NEW) In D', when the combined CandidateList (Library 1 + Pass 2) is passed to `refine_merge_families`, which family becomes the union-find representative when a Pass 2 candidate matches a Library 1 family? The answer matters for correctness of the "exclude via merge" strategy. Library 1 families typically have more instances (full-genome refine) than Pass 2 candidates (discovered from scratch on the same full genome with a screen). This should hold in practice but needs a guard.

