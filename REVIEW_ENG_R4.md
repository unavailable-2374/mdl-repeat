## Round 4 Engineering Review

_Reviewer: independent bioinformatics code review_
_Files read: REVIEW_ENG_R1.md, REVIEW_ENG_R2.md, REVIEW_ENG_R3.md,
QUALITY_PROPOSAL_v3.md, QUALITY_PROPOSAL_v4.md,
src/refine.c (all 2554 lines), src/mdl.c (all lines)_

---

### Integration check (R3 items in v4)

**ENG-N11 (check_instance_overlap fix spec): correctly integrated with one documentation error.**

The five-site table (lines 268, 275, 284, 307, 311) is correctly reproduced in v4. The fix direction is correct: guard at the inner forward-scan loop body (line 283, sorted path) and inside the linear-fallback inner loop (line 310/311), using `seq_index` from the already-populated `Instance` struct (candidates.h:22). The `seq_index` field is confirmed present in the struct and requires no struct change.

One documentation error carried from R3 into v4 verbatim: v4 line 93 states "The function is called at the merge gate (refine.c:557 + 748)." Those lines are where `nested_containment_fraction` is called (the B+Q6 bug), not `check_instance_overlap`. The actual call sites for `check_instance_overlap` are four, not two: lines 496, 532, 696, and 787. This is a citation error only; the fix spec (guards inside the function body) is unaffected and remains correct.

Minor precision: the table labels refine.c:307 as "linear-fallback outer loop — UNSAFE." Line 307 reads `a->instances[i].position` for family `a`'s outer loop; the comparison against `b` instances has not yet happened. The actual unsafe comparison is between the outer-loop `a_start` and the inner-loop `b_start` at line 311. Labeling 307 as "UNSAFE" is slightly misleading — the unsafe event is the combination of 307 (outer) and 311 (inner). The fix correctly targets line 310 (inner loop body), which is right. The table labeling is imprecise but does not cause a wrong fix.

Conclusion: ENG-N11 correctly integrated. One call-site citation is wrong (557+748 should be 496+532+696+787). Fix spec is correct.

**ENG-N8 clarification (runtime guard alone): correctly integrated.**

v4 explicitly resolves: "ship the runtime guard only (`if (cur->seq_index != active->seq_index)`) without changing the sort key." This is correct — after the guard fires and `active = cur; continue;`, subsequent iterations from the new chromosome test against the new `active` correctly. The sort by raw position is safe here because the guard intercepts cross-chromosome adjacency at runtime regardless of sort order. The clarification is unambiguous and the decision is sound.

**ENG-N9 clarification (sort change alone): correctly integrated.**

v4 explicitly resolves: "sorting `cmp_instance_entry` by `(seq_index, start)` alone is sufficient — the existing `break` at refine.c:2005 fires naturally at the chromosome boundary because `entries[j].start - entries[i].start` becomes the full chromosome offset (>> D)." Confirmed from the code: `cmp_instance_entry` (refine.c:1888-1895) currently sorts by `.start` only, and `InstanceEntry` (refine.c:1864-1870) has no `seq_index` field — both still require the changes ENG-N9 specifies. The clarification correctly eliminates the need for an in-body guard in the sweep-line loop body, simplifying implementation. Correctly integrated.

**ENG-N9.aux (10-site enumeration): correctly integrated.**

v4 closes the audit: 10 `.position` sites enumerated, 5 unsafe (all in `check_instance_overlap`, all captured by ENG-N11), 2 safe by design (refine_prune_families, cov[] array), 1 covered by ENG-N8 (coalesce sort key), 2 covered by B+Q6 (nested_containment_fraction). Confirmed against live `grep -n '\.position' src/refine.c` output — the 10 sites match exactly. Audit is complete and accurate.

---

### NEW action items (Round 4 only)

**ENG-N12: `mdl_select_library` allocates O(genome_len/8) bitmap — same scalability class as ENG-N2, not covered by ENG-N2's scope.**

mdl.c:247-248:
```c
size_t bitmap_bytes = ((size_t)genome_len + 7) / 8;
uint8_t *covered = calloc(bitmap_bytes, 1);
```

ENG-N2 covers `refine_prune_families` at refine.c:1769 (`calloc(genome_len, 1)` — 1 byte per base). `mdl_select_library` has an independent allocation of `genome_len/8` bytes (bitpacked). For wheat (17 Gb): 17e9/8 = 2.1 GB in a single `calloc`. This is 8x smaller than refine.c's allocation for the same genome, but it is still a fatal allocation on 17 Gb inputs (2.1 GB requested in one `calloc` call will fail on many HPC nodes with per-process memory limits, and the error path at mdl.c:249 returns an empty result that will silently discard all families). For 2.3 Gb maize it is 287 MB — at the edge of acceptability. v4's ENG-N2 fix (sorted instance sweep-line) applies identically here: walk instances sorted by (seq_index, position), track exclusive coverage by interval sweep without a per-base array.

The v4 test-coverage table lists "malloc cap regression test" for ENG-N2 but not for mdl_select_library. Any large-genome scalability test must cover both allocations.

Severity: same class as ENG-N2. For the stated goal of "universal tool" covering wheat and maize, this is a blocker. Note it is 8x less severe than ENG-N2 in memory consumption (bitmap vs byte array), so it triggers at larger genome sizes, but the fix is the same sweep-line replacement.

Fix direction: replace the genome_len/8 bitmap in `mdl_select_library` with the same sorted-interval sweep-line approach proposed for ENG-N2. Bundle both into a single "O(genome_len) coverage array elimination" patch. Effort: add 0.5 day to ENG-N2's existing 1-day estimate (total 1.5 days for both sites).

---

### Open questions remaining

1. (Unchanged from R1/R2/R3) At which pipeline stage do the 24 missed chr4 families drop out? Entire priority ordering of B through E depends on this.

2. (Unchanged) Is the 0-split observation on chr4 explained by the REFINE_MIN_SPLIT_INSTANCES=10 floor or genuine unimodality?

3. (Unchanged) Run-to-run variance with -threads 4? ENG-N4 must run before any recall-delta claim is trusted.

4. (Unchanged) In D', which side becomes union-find representative when Library 1 + Pass 2 are merged? v4 specifies the fix (id-based tiebreak); implementation must commit to one approach before coding.

5. ENG-N11 call-site citation: v4 cites "refine.c:557 + 748" as call sites for `check_instance_overlap` — these are wrong (they are `nested_containment_fraction` call sites). The actual call sites are 496, 532, 696, 787. This should be corrected in any future version to avoid confusion for the implementer.

---

### Final assessment

ENG-N11 is correctly integrated with one non-critical citation error (call-site line numbers). Both spec clarifications (ENG-N8 runtime-guard-only, ENG-N9 sort-change-only) are correctly integrated and technically sound. The ENG-N9.aux 10-site audit is confirmed complete and accurate against the live code.

One new item (ENG-N12) was found: `mdl_select_library` at mdl.c:247-248 has an O(genome_len/8) bitmap allocation that is the same scalability class as ENG-N2 but not covered by ENG-N2's scope. This should be bundled with the ENG-N2 patch.

This is round 1 of "new items" (ENG-N12 found). R5 is required to reach the 2-consecutive-round convergence criterion. ENG-N12 is a scope addition to an existing item (ENG-N2), not a new structural concern, so if R5 finds nothing new, convergence is reached.
