## Round 3 Engineering Review

_Reviewer: independent bioinformatics code review_
_Files read: REVIEW_ENG_R1.md, REVIEW_ENG_R2.md, QUALITY_PROPOSAL_v3.md,
QUALITY_PROPOSAL_v2.md, src/refine.c (all lines), src/types.h,
src/candidates.h, src/align.h, src/genome.h_

---

### Integration check (R2 items in v3)

**ENG-N8 (coalesce seq_index): correctly integrated.**
v3 §ENG-N8 correctly identifies the bug at refine.c:2487-2496, correctly
specifies the fix (`if (cur->seq_index != active->seq_index) { active = cur;
continue; }` before the gap test at line 2494), and correctly marks it a
BLOCKER for multi-chromosome inputs. One precision note: the v3 fix spec says
"add seq_index to the sort key" but the current sort (cmp_inst_by_start,
line 2446-2453) only compares `.start`. If seq_index is added to the sort key
the fix is trivially correct because cross-chromosome pairs will never be
adjacent in the sorted order. If the seq_index check is added as a runtime
guard WITHOUT changing the sort key, the guard still works correctly (it fires
whenever cur and active are from different chromosomes, regardless of
adjacency). Either approach is sound; the proposal should pick one and be
explicit. This is a specification ambiguity, not a correctness error in the
proposed fix.

**ENG-N9 (assemble seq_index): correctly integrated, with one residual
concern about the break vs. continue choice.**
v3 §ENG-N9 correctly identifies the sweep-line bug at refine.c:2000-2007 and
correctly notes that `InstanceEntry` (line 1864-1870) lacks a `seq_index`
field. The fix direction (add seq_index to InstanceEntry, populate at
line 1946, update cmp_instance_entry to sort by (seq_index, start)) is
correct. The break-vs-continue question raised in R2 is relevant: v3 says
"guard inner loop" without resolving it. Looking at the actual code: the inner
loop at line 2004-2007 `break`s when `entries[j].start - entries[i].start > D`
because entries are sorted by `.start`. If the sort is changed to
`(seq_index, start)`, then all entries from chromosome A come before all
entries from chromosome B. When the inner loop reaches the first entry on the
next chromosome, `entries[j].start - entries[i].start` will be the full
chromosome-boundary offset — which will be >> D for any real genome. So the
existing `break` will fire naturally at the chromosome boundary, and no
explicit `seq_index != seq_index` guard inside the inner loop is needed once
the sort is updated to `(seq_index, start)`. v3 should clarify this: the fix
is entirely in the sort comparator; no new guard inside the loop body is
needed.

**ENG-N9.aux (audit remaining position-comparison sites): PARTIALLY CORRECT —
one site is missing.**
v3 states "grep for all raw-`.position` comparisons in refine.c; for each,
confirm whether cross-chromosome safety is needed."

Running `grep -n '\.position' src/refine.c` returns exactly 10 sites:

- refine.c:225 — nested_containment_fraction outer loop (covered by B+Q6 fix)
- refine.c:230 — nested_containment_fraction inner loop (covered by B+Q6 fix)
- refine.c:268 — check_instance_overlap sorted path, outer loop
- refine.c:275 — check_instance_overlap sorted path, binary search (b_end)
- refine.c:284 — check_instance_overlap sorted path, scan forward (b_start)
- refine.c:307 — check_instance_overlap linear path, outer loop
- refine.c:311 — check_instance_overlap linear path, inner loop
- refine.c:1776 — refine_prune_families coverage array inner loop
- refine.c:1834 — refine_prune_families coverage decrement inner loop
- refine.c:2474 — refine_coalesce_tandem_instances sort key (covered by ENG-N8)

**`check_instance_overlap` (lines 252-334) is NOT listed in v3's ENG-N9.aux
audit and it uses raw `.position` comparisons on instances from two different
families, without any `seq_index` guard, across both the sorted path
(lines 267-301) and the linear-scan fallback (lines 306-328).**

This function is called in the merge gate. If two families each have instances
on chromosome A near its end AND chromosome B near its start, those instances
are adjacent in raw-position-sorted order. In the sorted path, the binary
search finds a cross-chromosome `b` instance whose `b_end` falls just above
`a_start` (because the padded raw coordinates are numerically close), then the
forward scan counts it as an overlap. In the linear-scan fallback the same
false overlap fires without any distance bound at all. The result is that the
merge gate returns "significant overlap" for two families that actually have
instances on different chromosomes — which suppresses a merge that might be
desirable, or triggers one that is not.

For large chromosomes (each >> PADLENGTH) this is unlikely to fire in
practice. For short contigs where PADLENGTH separates real sequence segments
by only a few hundred bases, it can fire. This is the same category as ENG-N8
and ENG-N9.

This site must be added to the ENG-N9.aux audit list as a new finding (see
ENG-N11 below).

**refine_prune_families at lines 1776 and 1834** uses `.position` as an index
into `cov[]` (the genome-wide coverage array). Since `cov` is allocated as
`genome_len` bytes and indexed by raw padded `.position`, these sites are
correct-by-design: the whole function is structured around the padded global
coordinate. No cross-chromosome bug here; both sites are safe as written.

**ENG-N10 (qsort prune): correctly integrated.**
v3 §ENG-N2 states "Bundle with ENG-N10 (NEW R2): replace the O(n^2) selection
sort at refine.c:1757-1766 with qsort. Both touch refine_prune_families so
combined patch." The bundling is sound and the effort increase (+0.5 day) is
reasonable. Confirmed: lines 1757-1766 contain the selection-sort loop as
described in R2.

**Q4 finding (multipass.sh = residual masking): correctly integrated.**
v3 §D' explicitly states "Q4 confirmed: existing multipass.sh uses RESIDUAL
masking (the OLD v1 design that bio rejected)" and correctly frames D' as "a
genuine rewrite, not a script tweak." The consequence (ENG-N1 prerequisite,
3-day effort) is accurately carried forward.

**Q6 finding (nested cross-chr bug): correctly integrated.**
v3 §B+Q6 correctly specifies the fix at refine.c:229 as `if
(shorter->instances[i].seq_index != longer->instances[j].seq_index) continue;`
inside the inner `j` loop. Confirmed against the actual code: line 229 is the
inner loop start (`for (int j = 0; ...`), so the guard would need to be
inserted as the first statement at line 230, immediately before the raw
`l_start` read. This is exactly right.

---

### Position-comparison audit (verify ENG-N9.aux scope)

`grep -n '\.position' src/refine.c` returns 10 results. Analysis of each:

- **refine.c:225**: `nested_containment_fraction`, outer loop reads `shorter->instances[i].position`. Unsafe across chromosomes. **Covered by B+Q6 fix.**

- **refine.c:230**: `nested_containment_fraction`, inner loop reads `longer->instances[j].position`. Unsafe across chromosomes. **Covered by B+Q6 fix.**

- **refine.c:268**: `check_instance_overlap`, sorted path, outer-loop reads `a->instances[i].position`. Unsafe: no seq_index guard. **NOT covered in v3. Raised as ENG-N11.**

- **refine.c:275**: `check_instance_overlap`, sorted path, binary-search reads `b_sorted[mid].position`. Unsafe: same function, same gap. **NOT covered in v3. Raised as ENG-N11.**

- **refine.c:284**: `check_instance_overlap`, sorted path, forward-scan reads `b_sorted[j].position`. Unsafe: same function. **NOT covered in v3. Raised as ENG-N11.**

- **refine.c:307**: `check_instance_overlap`, linear-fallback, outer loop reads `a->instances[i].position`. Unsafe. **NOT covered in v3. Raised as ENG-N11.**

- **refine.c:311**: `check_instance_overlap`, linear-fallback, inner loop reads `b->instances[j].position`. Unsafe. **NOT covered in v3. Raised as ENG-N11.**

- **refine.c:1776**: `refine_prune_families`, index into genome-wide `cov[]` array allocated over the full padded genome. **Safe by design** — the entire function uses padded global coordinates intentionally.

- **refine.c:1834**: Same function as above, decrement path. **Safe by design.**

- **refine.c:2474**: `refine_coalesce_tandem_instances`, sort key construction. **Covered by ENG-N8 fix.**

**Summary**: v3's ENG-N9.aux audit list is incomplete. `check_instance_overlap`
(refine.c:252-334) has five raw-position comparison sites across two code
paths, none guarded by seq_index, and none listed in the v3 audit.

---

### Verification of revised items

**D' (representative-direction guard at refine.c:886-892 union-find): accept.**
v3 correctly identifies the risk and adds "explicitly mark Library 1 families
as 'protected' (preferred representative) during the D' merge call OR enforce
by family-id ordering (Library 1 IDs < Pass 2 IDs)." Both approaches are
sound. The family-id-ordering approach is simpler (no new struct field) and
should be the default recommendation. The current union-find at refine.c:886-892
compares `num_instances` only; if Library 1 families are pre-labeled with IDs
in [0, N1) and Pass 2 families with IDs in [N1, N1+N2), a tiebreak on `id`
(smaller id wins) gives the desired protection with zero structural change.
v3's spec is acceptable.

**F' (BLAST minus-strand cons_start/cons_end + PADLENGTH offset): accept with
one remaining unresolved item.**
v3 correctly adds "Test required to assert correct mapping" for minus-strand
hits and correctly states "Must add PADLENGTH or silent off-by-N bug" for the
position field. Both are now captured as spec items rather than gaps. What
remains unresolved (flagged in R2, still not resolved in v3): the proposal
does not state WHERE the PADLENGTH addition occurs in the integration. The
only correct place is at the point where `inst->position` is set from the
BLAST-reported genomic coordinate, i.e., in the BLAST-output parser, not in
the BED writer. If the parser is written first without the PADLENGTH offset and
the BED writer is tested against the trace-dir path (main.c:104: `gpos_t
raw_start = inst->position - PADLENGTH`), the bug would appear as zero-based
coordinates that are off by PADLENGTH in both the BED and the internal
Instance. This is a subtle ordering risk in test design, not in the fix itself.
The v3 test requirement "Short-element BLAST: 0-based BED, PADLENGTH offset,
minus-strand cons_start/cons_end" is correct. Accept with the note that the
test must explicitly compare BED chromosome-local coordinates against the
expected genomic position, not just check that `inst->position` is non-zero.

**J' (unconditional ALIGN_MAX_EXTENSION bump): accept.**
v3 §J' correctly changes the R2 "skip if Arabidopsis doesn't need it" to
"ship unconditionally for universal-tool requirements." The rationale (maize
and wheat have routine 10-19 kb retrotransposons) is sound. The chr4 +
TAIR10 non-regression requirement is carried forward. Accept.

---

### NEW action items (Round 3 only)

**ENG-N11: `check_instance_overlap` contains five raw-position comparisons
across two code paths with no seq_index guard — must be added to ENG-N9.aux
audit.**

refine.c:252-334. Both the sorted path (lines 267-301) and the linear-fallback
(lines 306-328) compare raw `.position` values from two different families'
instance arrays. There is no `if (a->instances[i].seq_index !=
b->instances[j].seq_index) continue;` guard in either path. On a
multi-chromosome genome where chromosome boundaries are separated by only
PADLENGTH bytes in the concatenated buffer, an `a` instance near the end of
chromosome A and a `b` instance near the start of chromosome B can produce a
spurious overlap count. The merge gate at the call sites (refine.c:557 and
748) decides whether to block a merge based on this overlap count: a false
positive suppresses a legitimate merge; a false negative admits an incorrect
merge. For large chromosomes (>10 Mb each) the practical risk is low; for
short-contig assemblies it is real.

Severity: same class as ENG-N8 and ENG-N9. For a paper claiming universality
across genome types, this must at minimum be documented in ENG-N9.aux. If B+Q6
and ENG-N8+N9 are being patched anyway, the same `seq_index` guard should be
added to `check_instance_overlap` in the same patch.

Fix: in the linear-fallback inner loop (line 310), add
`if (a->instances[i].seq_index != b->instances[j].seq_index) continue;`
as the first body statement. In the sorted path, `b_sorted` is a copy of
`b->instances` sorted by raw position, so the seq_index field is preserved;
add the same guard at line 283 (the `for (int j = lo; ...)` body) before the
`b_start >= a_end` break check.

Effort: 15 minutes if done as part of ENG-N8+N9 combined patch; already
budgeted infrastructure (seq_index is in Instance struct).

---

### Open questions remaining

1. (Unchanged from R1/R2) At which pipeline stage do the 24 missed chr4
   families drop out? Entire priority ordering of B through E depends on this.

2. (Unchanged) Is the 0-split observation on chr4 explained by the
   REFINE_MIN_SPLIT_INSTANCES=10 floor or genuine unimodality? C
   instrumentation answers this.

3. (Unchanged) Run-to-run variance with -threads 4? ENG-N4 must run before
   any recall-delta claim is trusted.

4. (Unchanged) In D', when Library 1 + Pass 2 merged, which side becomes
   union-find representative? v3 specifies the fix (id-based tiebreak or
   protected flag); the implementation must choose one and document it before
   coding.

5. (Unchanged) ENG-N9 break-vs-continue: v3 says "guard inner loop" without
   resolving whether a seq_index check is needed INSIDE the loop body or
   whether sorting by (seq_index, start) alone makes the existing break
   sufficient. As analyzed above, sorting by (seq_index, start) alone IS
   sufficient; no new in-body guard is needed. This should be clarified in
   the v3 spec to avoid double work.

