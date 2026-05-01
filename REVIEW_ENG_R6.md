## Round 6 Engineering Review

_Reviewer: independent bioinformatics code review_
_Files read: REVIEW_ENG_R1.md through REVIEW_ENG_R4.md, QUALITY_PROPOSAL_v5.md, QUALITY_PROPOSAL_v6.md, src/refine.c (2554 lines), src/mdl.c (393 lines), diff v5 vs v6_

---

### Doc-fix verification (R5 items in v6)

**Doc-1 (Phase 2 execution order — add ENG-N11 to step 6, add ENG-N12 to step 8, bump effort 1 to 1.5 days, update Phase 2 total 2 to 2.5 days): FIXED**

v6 Phase 2 execution order:

- Step 4 (ENG-N9.aux, 0.5 day) removed; steps renumbered; note appended: "ENG-N9.aux audit COMPLETED in R3 — removed from execution order"
- Step 6 now reads "B + Q6 + ENG-N11 (1 day): fix all multi-chr unsafe sites in refine.c (`nested_containment_fraction` + `check_instance_overlap`); tighten nested gate predicate; both threading paths"
- Step 8 now reads "ENG-N2 + ENG-N10 + ENG-N12 (1.5 days): refine_prune + mdl_select_library coverage memory (both O(genome_len)) + qsort prune order"
- Phase 2 header changed from "2 days" to "2.5 days"
- Test-coverage table row updated from "ENG-N2 | malloc cap regression test (large genome simulation)" to "ENG-N2+N10+N12 | malloc cap regression test covering BOTH refine.c:1769 AND mdl.c:248 allocations; n_accepted >= 10k families to verify qsort scales"
- Separate "ENG-N10" and "ENG-N9.aux" rows removed from the test table; "ENG-N11" row added

All of Doc-1's required changes are present and internally consistent.

**Doc-2 (ENG-N9.aux removed from execution order): FIXED**

v5 Phase 1 listed "4. ENG-N9.aux (0.5 day): grep + audit all positional-comparison sites in refine.c" as a live work item. v6 removes this step number entirely and adds the inline note "ENG-N9.aux audit COMPLETED in R3 — removed from execution order" at the end of Phase 1. Confirmed by diff: the ENG-N9.aux step line is absent from v6's Phase 1. No trace of ENG-N9.aux remains in the execution order or test-coverage table.

**Doc-3 (audit count corrected in ENG-N9.aux section): FIXED — with one residual clarity nit**

v5 ENG-N9.aux body read: "5 unsafe: B+Q6 (lines 225, 230), ENG-N11 (lines 268, 275, 284, 307, 311)". This was wrong: it conflated 7 unsafe sites under the label "5 unsafe." v6 splits the line into two separate bullets (2 covered by B+Q6 + 5 covered by ENG-N11) and adds: "Total: 8 unsafe sites pre-fix (now all addressed by B/Q6/N8/N11) + 2 safe-by-design."

The arithmetic is now correct: 2 (B+Q6) + 5 (ENG-N11) + 1 (ENG-N8) = 8 unsafe, plus 2 safe-by-design = 10 total. Confirmed against live `grep -n '.position' src/refine.c` which returns exactly 10 sites at lines 225, 230, 268, 275, 284, 307, 311, 1776, 1834, 2474.

Residual clarity nit: the sentence "Total: 8 unsafe sites pre-fix (now all addressed by B/Q6/N8/N11)" elides ENG-N9's role. ENG-N9 does not fix a `.position`-comparison site; it fixes the sweep-line sort key in `InstanceEntry` (refine.c:1888-1895). The parenthetical is technically correct as written (B, Q6, N8, and N11 do cover all 8 unsafe `.position` sites), but omitting N9 from the parenthetical could confuse an implementer who asks "where does N9 fit in the audit?" This is a clarity concern, not an arithmetic error. The count 8+2=10 is accurate.

---

### NEW action items (Round 6 only)

**ENG-N13 [QUALITY_PROPOSAL_v6.md:1]: File title heading still reads "v5"; the change-log header describes "Changes from v4", not "Changes from v5"**

The file on disk is named QUALITY_PROPOSAL_v6.md but its H1 heading reads "mdl-repeat Output Quality Enhancement Proposal — v5" (line 1) and its change-log block opens with "Changes from v4" (line 3). The three Doc-1/Doc-2/Doc-3 fixes are correctly applied to the body, but there is no "Changes from v5" entry summarizing what R5 found and what v6 fixed. If the file is later referenced by reviewers or cited as v6, the header will cause confusion about which version is authoritative.

Fix: update line 1 to "v6"; add a "Changes from v5" paragraph above the existing "Changes from v4" block listing the three doc fixes (Doc-1: execution order updated with ENG-N11 bundled into step 6, ENG-N12 bundled into step 8, effort bumped to 1.5 days, Phase 2 total to 2.5 days; Doc-2: ENG-N9.aux removed from execution order; Doc-3: audit count corrected to 8 unsafe + 2 safe-by-design).

Severity: documentation self-labeling only. Does not affect any implementation spec, code-line citation, fix direction, test requirement, or effort estimate. The proposal body is internally correct.

---

### Final convergence assessment

R5 raised three doc-only items (Doc-1, Doc-2, Doc-3). All three are correctly applied in v6's body. One new item (ENG-N13) is found in Round 6: the file's title heading and change-log block were not updated from "v5" / "Changes from v4" to "v6" / "Changes from v5." This is a header self-labeling error only.

ENG-N13 has no effect on any implementation spec and requires no code review to close — it is a two-line header edit. The author should decide:

- If the convergence criterion requires zero new items of any kind: R7 is needed solely to confirm the ENG-N13 header edit was applied. No technical re-review of any code or spec is required.
- If the convergence criterion allows "zero new engineering-correctness items" (header label errors excluded): **Engineering review CONVERGED — R5 and R6 both found 0 new implementation-correctness items.**
