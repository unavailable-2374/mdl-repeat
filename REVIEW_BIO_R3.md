## Round 3 Biology Review

**Reviewer perspective**: TE biology / plant genomics
**Review date**: 2026-04-28
**Source documents read**: QUALITY_PROPOSAL_v3.md, QUALITY_PROPOSAL_v2.md, REVIEW_BIO_R1.md, REVIEW_BIO_R2.md

---

### Integration check (R2 items in v3)

**BIO-N7 (test params 20+80, internal-l-mer seeding)**: Correctly integrated.

v3 Tier 2 BIO-N2 revised spec precisely reflects the R2 requirement: 20 full-length elements, 80 solo LTR copies, and — critically — the internal-sequence l-mer seeding requirement is explicitly stated: "full-length element must be seeded by an INTERNAL-sequence l-mer (distinct from LTR sequence), otherwise pipeline produces only solo-LTR consensus." This is the core structural safeguard R2 asked for. The biological rationale (solo:full ratio 3:1 to 6:1 in plant gypsy elements per Wicker 2007 and Heitkam & Schmidt 2020) is preserved in the changelog header. No distortion.

One minor verification: v3 says "≥20 full-length copies" and "80 solo LTR copies" without a stated minimum for solo. R2 asked for "≥60 solo LTR copies" and ≥20 full-length. v3 states 80 solo (which exceeds 60) and 20 full-length (exactly at the minimum). This is consistent with R2's intent. Accept.

**BIO-N8 (Stage B Fix 1 baseline status)**: Correctly integrated.

v3 locked baseline section explicitly states: "Stage B Fix 1 (MDL standalone-fallback admit gate) ✅ shipped — accepts families with standalone_savings > 0 AND consensus_length ≥ 50 AND num_instances ≥ 3." The sentence further clarifies that the −0.241 bp_rec drop in REFINE_TRACE_REPORT was BEFORE this fix, and the current baseline already includes it. This directly resolves the R2 ambiguity: the −0.241 drop is now confirmed post-fix residual, meaning it represents the hard structural greedy-coverage problem, not a pending implementation gap. BIO-N8's core ask (clarify baseline status so proposer does not skip the biggest bottleneck) is satisfied. Accept.

**BIO-N9 (maize/wheat advisory)**: Correctly placed and expanded.

v3 adds a dedicated Tier 1.7 "Cross-genome calibration documentation" section that addresses all three BIO-N9 sub-points:
- MINTHRESH scaling advisory (3 may be too low for maize 2.3 Gb / wheat 17 Gb; suggests `max(3, log10(N) - 5)`) — present.
- ALIGN_MAX_EXTENSION unconditional bump for maize/wheat — present as J' amendment with explicit "ship the bump unconditionally."
- Coalesce factor 20× too aggressive for nested-stack genomes; recommend 5× for wheat/maize via ENG-N7 flag — present.
- E'' LTR-pair detection window configurable 15-30 kb — present.

The J' amendment ("ship the bump unconditionally") correctly reverses the v2 conditional logic. The biological reasoning is sound: even if Arabidopsis consensuses don't hit the 10 kb cap, Arabidopsis is not the design target's worst case; maize/wheat element sizes make the unconditional bump a universally safe +0.5-day insurance policy with no known downside for Arabidopsis. Accept unconditional J'.

**BIO-N10 (post-MDL assertion in BIO-N2 spec)**: Correctly placed.

v3 BIO-N2 revised spec includes: "Assertion (per BIO-N10): BOTH families must be in the FINAL post-MDL library, not just pre-MDL candidates. This simultaneously validates B's gate fix AND Stage B Fix 1 (standalone-fallback admit)." The failure-mode disambiguation is also present: "If only one family appears, distinguish whether it's bug B (merge collapse) or MDL select (fallback didn't fire)." This is precisely the assertion R2 asked for. Accept.

---

### Verification of revised items

**D' (rep-direction guard)**: Accept.

v3 explicitly adds a "representative-direction concern" to the D' spec: "explicitly mark Library 1 families as 'protected' (preferred representative) during the D' merge call. OR enforce the absorption direction by family-id ordering (Library 1 IDs < Pass 2 IDs)." This is biologically sound — a Pass 2 candidate with slightly more instances than its Library 1 counterpart should never absorb the already-refined Library 1 consensus, because the Library 1 consensus has accumulated multiple refinement iterations; replacing it with an unrefined Pass 2 seed degrades quality. The two implementation options given (protection flag vs. ID ordering) are both valid; the ID-ordering approach is simpler and less likely to introduce new bugs. No further biological concern. Accept.

**F' (BLAST minus-strand + PADLENGTH)**: Accept with one standing note.

v3 correctly inherits the R2 Eng concerns (minus-strand cons_start/cons_end mapping; PADLENGTH offset) as spec clarifications. The R2 bio note about MITE TIR boundary misrepresentation via BLAST HSP chaining is present and correctly labeled "Future concern; not blocking." This is the right disposition — it should not be elevated until A+BIO-N1 confirms MITEs are a major missed class. No new biological concern beyond what is already tracked. Accept.

**J' (unconditional bump)**: Accept.

The "conditional → unconditional" reversal is biologically justified by BIO-N9 (maize/wheat baseline). The specific framing in v3 ("even if Arabidopsis doesn't need it, ship the bump unconditionally") is correct. The verification gate ("verify chr4 + TAIR10 don't regress") is the right safeguard. No new concern.

---

### NEW action items (Round 3 only)

**No new biology action items in Round 3. Proposal is biologically converged from my perspective.**

Rationale: All four R2 items (BIO-N7 through BIO-N10) are faithfully integrated. The rescoped items (D', F', J') are biologically sound as written. The Tier 1.7 cross-genome documentation section fully addresses the non-Arabidopsis concerns raised in BIO-N9. The BIO-N2 revised spec is rigorous enough to catch the failure modes that R2 identified. I have no new biological concerns with the v3 wording or design.

---

### Open questions remaining

The following open questions from R2 remain unanswered because they require empirical data from the project, not from further review:

1. **What TE class are the 24 missed families?** (gates A+BIO-N1; unanswerable without BLAST against Dfam/Repbase on the extracted sequences)
2. **Does any chr4 library consensus reach 10 kb?** (Q2 awk; one command; still unlocked)
3. **Are ATHILA intervals in the recall denominator or excluded with CEN180?** (BIO-N6; positional verification against truth BED stratification)
4. **Run-to-run variance at -threads 4?** (ENG-N4; must be measured before any 1–2 pp delta claim is trusted)
5. **What is the solo:full-length ratio for missed families if any are solo LTRs?** (needed to confirm BIO-N2 test parameters are calibrated to real data after A+BIO-N1 result)

These are not items requiring further review rounds — they are data collection tasks that the proposer must execute. The proposal correctly lists them in its "Open questions remaining after R2" section.

**Stopping criterion assessment**: This is Round 3 with no new items from me. Combined with the Round 2 summary ("close to convergence; BIO-N7/N8/N10 are low-effort clarifications"), the biology review track has reached the 2-consecutive-round criterion for convergence, assuming the engineering reviewer also finds no new items in their Round 3. The proposal is ready to execute from a biological standpoint.
