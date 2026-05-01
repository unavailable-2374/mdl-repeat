## Round 4 Biology Review

**Reviewer perspective**: TE biology / plant genomics
**Review date**: 2026-04-28
**Source documents read**: REVIEW_BIO_R1.md, REVIEW_BIO_R2.md, REVIEW_BIO_R3.md, QUALITY_PROPOSAL_v4.md, REVIEW_ENG_R3.md

---

### Verification: bio side unchanged in v4 vs v3

v4's header states explicitly: "Bio R3 = 0 new items (1st converged round). Eng R3 = 1 new item ENG-N11 + 2 spec clarifications." The change description is narrow and purely engineering.

Checking each biology-relevant section in v4 against the v3 content I reviewed in R3:

- **A+BIO-N1, A.aux**: identical wording.
- **BIO-N2 revised spec** (Tier 2): identical. 20 full-length / 80 solo, internal-l-mer seeding requirement, and the BIO-N10 post-MDL assertion are all present and unchanged.
- **J' (unconditional bump)**: identical. The "ship unconditionally" decision from BIO-N9 is intact.
- **Tier 1.7 cross-genome calibration** (BIO-N9 items): identical. MINTHRESH scaling note, coalesce-factor advisory, E'' window configurability, all present and unchanged.
- **BIO-N3, BIO-N5**: identical dispositions (conditional deferral and Tier 5 respectively).
- **G, H, I reject dispositions**: unchanged.
- **F' bio note** (MITE TIR boundary future concern): unchanged, still correctly labeled "Future concern; not blocking."
- **D' representative-direction guard**: unchanged.
- **Open questions 1–6**: unchanged. Question 7 (cross-chr sites beyond N8/N9/Q6) is now answered by ENG-N9.aux inline in v4 — this is an engineering resolution, not a biological change.

The only material differences between v3 and v4 are:
1. Addition of ENG-N11 (the `check_instance_overlap` five-site audit and fix spec) — engineering only.
2. ENG-N8 spec clarification (runtime guard alone, no sort-key change) — engineering only.
3. ENG-N9 spec clarification (sort by `(seq_index, start)` alone is sufficient; no in-body guard) — engineering only.
4. ENG-N9.aux section replaced with an inline completed audit table — engineering bookkeeping only.

None of these touch biology. **Bio sections are verified unchanged from v3.**

---

### NEW action items (Round 4 only)

**No new biology action items in Round 4. This is round 2 of "no new", satisfying the bio convergence criterion.**

Rationale: ENG-N11 is entirely within the merge-gate cross-chromosome correctness theme. Its biological consequence (false-positive or false-negative merge decisions for multi-chromosome inputs) was already covered in general terms by the B+Q6 and ENG-N8/N9 discussion in earlier rounds. The ENG-N11 fix closes the same gap in a third code path; it does not raise any new biological question about TE biology, element classification, benchmark interpretation, or evaluation scope. All biology-relevant items (BIO-N1 through BIO-N10) remain correctly represented in v4 with no distortion.

---

### Final assessment

Biology review is converged. Round 3 produced 0 new items; Round 4 produces 0 new items. The 2-consecutive-round stopping criterion is satisfied.

The biological foundations of the proposal are sound:
- TE class diagnosis gates the entire roadmap (A+BIO-N1) — correctly positioned as Tier 1.
- Solo LTR / full-element structural test (BIO-N2) is specified with biologically realistic copy-count ratios and the correct post-MDL assertion.
- MITE-seeding concern (BIO-N3) is correctly held conditional on A results.
- Cross-genome calibration notes (BIO-N9 / Tier 1.7) correctly flag known failure modes for large-genome use.
- All rejects (G, H, blanket-I) remain biologically justified.

Five open questions from R2/R3 still require empirical data (TE class of 24 missed families, chr4 max consensus length, ATHILA denominator verification, run-to-run variance, solo:full ratio). These are data collection tasks for the proposer, not items requiring further review.

**Biology review: CLOSED.**
