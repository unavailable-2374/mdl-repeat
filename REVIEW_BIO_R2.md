## Round 2 Biology Review

**Reviewer perspective**: TE biology / plant genomics
**Review date**: 2026-04-28
**Source documents read**: QUALITY_PROPOSAL_v2.md, QUALITY_PROPOSAL_v1.md, REVIEW_BIO_R1.md, HANDOFF.md, STAGE_A_REPORT.md, REFINE_TRACE_REPORT.md

---

### Integration check

**BIO-N1 in v2**: Correctly integrated. v2 collapses A and BIO-N1 into a single combined item (A+BIO-N1), lists the correct Wicker-class categories (solo LTRs, MITEs, TSI/pericentric tandem variants, TIR/Helitron, full-length low-chr4-copy elements), and cites Heitkam & Schmidt 2020 for copy numbers. The framing that this gates everything downstream is retained faithfully. BIO-N4 (TSI/CEN180 purity check) and BIO-N6 (ATHILA denominator verification) are correctly absorbed into the A.aux sanity block. The ENG-N1 controlled small-genome diagnostic is added here as a companion, which is reasonable. One minor weakening: v2 says "Effort: 1 day" for A+BIO-N1, where R1 said "1 day" for A+BIO-N1 combined — no distortion.

**BIO-N2 in v2**: Correctly integrated. v2 reproduces the synthetic ATHILA test design and frames it explicitly as a go/no-go gate for the 500 bp – 12 kb length-bin recovery. The framing that the two expected output families are (a) ~11 kb internal consensus and (b) ~600 bp LTR consensus is retained. One concern with the test parameters is discussed in the "test specs" section below (see BIO-N7).

**BIO-N3 in v2**: Correctly integrated, and appropriately demoted to Tier 4 (conditional on BIO-N1 identifying MITEs as a major missed class). The diagnostic description (check whether palindromic-self-RC seeds are chosen; look for high-frequency reverse-complement twin of the seed l-mer at similar genomic positions) accurately captures the R1 concern. The "deprioritize palindromic-self-RC seeds" fix approach is correct.

**BIO-N4 in v2**: Correctly integrated as part of the A.aux sanity block. The specific size concern (358 bp ≈ 2× CEN180 monomer; 470 bp ≈ 2.5×) is not restated explicitly, but the action "verify any of the 24 are CEN180/TSI fragments that should be removed from the recall denominator" captures the intent. No distortion.

**BIO-N5 in v2**: Correctly integrated as a separate Tier 5 polishing item ("Library-utility correctness" section). The biological rationale (canonical strand matters for K-divergence dating, CpG decay analysis, RepeatMasker annotation direction) is stated correctly. The implementation approach (BLAST against UniProt protein domains, flip if INT-RH-RT order observed in 5'→3' direction) is correct and faithful to R1. The 2–3 day effort estimate is reasonable for the domain-search integration step.

**BIO-N6 in v2**: Correctly integrated in the A.aux sanity block as "Verify ATHILA intervals are correctly inside the recall denominator (not excluded with CEN180)." The distinction between ATHILA (interspersed LTR retrotransposon, in-scope) and CEN180 (tandem satellite, out-of-scope) is preserved. The open question "Are ATHILA intervals counted in recall denominator or excluded with CEN180?" remains listed as an open question, which is appropriate since it has not yet been answered empirically.

---

### NEW action items (Round 2 only)

**BIO-N7: The BIO-N2 synthetic test parameters require revision — 8 + 12 is not the right copy count distribution for testing solo LTR handling.**

The v2 BIO-N2 spec says "8 copies of a full-length ATHILA-like element (~11 kb internal + 600 bp LTR each side = 12 kb total)" and "12 copies of the same 600 bp LTR sequence as solo elements." This ratio (8 full-length, 12 solo) is biologically inverted for a tool designed to target high-copy families: the solo LTR class at 12 copies and the full-length class at 8 copies will both fall near the MDL acceptance floor. The concern is that the test may pass or fail for the wrong reason.

The biological context matters here. In Arabidopsis, solo LTRs arise from intrastrand recombination between the 5' and 3' LTRs of a full-length element. For every full-length ATHILA element that undergoes this recombination, one solo LTR is left behind. Empirically, solo LTRs substantially outnumber full-length elements in plant genomes: Wicker et al. 2007 (PMID 17984973) and subsequent analyses of Arabidopsis report solo-LTR-to-full-length ratios of approximately 3:1 to 6:1 for gypsy elements. For ATHILA specifically, the Heitkam & Schmidt 2020 survey (PMC7385966) reports ATHILA2 with 413 total copies — if the solo:full ratio is ~4:1, that is approximately 80–100 solo LTRs and 20–25 full-length copies.

A biologically realistic test design for the BIO-N2 synthetic should therefore be:
- 20 copies of the full-length ATHILA-like element (12 kb total)
- 80 copies of the 600 bp solo LTR alone

At these copy counts, both families are well within the MDL acceptance zone and above the MINTHRESH floor. The 8/12 design risks passing the go/no-go gate vacuously because the tool may simply find whichever family happens to cross MDL acceptance by chance, not because the pipeline correctly handles the structural distinction.

An additional complication: the v2 spec says the full-length element is "~11 kb internal + 600 bp LTR each side = 12 kb total." For the purpose of testing pipeline behavior, the critical structural feature is that the 600 bp LTR sequence is SHARED between the full-length element (at both ends) and the solo LTR. If the solo LTR frequency (80 copies) greatly exceeds the full-length LTR occurrence (40 copies, 2 per full-length × 20 elements), the seeding step will preferentially choose an l-mer from the 600 bp LTR region, and the pipeline may generate a solo-LTR consensus only and never build a 12 kb full-length consensus. This is actually a real algorithmic failure mode that the test should be designed to catch — not hide. Set up the test so that the full-length internal sequence has its OWN high-copy l-mer distinct from the LTR, to ensure the full-length element is independently seeded.

Action: Revise BIO-N2 spec to use ≥20 full-length copies (seeded by an internal-region l-mer that does not appear in the LTR) and ≥60 solo LTR copies, verifying that both seed independently in the discovery log.

**BIO-N8: The REFINE_TRACE shows MDL select as the dominant bottleneck (−0.241 bp_rec), but v2 has no item addressing greedy-accept order bias for overlapping same-region families.**

The REFINE_TRACE_REPORT stage-by-stage recall curve shows:

- After discovery: bp_rec = 0.580
- After MDL select: bp_rec = 0.280 (loss of 0.241)
- After coalesce: bp_rec = 0.397

The v2 proposal addresses the merge-gate bug (B), the split firing (C), and instance recruitment (F'). None of these addresses the largest single recall loss step: MDL select killing 28% of families whose instances happen to overlap an already-accepted family. The Stage B Fix 1 (standalone MDL gate fallback) addresses this in part — accept if standalone_savings > 0 AND num_instances ≥ 3. This is listed in the REFINE_TRACE_REPORT and is presumably already in the baseline the v2 proposal is building on top of (the v2 locked baseline includes "Stage B" results). But the v2 proposal does not mention this fix explicitly, and its status relative to the locked baseline is ambiguous.

The biological concern is: when two biologically distinct but sequence-overlapping TE families (e.g., an LTR retrotransposon whose 3' LTR overlaps the next element's 5' LTR insertion site) are evaluated by the greedy unique-coverage MDL select, the family accepted first receives credit for the shared bp, and the second family's marginal savings are near zero even though both are valid. This is not a bug in the tool's design — it is an inherent property of the greedy formulation — but it is the primary recall gap.

The v2 proposal should explicitly state whether the Stage B standalone-MDL fallback (Fix 1 from REFINE_TRACE_REPORT) is already in the locked baseline or still needs to be implemented. If already in baseline, the residual −0.241 bp_rec at MDL select is post-Fix-1 and represents a harder structural problem. If not yet in baseline, implementing it is likely the single highest-impact remaining fix and should appear as a Tier 1 item in v2.

This is not a new algorithmic concern — it was present in the REFINE_TRACE_REPORT — but v2's silence on it is an integration gap that could cause the proposer to implement lower-impact items (C, F') while the largest bottleneck remains unaddressed.

Action: Clarify in v2 whether Stage B Fix 1 (standalone MDL fallback) is in the locked baseline. If not, elevate it to Tier 1 alongside B.

**BIO-N9: Non-Arabidopsis genome considerations — maize and wheat TE landscapes differ structurally in ways that affect several v2 design choices.**

v2 mentions the wheat genome (17 Gb) and maize only in the context of the ENG-N2 memory fix (calloc fatal on wheat). But several of the biological design choices in v2 were calibrated to Arabidopsis and may not generalize:

1. The MINTHRESH floor of 3 copies (item H, correctly rejected) is calibrated to the Arabidopsis genome (~135 Mb), where a 3-copy family with 500 bp consensus provides non-trivial MDL savings. In maize (2.3 Gb) and wheat (hexaploid, ~17 Gb), many TE families have hundreds to thousands of copies, and the minimum meaningful family is well above 3 copies. The 3-copy floor is not just a precision question — it is a library bloat question. In a 17 Gb genome, accepting 3-copy families would generate thousands of spurious consensuses from sequence-shared regions. This is correctly rejected in v2, but the rationale given (design scope) is weaker than the genome-size-scaled rationale: MINTHRESH should scale with genome size, and the current 3-copy floor may already be too low for large genomes even by the original design philosophy.

2. The ALIGN_MAX_EXTENSION discussion (J') focuses on ATHILA at ~11 kb. Maize LTR retrotransposons (CINFUL-Zeon, Grande, Opie-Ji) are commonly 10–18 kb internal region plus 0.2–0.5 kb LTRs for a total of 10–19 kb. Wheat TEs include elements (Wis, Barbara, Fatima) similarly large. If J' is being deferred pending the awk check of maximum consensus length in the chr4 output, the awk result showing "no consensus ≥ 10 kb on Arabidopsis" would incorrectly suggest J' is irrelevant — when in fact ALIGN_MAX_EXTENSION may be a hard ceiling for large-genome TE families regardless of Arabidopsis output. The J' conditional logic should include a statement that even if Arabidopsis does not trigger the limit, maize/wheat testing should use ALIGN_MAX_EXTENSION ≥ 20000 by default.

3. The coalesce factor ENG-N7 item (expose -coalesce-factor CLI) is motivated by "compact genomes (fungal, protozoan) where adjacent same-family copies are NOT in tandem arrays." But the more important motivation is the opposite direction: in maize and wheat, LTR retrotransposons form nested stacks (LTR inside LTR inside LTR), and a 20× coalesce factor would spuriously merge adjacent elements from different families into one mega-interval. The default of 20× should be flagged as potentially too aggressive for large nested-TE genomes, not just too conservative for small genomes. A suggested default for wheat/maize might be 5× (within-element gap closure) rather than 20× (cross-element bridging).

None of these are blockers for the Arabidopsis benchmark, but if the "universal tool" framing is genuine, they should be listed as design constraints.

Action: Add a brief "Genome-scope caveats" section to v2 noting the three points above. No new implementation required; these are design notes for when non-Arabidopsis benchmarking begins.

**BIO-N10: The interaction between Stage B standalone-MDL gate and the BIO-N2 solo LTR test is not specified — there is a likely failure mode that the test must cover.**

The BIO-N2 test expects TWO output families: one ~12 kb full-length consensus and one ~600 bp solo LTR consensus. But consider the expected MDL accounting:

- Full-length ATHILA family: 20 copies × 12 kb = 240 kb of genome sequence explained. Standalone MDL savings are large — this family clearly passes.
- Solo LTR family: 60–80 copies × 600 bp = 36–48 kb of genome sequence explained. But: every full-length ATHILA copy ALSO contains a 600 bp LTR at each end. That is 40 occurrences of the 600 bp LTR sequence that are already claimed by the full-length family's consensus. After the greedy MDL select accepts the full-length family, the solo LTR family's marginal savings = savings from the 60–80 SOLO copies ONLY = 36–48 kb. Its standalone savings are computed over all 100+ occurrences, but the marginal (post-greedy) savings are over 60–80. This is exactly the regime where the Stage B standalone fallback gate matters.

If the Stage B standalone fallback is NOT in the locked baseline, the BIO-N2 test will fail (the solo LTR family will be rejected as having near-zero marginal savings after the full-length family is accepted), and the failure will be misattributed to the merge gate (B) rather than the MDL select gate. The fix is not in B — it is in the MDL select logic.

This interaction must be explicitly acknowledged in the BIO-N2 test specification: the test should check not only whether two families are produced but also whether both PASS MDL select independently, and the test expectation must be informed by which fallback gates are active.

Action: Add an assertion to the BIO-N2 test: both families must appear in the final library (post-MDL-select), not just in the pre-select candidate list. The spec should note that this test simultaneously validates the merge gate fix (B) AND the MDL select standalone fallback gate.

---

### Verification of rescoped items

**D' (multi-pass full-genome seeding)**: Accept with one remaining concern. v2 correctly addresses R1's masking-boundary artifact concern: it specifies that Pass 2 seeds on the FULL genome (not the residual), and uses a BLAST screen at the merge stage to exclude Pass-2 candidates contained in Library 1. This is structurally identical to RepeatModeler2's multi-pass logic. The pre-condition (ENG-N1 controlled diagnostic must confirm seed-competition is a real bottleneck) is appropriate. One residual concern: the "BLAST screen at merge stage" to exclude Library-1-contained candidates requires a functional BLAST integration that does not currently exist in the codebase. v2 should clarify whether this step uses the existing 80-80-80 merge logic (which is in-process and sequence-based) or a separate BLAST call. If the latter, the effort estimate of 2–3 days is likely understated by approximately 1 day for the BLAST integration plumbing. This is an engineering concern, not a biology concern — accept the biological framing.

**E'' (structural features — TIR pair detection, LTR pair detection)**: Accept. The revised framing avoids the circularity problem identified in R1 (profile-from-Pass-1 reinforces existing families). Structural features (palindromic TIR windows, long-range LTR identity within ~15 kb) are genomic properties of the sequence itself, not derived from a previously discovered library, and therefore do not suffer from the R1 circularity concern. The deferral to "after Tier 1-3 results are in" is appropriate given the multi-week implementation cost. One note: the 15 kb window for LTR-pair detection is reasonable for Arabidopsis (ATHILA full elements are ~12 kb), but for maize/wheat with larger elements this window may need to be 25–30 kb. Flag for when implementation begins.

**F' (BLAST for short elements)**: Accept with one biological clarification. The revised framing (replace banded DP only for consensus_length < 500 bp; BLAST helps with short elements where fixed bandwidth handles indels poorly) correctly implements R1's reframing from "divergence" to "element length." The 500 bp cutoff is a reasonable first approximation. A biological note: for MITE-class elements specifically (100–500 bp, TIR-flanked), BLAST's HSP chaining tends to produce alignments that span both TIRs as a single hit, which can misrepresent the element boundary. If MITEs are confirmed as a major missed class by BIO-N1, the F' BLAST implementation should separately handle TIR-flanked elements to avoid generating artificially short or internally-truncated consensuses. This is a future concern, not a current blocker.

**J' (ALIGN_MAX_EXTENSION conditional)**: Accept, but with the caveat noted in BIO-N9 above. The "check if any chr4 consensus reaches 10 kb" gate is appropriate for Arabidopsis but should not be used to conclude J' is irrelevant for larger genomes. The conditional logic is correct for the current evaluation scope.

---

### Open questions remaining after v2

1. **Is the Stage B standalone-MDL fallback (REFINE_TRACE_REPORT Fix 1) already in the locked v2 baseline, or still pending?** This determines whether the −0.241 bp_rec at MDL select is already post-fix (hard structural problem) or a remaining high-impact fix. The v2 proposal is silent on this and the REFINE_TRACE_REPORT's "Stage B target" framing suggests it was a proposal, not a completed fix. If not yet implemented, this is the highest-priority remaining biology item, outranking all others.

2. **What is the TE class of the 24 missed families?** Still unanswered — this is the gate for the entire roadmap. No change from R1 open question 1.

3. **Does any chr4 library consensus reach 10 kb?** One `awk` command on the current output library. Still listed as open question Q2 in v2. Resolve before spending 0.5 day on J'.

4. **Are ATHILA intervals in the recall denominator or excluded with CEN180?** Still listed as open question 3 in v2. Affects whether the 0.465 in-scope recall is accurately measured.

5. **Run-to-run variance with `-threads 4`?** ENG-N4 is specified but not yet executed. Until this is measured, the significance of any 1–2 pp improvement claim is unverifiable.

6. **What is the solo:full-length copy ratio for the missed families if they are solo LTRs?** This directly determines the BIO-N2 test design (addressed in BIO-N7 above). If BLAST classification in A+BIO-N1 confirms solo LTRs, the ratio should be extracted from the chr4 RepeatMasker annotation (TAIR10 soft-mask source) before the synthetic test is constructed.

7. **Does the 20× tandem coalesce factor generate spurious mega-intervals in pericentromeric regions where ATHILA and TSI elements are adjacent?** This was raised in R1 open question 5 and is not addressed in v2. ENG-N7 adds a CLI flag for the coalesce factor but does not specify a test for the pericentromeric false-merge scenario. The benchmark stratification should include "coalesced vs non-coalesced" comparison for the inside-centromere precision numbers.

---

### Summary judgment on v2 convergence

All six BIO-N1..N6 items were integrated correctly and without distortion. The four rescoped items (D', E'', F', J') faithfully reflect the biological arguments from R1. Four new biology action items are raised (BIO-N7 through BIO-N10). Of these:

- **BIO-N7** (BIO-N2 test parameters) and **BIO-N10** (MDL-select / BIO-N2 interaction) are directly blocking: the BIO-N2 test as written may pass or fail for the wrong reason. These should be addressed before executing BIO-N2.
- **BIO-N8** (MDL select bottleneck status) is a tracking/clarity issue, not a new algorithm concern; it requires one clarification in v2 to resolve.
- **BIO-N9** (non-Arabidopsis genome considerations) is advisory for future work and does not block Arabidopsis benchmarking.

The proposal is close to convergence on biology. Resolving BIO-N7, BIO-N8, and BIO-N10 (all low-effort clarifications) would make it ready to execute.
