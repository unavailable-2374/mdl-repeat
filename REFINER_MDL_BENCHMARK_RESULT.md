# Refiner_mdl Benchmark Result — TAIR10 nuclear

**Date**: 2026-05-01
**Scope**: End-to-end empirical test of Refiner_mdl as a downstream polish for mdl-repeat
**Status**: **Refiner_mdl in current form NET HURTS family-level recall.** Recommended status: NOT a standard downstream for mdl-repeat.

---

## 0. TL;DR

We ran the full pipeline `mdl-repeat → Refiner_mdl` on Arabidopsis TAIR10 nuclear (119 Mb) and measured family-level recall against the project's standard `tair10_truth_clusters.fa` set (1517 in-scope ≥3-copy truth families).

| Pipeline | seqs | 80×50 | 80×80 | 90×80 | Δ vs baseline (80×80) |
|---|---:|---:|---:|---:|---:|
| **mdl-repeat baseline (v6_final)** | 16760 | **0.8583** | **0.8207** | **0.6605** | — |
| mdl-repeat v6.1 (R8 header only) | 16768 | 0.8583 | 0.8214 | 0.6618 | +0.0007 (noise) |
| min: hard_filter only | 13770 | 0.8438 | 0.8069 | 0.6467 | -1.4 pp |
| Phase1-only (R5+R6+R7+R8 polish, no te_struct) | 12522 | 0.8438 | 0.7963 | 0.6427 | -2.4 pp |
| **Refiner_mdl analysis lib (full pipeline)** | 854 | 0.2624 | 0.2307 | 0.1819 | **-59 pp** |
| Refiner_mdl masking lib (full pipeline) | 754 | 0.2512 | 0.2215 | 0.1701 | -60 pp |

**Refiner_mdl is NOT a recall-improving step.** All variants degrade recall; the full pipeline collapses it by ~60 pp.

---

## 1. Background

The previous session (FINAL_REPORT.md) closed with three open items at TAIR10 nuclear scale:
- 90×80 = 0.661 (40% of in-scope families lack a full-length 90%-identity consensus)
- BIO-N5 canonical strand orientation deferred (UniProt domain DB needed)
- F' BLAST short-element recruitment reverted (assembly O(n²) blowup)

Refiner_mdl (Pan_TE/bin/Refiner_mdl/) was a candidate downstream that, on paper, should have addressed all three by:
- Phase 1: BLASTN+MAFFT consensus rebuild → close 90×80 gap
- Phase 2 te_structure_filter: protein-domain detection → BIO-N5 substitute
- Phase 0.4: BLAST-based fragment assembly → F' substitute (without the assembly cliff)

This benchmark tests whether Refiner_mdl actually delivers any of those gains.

---

## 2. Setup

**Input**: TAIR10 nuclear, 5 chromosomes, 119,157,348 bp
**Truth set**: `/tmp/ath_bench/family_eval/tair10_truth_clusters.fa` — RepeatModeler2-style cd-hit clusters of TAIR10 RepBase TE annotations (31,206 clusters; 1,517 in-scope ≥3-copy)
**Eval**: `family_recall_tair10.py` — dc-megablast truth representatives → library, then per-cluster identity×coverage thresholds

**Two mdl-repeat outputs tested**:
1. `tair10_v6f.fa` (16760 families) — pre-existing v6_final, FASTA header lacks `div=` / `topo=`
2. `tair10_v61.fa` (16768 families) — built with the v6.1 mdl-repeat (R8 header expansion: `div=%.3f topo=%s` added, fully backward compatible)

**Refiner_mdl** invoked through `python3 -m Refiner_mdl.main` with PGTA conda env (blastn 2.16, MAFFT 7.520, cd-hit-est 4.8.1, samtools, RepeatPeps.lib + RepeatMasker.lib from RepeatMasker 4.1.8).

---

## 3. Per-stage marginal contribution to recall

```
mdl-repeat output ............................ 16768 (baseline 0.8214 at 80×80)
  ↓ hard_filter (length / copies / N% / entropy / DUST / cyclic)
                                                13770 (-1.4 pp)
  ↓ Phase 0.4 fragment_assembly (directed greedy chain, BED-based)
                                                13750 (≈0)         16 scaffolds, only 36 fragments merged
  ↓ Phase 0.5 cd-hit-est dedup (95% / 90%)
                                                12523 (-1.0 pp)
  ↓ Phase 1 BLASTN→MAFFT consensus polish + chimera split
                                                12524 (≈0)         3304/12524 = 26.4% sequences got MAFFT consensus, but no measurable recall change
  ↓ Phase 2 QC (length / entropy / DUST)
                                                12522 (≈0)
  ↓ Phase 2 te_structure_filter (protein / TIR / LTR / known TE)
                                                ≈ 854 (-57 pp) ★ catastrophic
                                                12522 → 854: 93% of post-Phase-1 sequences fail all four structure tests
```

---

## 4. Why te_structure_filter is incompatible with mdl-repeat

`te_structure_filter.py` keeps a sequence iff it has at least one of:
- protein domain hit (blastx vs RepeatPeps.lib)
- TIR (terminal inverted repeat, self-blast minus strand)
- LTR (long terminal repeat, self-blast plus strand)
- known TE homology (blastn vs RepeatMasker.lib)

This silently assumes every true TE shows one of these signals. But mdl-repeat's design scope (CLAUDE.md) explicitly includes:
- **Helitrons** — no classical TIR
- **MITEs** — short, no protein domain, weak TIR signal
- **SINEs** — no protein domain
- **Novel/diverged families** — no hit in RepeatMasker.lib

These are exactly the family classes that fall through all four filters.

The two TE definitions don't overlap:

| Criterion | mdl-repeat | Refiner_mdl te_structure_filter |
|---|---|---|
| copies ≥ 3 | required | ignored |
| MDL savings > 0 | required | ignored |
| has protein/TIR/LTR/known homology | ignored | required |
| length ≥ 50 bp | required | required |

---

## 5. Verdict on the multi-round-iterated R5/R6/R7/R8 contract changes

These four changes were implemented (see this session's git diff for `src/output.c` and `Refiner_mdl/{config,phase0_triage,phase1_consensus}.py`) and are correct in their own right, but on TAIR10 they produced no measurable recall improvement:

| Change | Implementation | Empirical effect on TAIR10 |
|---|---|---|
| R5: Phase 1 BLASTN min_pident 80→70 | Refiner_mdl/config.py | Phase 1 itself contributes ≈0 to recall, so R5 is silent |
| R6: cd-hit dedup sort key (mdl, length) | Refiner_mdl/phase0_triage.py | Operates inside dedup, which costs -1.0 pp; R6 is a sub-optimization with no measurable contribution |
| R7-cyclic: drop topology=cyclic in hard_filter | Refiner_mdl/phase0_triage.py | TAIR10 nuclear has 0 cyclic families (no plasmid/organelle); R7-cyclic does not fire |
| R7-complex: force chimera detection on topology=complex | Refiner_mdl/phase1_consensus.py | TAIR10 has 0 complex families; does not fire |
| R7-div ceiling: T1→T2 demotion when div>0.20 | Refiner_mdl/phase0_triage.py | Demoted 3062 families. Phase1-only old (R5+R6) → 0.7943; Phase1-only new (R5+R6+R7+R8) → 0.7963. Δ=+0.0020 (+3 truth families), within noise |
| R8: mdl-repeat FASTA header `div=` `topo=` | mdl-repeat/src/output.c | Backward compatible; sanity-check confirms identical-to-baseline recall (0.8214 vs 0.8207, +1 family in noise) |

The R5/R6/R7/R8 work fixed real interface problems — but those problems didn't matter on TAIR10 because the dominant negative was downstream (te_structure_filter).

**Decision**: keep R8 (the mdl-repeat header expansion costs nothing, helps if Refiner_mdl is ever fixed) and keep R5/R6/R7 in Refiner_mdl (they're correct configurations if someone uses Refiner_mdl), but stop treating "Refiner_mdl as standard downstream" as the goal.

---

## 6. What this benchmark closes (and doesn't)

**Closes**:
- "Should mdl-repeat ship Refiner_mdl as the recommended downstream?" → **No.**
- "Will Refiner_mdl Phase 1 polish close the chr4 90×80 = 0.592 gap?" → **No.** Phase 1 contributes 0 ± noise on TAIR10 90×80 (baseline 0.6605 → 0.6427). The 90×80 gap on chr4 will not be closed by MAFFT consensus rebuild.
- "Are the contract changes (R5/R6/R7/R8) correct?" → **Yes**, all behave as designed. **Useful?** → only if someone fixes the te_structure_filter destructiveness; until then they're correct-but-silent.

**Doesn't close**:
- Whether a *different* downstream polish (not Refiner_mdl) could close the 90×80 gap. mdl-repeat's column-majority consensus is empirically already near-optimal vs MAFFT — but that's a 16 k-family aggregate observation; the 40% of families that fail 90×80 might still benefit from something more targeted (e.g., per-family banded re-extension with stricter stopping; conditional consensus rebuild only for div>0.15 families with copies>10).
- The genuine cost-benefit of Refiner_mdl when the metric is **library compactness** (not recall): a 754-sequence masking library can still mask 80%+ of the genome's TE bases if the families are well-chosen — that's a different evaluation we did not run.
- Whether Refiner_mdl is appropriate for a different upstream tool (RepeatModeler2 / EDTA / REPET output) where the input is already protein-validated. Likely yes — Refiner_mdl was built for that profile, not for mdl-repeat's broad inclusion of MITEs / Helitrons.

---

## 7. Recommendations (replaces the previous session's "next-step proposals")

### Stop doing
- Do not promote Refiner_mdl as a downstream for mdl-repeat.
- Do not redo "F' BLAST short-element recruitment" inside mdl-repeat assuming Refiner_mdl Phase 1 substitutes for it — Phase 1 doesn't substitute, it just doesn't do anything measurable on this dataset.
- Do not pursue R4 (disable Refiner_mdl Phase 0.4) or R9 (merge log emission) — they presumed Refiner_mdl was the right downstream, which is now disproven.

### Start doing (or consider)
- **Document in CLAUDE.md / FINAL_REPORT.md**: mdl-repeat's TAIR10 0.821 / 0.661 numbers are the recommended final library; do not pipe through Refiner_mdl by default.
- **If a "cleaner" library is needed for a specific downstream** (RepeatMasker masking, TE classification), apply only `Refiner_mdl.phase0_triage.hard_filter` (or its equivalent — DUST + entropy + min-length + min-copies). Cost: -1.4 pp at 80×80, removes 18% of marginal sequences.
- **For 90×80 = 0.66 → ~0.75 gap closure**: design a **mdl-repeat-internal** consensus refinement — conditional re-extension with banded DP at lower stopping threshold for high-divergence (>0.15) high-copy (≥10) families only. Test G fixtures (retained in tests/data/testG.{fa,_truth.bed}) provide a starting point.
- **If Refiner_mdl is to be made compatible with mdl-repeat output**, the project owner of Refiner_mdl should rewrite te_structure_filter to be a soft labeling step (annotate `te_signal=protein|tir|ltr|known|none`) rather than a destructive filter, and possibly add Helitron / MITE / SINE-specific structural probes.

### Already done in this session (locked, no rollback recommended)
- mdl-repeat `src/output.c` — FASTA header expanded with `div=%.3f topo=%s` (R8)
- Refiner_mdl `config.py / phase0_triage.py / phase1_consensus.py` — R5+R6+R7 implementation, regression tests pass (4 behavior assertions), 22+34+17 mdl-repeat tests pass
- This benchmark report

---

## 8. Artifacts

**mdl-repeat outputs**
- `/tmp/ath_bench/v6_final/tair10_v6f.{fa,bed,tsv}` — pre-existing baseline (FINAL_REPORT.md numbers)
- `/tmp/ath_bench/v61/tair10_v61.{fa,bed,tsv,log}` — v6.1 (R8 header) build, this session

**Refiner_mdl outputs**
- `/tmp/ath_bench/refiner_oldhdr/` — Refiner_mdl on tair10_v6f.fa (R5+R6 only)
  - `consensus_masking.fa` 754, `phase3_analysis_library.fa` 854, `phase1_only.fa` 12522, `refiner_mdl_stats.json`
- `/tmp/ath_bench/refiner_newhdr/` — Refiner_mdl on tair10_v61.fa (R5+R6+R7+R8)
  - `consensus_masking.fa` 724, `phase3_analysis_library.fa` 856, `phase1_only.fa` 12522
- `/tmp/ath_bench/refiner_min/` — minimal-pipeline variant (hard_filter only)
  - `all_post_hardfilter.fa` 13770, `consensus_masking.fa` 4993, `phase3_analysis_library.fa` 4993

**Family-level eval cache**
- `/tmp/ath_bench/family_eval/eval_<label>.tsv` — per-pipeline cached BLAST output

---

## 9. Bottom line

**Refiner_mdl is not the path to closing mdl-repeat's open recall gaps**. The empirical answer to "should mdl-repeat be paired with Refiner_mdl?" on Arabidopsis TAIR10 is a clear **no** under the family-level-recall metric. The contract optimizations (R5/R6/R7/R8) were correctly implemented but cannot compensate for the te_structure_filter mismatch. Future work on mdl-repeat's open items (chr4 90×80, big-genome benchmarks) should not depend on Refiner_mdl as a downstream.

The R8 header change is locked into mdl-repeat output as a low-cost metadata enhancement; it remains useful if some future downstream tool wants per-family `divergence` and `topology` without re-parsing the stats TSV.
