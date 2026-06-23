# v7-minimal Empirical Result on TAIR10

## Filter logic tested

- HelitronScanner produces 578 filtered candidates
- v7-minimal filter: drop candidates whose consensus matches mdl-repeat library at ≥80% identity × ≥80% coverage
- 504 covered + 74 novel; v7-minimal library = mdl + 74 novel

## Result (bp-level via RepeatMasker)

```
                              recall    precision    F1
mdl-only (baseline):          0.6357    0.7469     0.6868
mdl + ALL 578 Heli (raw):     0.6800    0.7000     0.6899  (+0.30 F1)
v7-minimal (mdl + 74 novel):  0.6394    0.7409     0.6864  (-0.04 F1)  ← FILTER OVER-AGGRESSIVE
mdl + UNION (4 tools):        0.6867    0.7015     0.6940  (+0.72 F1)
```

## Why v7-minimal underperforms

Instance-level analysis of 578 HelitronScanner candidates:
- Total Heli instance bp: 3.60 Mb
- Heli bp NOT in mdl-mask (genuinely "novel" by genome position): 2.09 Mb
- Of which intersects truth TE: 0.46 Mb (22% of novel bp is real)
- Of which is FALSE POSITIVE: 1.63 Mb (78% of novel bp is non-TE)

→ HelitronScanner has ~78% FP rate at the bp level for novel territory
→ Consensus-level filter (R=80%×80%) drops 504 "covered" candidates that STILL contributed truth-bp via their distinct boundaries (instance-level novelty ≠ consensus-level novelty)
→ The 74 "novel" candidates are mostly FP, not real Helitrons

## Bottom line

v7-minimal with HelitronScanner is a NEGATIVE result:
- F1 gain: -0.04pp (worse than mdl-only)
- The structural-tool architecture (v7 / v7-minimal) is sensitive to FP rate
- HelitronScanner 2014 has ~30% FP per literature; HELIANO 2024 has ~17% FP
- Without being able to install HELIANO (conda env creation + git clone of external code both denied), v7-minimal cannot be properly validated with a modern Helitron tool

## Path forward

Drop the v7-minimal HelitronScanner path. Three remaining options:

1. **Install HELIANO 2024** — requires user permission for either conda env creation OR git clone of external code
2. **chr4 90×80 internal fix in mdl-repeat** — per FINAL_REPORT.md §7 #1 priority; algorithmic work in refine.c (conditional banded-DP re-extension for high-divergence high-copy families)
3. **Run mdl-repeat on bigger genomes** (maize 2.3Gb / wheat 17Gb) — validate scaling claim from FINAL_REPORT.md §6

Recommendation: option 2. mdl-repeat-internal improvements directly target the metric (90×80 = 0.66 → potentially 0.72) without dependency on external tool installations.
