# mdl-repeat Output Quality Enhancement Proposal — v1

**Goal**: Improve mdl-repeat library output quality beyond Stage B + Option K baseline.

**Current state**:
- chr4 family-level recall: 80×50=0.823, 80×80=0.769, 90×80=0.599 (with K)
- TAIR10 nuclear: 80×50=0.811, 80×80=0.767, 90×80=0.601 (with K)
- 24 in-scope (≥3 copy) families on chr4 still missed at 80×50
- ~50 missed at 80×80

**User constraint**: Quality first, runtime/implementation cost is acceptable.

**Design philosophy** (from CLAUDE.md / memory):
- Target: high-copy ACTIVE TEs (3+ copies, low divergence)
- OUT OF SCOPE: 1-2 copy fragments (RECON), tandem repeats (TRF), low-complexity
- Architecture: RepeatScout-derived seed-and-extend + MDL selection (NOT cdBG)

---

## Tier 1: Diagnostic-driven (highest leverage, evidence-based)

### A. Analyze 24 missed in-scope families on chr4

Sample failures from strict_recall.py:
```
size=17  chr4:10605713-10606425  (713 bp)   ← 17 copies, 700bp - should be trivial
size=16  chr4:11091152-11091593  (441 bp)
size=15  chr4:11049144-11049614  (470 bp)
size=15  chr4:11048649-11049007  (358 bp)
```

Cluster around chr4:10.5-11.1 Mb (pericentric heterochromatin). 358-713 bp, 15-17 copies. Squarely in design scope but missed.

**Action**: extract truth interval sequences for the 24 missed families, characterize them (sequence composition, GC, simple-repeat content, similarity to library entries), determine WHY each fails. Likely reveals a single common cause that targets a high-leverage fix.

**Expected effort**: 1-2 hours
**Expected gain**: unknown, depends on diagnosis. If common cause found: +5-10pp possible.

### B. Fix nested_containment_fraction over-merge bug

Documented in STAGE_A_REPORT §C and REFINE_TRACE_REPORT. Genome C (synthetic) shows nested LINE+LTR composite (2897bp) + plain LINE (2194bp) collapsed into one under-populated family by `refine.c::nested_containment_fraction`. Pre-existing bug, untouched by Stage A/B/K.

**Action**: tighten the gate — require longer/shorter to NOT be a strict prefix/suffix relationship; require at least M instances of the inner element have non-overlapping outer-element neighbors.

**Expected effort**: 0.5 day
**Expected gain**: +1-3pp on chr4 (small, since few obviously-nested families); larger on big plant genomes with many LTR-LINE composites.

### C. Improve subfamily splitting (Otsu currently fires 0 times on chr4)

`refine_split_families` uses Otsu's method on per-family divergence histogram. Should detect bimodal families (e.g., LTR with two divergence clusters → 2 subfamilies). Currently fires 0 splits on chr4 — either no bimodal families exist (unlikely) or thresholds are too conservative.

**Action**: 
1. Instrument `refine_split_families` to log how many families were CONSIDERED for splitting and why each was rejected
2. Based on log, relax thresholds (REFINE_BIMODALITY_THRESH 0.20 → 0.10? REFINE_MIN_DIV_GAP 0.05 → 0.03?) or replace Otsu with k-means / GMM
3. Validate that split families don't get re-merged in subsequent merge stage

**Expected effort**: 1-2 days
**Expected gain**: +2-5pp on 80×80 / 90×80 (better consensus per subfamily → matches truth more strictly)

---

## Tier 2: Architectural improvements

### D. Multi-pass discovery (Stage A unblocked)

HANDOFF §4 reports failed multi-pass attempt: "Pass 2 produced 0 families. Combined library = Pass 1." Tested under single-thread hash bottleneck — Pass 2 was killed at 16 min. Now with Stage A's parallel hash, Pass 2 finishes in reasonable time.

Pass 1: standard discovery → Library 1, mask coverage of accepted families
Pass 2: discovery on residual genome (already-masked positions excluded from seeding) → Library 2
Combined library: Library 1 ∪ Library 2 (deduplicated via merge-stage 80-80-80)

**Hypothesis**: Pass 1 exhausts MAXR slots on dominant high-copy elements. Pass 2 picks up rare/divergent families that Pass 1's seed competition starved. With Stage B's fallback admit gate, Pass 2 families can survive MDL even with weak standalone savings.

**Expected effort**: 1 day (wire up + test)
**Expected gain**: +0-5pp (uncertain — HANDOFF's earlier negative result is informative but inconclusive)

### E. Repeat-discriminative seed scoring

Replace pure-frequency seed selection with `freq × repeat_likelihood`:
- repeat_likelihood = local k-mer profile similarity to known-repeat regions
- Bootstrap: Pass 1 finds repeats by frequency → use those regions to define "repeat-like profile" → Pass 2 uses profile-aware seeding

Or: minimizer/syncmer-based seeding (more uniform spatial sampling).

**Expected effort**: 4-7 days
**Expected gain**: uncertain, potentially significant on dense-TE genomes

### F. BLAST-based instance recruitment

Replace `align_collect_instances`'s k-mer profile + banded DP with BLAST (dc-megablast). BLAST's affine gap model and HSP chaining handle highly diverged instances better.

**Expected effort**: 2-3 days
**Expected gain**: +1-3pp on 90×80 (stricter criteria benefit from more diverged instances being properly assigned to families)

---

## Tier 3: Parameter tuning (small, additive)

- **G. Loosen merge length-ratio guard**: 0.7 → 0.5. Cost: possible spurious cross-family merges.
- **H. Loosen Stage B fallback**: num_instances ≥ 3 → ≥ 2. Cost: more spurious low-copy families admitted.
- **I. Loosen align max-divergence**: 0.30 → 0.40. Cost: more diverged instances captured but also more noise.
- **J. Increase ALIGN_MAX_EXTENSION**: 10000 → 20000. Cost: slower refine, more memory.

Expected each: +0.5-2pp. Combined: +2-5pp. Need to test independently to avoid over-tuning.

---

## Tier 4: Out of scope (not recommended)

- RECON-style self-alignment for fragment recovery — conflicts with design philosophy
- TRF integration for tandem repeats — same
- ML-based seed scorer — overkill, unproven on TE finding

---

## Recommended execution order

1. **A first** (1-2 hours): may directly reveal the highest-leverage fix
2. **B next** (0.5 day): known bug, bounded fix
3. **D third** (1 day): leverage Stage A's parallel hash
4. **C** (1-2 days): more involved
5. **F** (2-3 days)
6. **E** (4-7 days, only if everything else falls short)
7. **G/H/I/J** as fine-tuning only after major fixes are landed

A's diagnosis result will likely re-prioritize this list.
