# mdl-repeat: Technical Documentation

## Design Principles and Algorithms

**Version**: March 2026
**Authors**: Implementation by Claude Code; design by Shuo Cheng

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Theoretical Foundation: MDL Framework](#2-theoretical-foundation-mdl-framework)
3. [Pipeline Overview](#3-pipeline-overview)
4. [Discovery Engine](#4-discovery-engine)
5. [K-mer Table and Position Index](#5-k-mer-table-and-position-index)
6. [Refinement Pipeline](#6-refinement-pipeline)
7. [MDL Scoring and Library Selection](#7-mdl-scoring-and-library-selection)
8. [Output Formats](#8-output-formats)
9. [Key Design Decisions](#9-key-design-decisions)
10. [Parameters Reference](#10-parameters-reference)

---

## 1. Introduction

### 1.1 Motivation

De novo repeat element discovery aims to identify repetitive sequences in a genome without prior knowledge of a repeat library. RepeatScout (Price, Jones & Pevzner, 2005) established the seed-and-extend paradigm: pick the highest-frequency l-mer, extend it into a consensus via banded pairwise alignment against all occurrences simultaneously, mask the found family, and repeat. This approach is effective but relies on hard-coded thresholds (`MINTHRESH`, `GOODLENGTH`) to decide which families to report.

mdl-repeat replaces these arbitrary thresholds with the **Minimum Description Length (MDL) principle**, an information-theoretic criterion that automatically determines how many families to report, how long each consensus should be, and the minimum copy number that justifies a family. The discovery engine faithfully replicates RepeatScout v1.0.7's N-sequence simultaneous banded DP extension, preserving its sensitivity while adding:

- **MDL-based model selection** — replaces all threshold-driven decisions
- **Iterative R convergence** — resolves the circular dependency between the number of accepted families and per-instance encoding costs
- **Refinement pipeline** — merge, split, fragment assembly, and prune stages to improve library quality
- **Multi-threading** — parallel position index construction and alignment refinement

### 1.2 System Requirements

- **Language**: C (C11 standard)
- **Dependencies**: Only `libm` (math library) and pthreads
- **Compiler**: GCC with `-O3 -Wall -Wextra`

---

## 2. Theoretical Foundation: MDL Framework

### 2.1 Core Principle

Given a genome $G$ of length $N$, a repeat library $L$ is a set of $R$ consensus sequences. The MDL principle states that the best library minimizes the total description length:

$$DL(G) = DL(L) + DL(G \mid L)$$

where $DL(L)$ is the cost to describe the library itself, and $DL(G \mid L)$ is the cost to describe the genome given the library. This is a two-part code in the sense of Rissanen (1978). No external thresholds are needed: the encoding scheme itself determines what constitutes a valid repeat family.

### 2.2 Rissanen's Universal Integer Code

To encode a positive integer $n$, we use Rissanen's log-star code:

$$L_{\text{int}}(n) = \log_2(c_0) + \log_2(n) + \log_2(\log_2(n)) + \log_2(\log_2(\log_2(n))) + \cdots$$

Only non-negative terms are summed. The normalizing constant $c_0 \approx 2.865064$ ensures the code satisfies the Kraft inequality ($\log_2(c_0) \approx 1.518$).

Reference values:

| $n$ | $L_{\text{int}}(n)$ (bits) |
|-----|---------------------------|
| 1   | 1.518                     |
| 2   | 2.518                     |
| 3   | 3.768                     |
| 10  | 7.364                     |
| 100 | 12.344                    |

Implementation iterates `x = log2(x)` while `x > 0`, accumulating the sum, then adds `log2(c0)`.

### 2.3 Library Encoding — DL(L)

$$DL(L) = L_{\text{int}}(R) + \sum_{r=1}^{R} \left[ L_{\text{int}}(\text{len}_r) + 2 \cdot \text{len}_r \right]$$

Each family's consensus is encoded as: its length (via $L_{\text{int}}$) plus 2 bits per base (uniform base distribution, intentionally conservative).

### 2.4 Per-Instance Encoding

A repeat instance with aligned length $a_i$ and $m_i$ edit operations is encoded as:

$$C_{\text{instance}} = L_{\text{int}}(a_i) + L_{\text{int}}(m_i + 1) + m_i \cdot \log_2(3) + \text{position\_encoding}(a_i, m_i)$$

The terms are:
- **$L_{\text{int}}(a_i)$**: length of the aligned region
- **$L_{\text{int}}(m_i + 1)$**: number of edits (the +1 avoids $L_{\text{int}}(0)$)
- **$m_i \cdot \log_2(3)$**: edit type specification (3 alternatives per edited position: substitution to one of 3 other bases)
- **position encoding**: specifies *where* the edits occur

The position encoding has three modes (selectable via `-mdl-mode`):
- **`exact`** (default): $\log_2 \binom{a_i}{m_i}$ via lgamma — exact binomial coefficient encoding of edit positions
- **`upper`**: $m_i \cdot \log_2(a_i)$ — upper bound (each edit position independently encoded)
- **`none`**: omitted — backward compatibility mode

The `exact` mode computes the binomial coefficient using the log-gamma function:

$$\log_2 \binom{a_i}{m_i} = \frac{\ln \Gamma(a_i+1) - \ln \Gamma(m_i+1) - \ln \Gamma(a_i - m_i + 1)}{\ln 2}$$

### 2.5 MDL Selection Criterion

For each candidate family with consensus length $\text{len}_r$:

$$\text{model\_cost} = L_{\text{int}}(\text{len}_r) + 2 \cdot \text{len}_r$$

$$\text{total\_savings} = \sum_{i=1}^{n} \left( 2 \cdot a_i - C_{\text{instance},i} \right)$$

$$\text{mdl\_score} = \text{total\_savings} - \text{model\_cost}$$

A family is accepted if and only if $\text{mdl\_score} > 0$. This means the compression gained by encoding its instances as references to the consensus exceeds the cost of storing the consensus itself. No minimum copy count or minimum length parameter is needed — the encoding scheme inherently penalizes short/low-copy families and rewards long/high-copy families.

### 2.6 Why MDL Works

The criterion creates a natural tension:

| Scenario | Consensus cost | Savings | MDL verdict |
|----------|---------------|---------|-------------|
| Alu-like (300 bp, 1M copies) | 600 bits | ~600M bits | Accept |
| Long TE (5 kb, 3 copies) | 10,000 bits | ~30,000 bits | Accept |
| Marginal (50 bp, 2 copies) | 100 bits | ~200 bits | Borderline |
| Noise (200 bp, 2 copies, 15% div) | 400 bits | ~250 bits | Reject |

---

## 3. Pipeline Overview

```
Input: genome.fa

Step 1: Load genome
    → Genome struct with padded sequence, boundary array, sequence IDs

Step 2: Discover consensus families (seed-and-extend)
    → CandidateList with consensus sequences and instances
    → Optionally reads/writes l-mer frequency table

Step 3: Build k-mer table + position index
    → KmerTable for refinement pipeline lookups

Step 4: Compact
    → Remove families with < 2 instances or short consensus

Step 5: Merge redundant families
    → 80-80-80 rule, k-mer Jaccard screening, union-find transitive

Step 6: Split bimodal families
    → Otsu's method on divergence, MDL-validated

Step 6b: Fragment assembly
    → Spatial co-occurrence sweep-line, nesting guard, MDL-validated

Step 7: MDL scoring and library selection
    → Iterative R convergence, greedy acceptance

Step 8: Prune marginal families
    → Exclusive coverage vs. model cost

Step 8b: Post-prune recovery pass
    → Re-score rejected families with reduced R

Step 9: Output
    → FASTA library, optional BED instances, optional TSV statistics
```

---

## 4. Discovery Engine

The discovery engine (`discover.c`, ~2200 lines) faithfully replicates RepeatScout v1.0.7's algorithm with two bug fixes:

1. **Sequence boundary check**: RepeatScout hard-codes `bStart=0; bEnd=50000000`, preventing correct multi-sequence handling. mdl-repeat computes boundaries from the actual FASTA input.
2. **Dynamic `value[]` allocation**: RepeatScout's `struct repeatllist.value[20]` overflows when l > 19. mdl-repeat allocates dynamically based on the actual l-mer length.

### 4.1 Genome Doubling

Before discovery, the genome is doubled to capture both strands:

```
[PADLENGTH N's] [forward genome] [PADLENGTH N's] [reverse complement genome]
```

- `PADLENGTH = 11000` (must be >= max extension distance `L` + slack)
- Forward copy: positions `[PADLENGTH, orig_length + PADLENGTH)`
- RC copy: positions `[orig_length + 2*PADLENGTH, 2*orig_length + 2*PADLENGTH)`
- Total length: `2 * orig_length + 2 * PADLENGTH`

This doubling is critical: it approximately doubles the number of seed occurrences, providing significantly more evidence for the simultaneous DP extension.

### 4.2 L-mer Frequency Table

L-mers are hashed symmetrically (forward and reverse complement map to the same bucket). The hash function:

```c
hash_forward  = sum(4 * ans + (base % 4)) % HASH_SIZE    // left to right
hash_revcomp  = sum(4 * ans + (3 - base)) % HASH_SIZE    // right to left
stored_hash   = max(hash_forward, hash_revcomp)
```

`HASH_SIZE = 16,000,057` (prime).

**Internal counting** scans only the forward copy of the doubled genome with per-strand TANDEMDIST filtering: an l-mer occurrence on the forward strand is only counted if its distance from the previous same-strand forward occurrence exceeds TANDEMDIST (default: 500 bp). The same applies independently for the reverse strand. This prevents inflated counts from tandem repeats.

The frequency table can be:
- Computed internally (`build_headptr_internal`)
- Read from a pre-computed file (`-freq <file>`, compatible with RepeatScout's `build_lmer_table` format)
- Written for reuse (`-freq-output <file>`)

### 4.3 Seed Selection

Seeds are selected greedily by frequency. The function `find_besttmp` performs a linear scan of the hash table for the l-mer with the highest frequency that has not been masked. A locality optimization remembers the last-found frequency and hash position to accelerate repeated scans.

### 4.4 N-Sequence Simultaneous Banded DP Extension

This is the core algorithmic contribution of RepeatScout. Given a seed l-mer with $N$ occurrences in the genome, **all $N$ sequences are extended simultaneously** using banded dynamic programming.

**DP state**: `score[2][MAXN][2*MAXOFFSET+1]`
- Dimension 1: current/previous row (rolling)
- Dimension 2: one entry per occurrence
- Dimension 3: band offsets from -MAXOFFSET to +MAXOFFSET

**Extension procedure**:

1. **Right extension**: For each consensus position (starting from the seed's right end), compute the best alignment score across all occurrences at each band offset. The consensus base at each position is chosen by majority vote weighted by alignment scores:
   - For each base (A, C, G, T), sum the scores of occurrences that have that base at the current position
   - The base with the highest total score becomes the consensus base
   - A position is "improving" if the total best score increases by at least MINIMPROVEMENT

2. **Left extension**: After right extension completes, restart from the seed's left end and extend leftward using the same algorithm. The left extension is asymmetric: it does not re-save offset checkpoints (intentional design choice from RepeatScout).

3. **Stopping criteria**:
   - Stop if no improvement for WHEN_TO_STOP consecutive positions (default: 100)
   - Stop if the total best score drops below CAPPENALTY (default: -20) relative to the current best

4. **Shannon entropy filter**: After extension, compute the Shannon entropy of the consensus. Reject low-complexity sequences where entropy exceeds MAXENTROPY (default: -0.70, using natural logarithm, yielding a negative value).

### 4.5 Masking

After a family is discovered, its occurrences in the genome are masked (set to DNA_N) to prevent re-discovery. Masking uses a separate 1-vs-1 banded DP alignment (not the simultaneous N-sequence DP) to determine the precise extent of each occurrence in the genome.

### 4.6 Instance Collection

Instances are collected from the occurrence positions that participated in the extension. Each instance records:
- Genome position (in padded coordinates)
- Aligned length and divergence
- Consensus coordinates (cons_start, cons_end)
- Edit count and alignment score
- Strand (+1 forward, -1 reverse)

A critical filter ensures that only instances mapping to the **forward copy** of the doubled genome are reported (position < orig_length + PADLENGTH), since forward and RC copies represent the same genomic locus.

---

## 5. K-mer Table and Position Index

After discovery, a separate k-mer table is built for the refinement pipeline. This uses a different data structure from the discovery hash table, optimized for position-based lookups.

### 5.1 Canonical K-mers

K-mers are packed into 64-bit integers (2 bits per base, maximum k = 31):

```c
canonical(kmer) = min(pack(kmer), revcomp(pack(kmer)))
```

The hash function uses Fibonacci hashing for excellent distribution:

```c
hash = (kmer * 11400714819323198485ULL) % table_size
```

Table size is scaled to `max(16000057, genome_length / 4)` to minimize chain lengths.

### 5.2 Position Index

After counting, `kmer_build_positions` constructs a position array for each k-mer. This enables O(1) lookup of all genome positions where a given k-mer occurs, replacing O(N) genome scans.

- Position arrays are sized at `frequency * 6` (with headroom), capped at `KMER_MAX_POSITIONS = 50000`
- Parallel construction: thread-safe atomic fetch-add for slot claiming
- Sign encoding: positive positions for forward strand, negative for reverse complement
- TANDEMDIST filtering is applied during counting (same as discovery)

### 5.3 Trimming

`kmer_trim` removes k-mers with frequency below MINTHRESH (default: 2) from the table. This reduces memory usage and speeds up lookups by eliminating singleton k-mers that cannot contribute to repeat discovery.

---

## 6. Refinement Pipeline

The refinement pipeline (`refine.c`, ~1800 lines) improves the raw discovery output through four stages.

### 6.1 Merge: 80-80-80 Rule

The merge stage detects and combines redundant families that represent the same biological repeat element.

**Stage 1: K-mer profile screening**
- For each pair of families, compute the 8-mer Jaccard index (bitset intersection/union)
- Pairs with Jaccard < 0.15 are skipped (fast rejection)

**Stage 2: Semi-global alignment**
- Align the shorter consensus against the longer one (free leading/trailing gaps on the shorter sequence)
- Both forward and reverse complement orientations are tested
- DP scoring: match +2, mismatch -3, gap -2
- For very long families where the DP matrix exceeds 10M cells, fall back to Jaccard + instance overlap

**Merge criteria (80-80-80 rule)**:
- Identity >= 80%
- Coverage >= 80% (of the shorter consensus)
- Aligned region >= 80 bp

**Relaxed criteria** (with additional instance overlap verification):
- Identity >= 70%
- Coverage >= 70%
- Plus >= 50% of the shorter family's instances overlap with the longer family's instances

**Stage 3: Transitive merges via union-find**
- If A merges with B and B merges with C, all three form one family
- The representative family is the one with the most instances

**Stage 4: Post-merge refinement**
- Instances from merged families are transferred to the representative
- The representative's consensus is re-refined using multi-k-mer seeded alignment

### 6.2 Split: Otsu's Method

The split stage detects families with bimodal divergence distributions, indicating that two distinct subfamilies have been conflated.

**Algorithm**:
1. Collect divergence values for all instances of a family
2. Build a 100-bin histogram over the divergence range [0, 0.5]
3. Apply Otsu's method to find the threshold that maximizes inter-class variance
4. Accept the split only if:
   - Bimodality score >= 0.20
   - Both clusters have >= 3 instances
   - The divergence gap at the threshold >= 0.05
   - Both sub-families have positive MDL scores (validation)

If accepted, the family is replaced by two sub-families, each with re-aligned instances and rebuilt consensus sequences.

### 6.3 Fragment Assembly

Transposable elements are often discovered as multiple short fragments when the seed-and-extend process finds different parts of the same long element independently. Fragment assembly detects and joins such fragments.

**Pre-screening: Spatial co-occurrence**

Unlike consensus k-mer Jaccard (which fails for non-overlapping fragments that share zero k-mers), this stage uses **spatial co-occurrence** in the genome:

1. Build a sorted array of `(position, end_position, family_id, instance_idx, strand)` for all instances
2. Sweep-line: for each pair of instances from different families within distance $D$, increment a co-occurrence counter
3. $D = \text{median\_consensus\_length} \times 2$, clamped to [500, 5000] bp
4. Families with > 10K instances are subsampled to 10K for efficiency

**Relationship classification**:

For each co-occurring pair (A, B):
- **NESTED** (do NOT merge): >= 50% of B's instances are spatially contained within A's instances, OR the size ratio $\text{len}_B / \text{len}_A < 0.10$
- **ADJACENT** (merge candidate): >= 3 co-occurrences with consistent orientation (>= 80% same direction) and acceptable gap statistics

**Gap analysis**:
- Compute gap distribution (distances between co-occurring instance pairs)
- Require median absolute deviation (MAD) <= 100 bp (consistent gap size)
- Negative gaps up to -20 bp are allowed (representing target site duplication overlap)

**Assembly**:
1. Concatenate consensus sequences: `A_consensus + gap_fill + B_consensus`
2. Gap fill: N-padding for positive gaps, overlap resolution for negative gaps
3. Merge instances with adjusted consensus coordinates (B family instances offset by A's length + gap)
4. **MDL validation**: accept only if $\text{mdl\_assembled} > \text{mdl}_A + \text{mdl}_B$

Transitive assembly uses union-find with a chain length cap of 10 to prevent runaway merging.

### 6.4 Prune: Exclusive Coverage

The prune stage removes marginal families whose unique genome contribution does not justify their model cost.

**Algorithm** (weakest-first iteration):
1. Sort families by MDL score (ascending — weakest first)
2. For each family, compute the exclusive coverage: genome bases covered only by this family and no other accepted family
3. If $\text{exclusive\_coverage\_savings} < \text{model\_cost}$, reject the family
4. Iterate until no more families are pruned

---

## 7. MDL Scoring and Library Selection

### 7.1 Iterative R Convergence

The number of accepted families $R$ and per-instance encoding costs are circularly dependent: some cost terms (e.g., library-level overhead $L_{\text{int}}(R)$) depend on $R$, but $R$ depends on which families have positive MDL scores. We resolve this via iteration:

1. Initialize $R = \text{num\_families}$ (upper bound)
2. Score all families with current $R$
3. Count accepted: families with score > 0 and >= 2 instances
4. If accepted count == $R$, converge; otherwise update $R$ and re-score
5. Maximum 3 iterations (empirically sufficient for convergence)

The floor $R \geq 2$ prevents degenerate cases.

### 7.2 Library-Level Cost

The total library cost includes an $L_{\text{int}}(R)$ term for encoding the number of families:

$$DL_{\text{library}} = L_{\text{int}}(R) + \sum_{r \in \text{accepted}} \text{model\_cost}_r$$

### 7.3 Greedy Selection

After convergence, families are sorted by MDL score (descending) and accepted in order while the score remains positive. The acceptance criterion is:

$$\text{mdl\_score} > 0 \quad \text{AND} \quad \text{num\_instances} \geq 2$$

### 7.4 Coverage Reporting

Unique genome coverage is computed via a bitmap (1 bit per genome position). For each accepted family's instances, the corresponding genome positions are marked. The total number of marked positions gives the bases covered. This is used for reporting only, not for selection decisions.

### 7.5 Compression Ratio

$$\text{compression\_ratio} = \frac{DL_{\text{total}}}{2 \cdot N}$$

where $DL_{\text{total}} = DL_{\text{library}} + DL_{\text{genome} \mid \text{library}}$ and $2N$ is the literal cost of encoding the genome without a library.

### 7.6 Post-Prune Recovery Pass

After pruning reduces $R$, the per-instance overhead decreases (e.g., $\lceil \log_2(R) \rceil$ becomes smaller). This may cause previously rejected families to become viable. A single recovery pass re-scores all rejected families with the final $R$; any that now have positive scores are re-accepted.

---

## 8. Output Formats

### 8.1 FASTA Library

The primary output. Only families with `mdl_score > 0` are written.

```
>R=42 length=312 copies=15 mdl=1204.3
ACGTACGT...
```

Header fields:
- `R=<id>`: family identifier
- `length=<len>`: consensus length
- `copies=<n>`: number of instances
- `mdl=<score>`: MDL score (bits saved)

### 8.2 BED Instances (`-instances`)

One line per instance of each accepted family:

```
chr1    10000    10312    R=42    850    +
```

- Column 4: family identifier
- Column 5: `int(1000 * (1 - divergence))` (0-1000 scale)
- Column 6: strand

### 8.3 TSV Statistics (`-stats`)

Per-family statistics for all accepted families:

```
family_id    consensus_length    num_instances    divergence_mean    mdl_score    model_cost    topology
```

---

## 9. Key Design Decisions

### 9.1 Genome Doubling vs. Separate Strand Processing

mdl-repeat uses genome doubling (forward + padding + reverse complement), matching RepeatScout's `add_rc()`. This approximately doubles the number of occurrences available for simultaneous DP extension, providing stronger statistical evidence for each consensus position. The alternative (processing strands separately) was found to produce significantly lower sensitivity.

### 9.2 Sensitivity-Preserving MDL Formula

The DESIGN_DOC specifies a complete per-instance encoding that includes type bit (1 bit), family identifier ($\lceil \log_2(R) \rceil$ bits), strand (1 bit), and consensus start pointer ($\lceil \log_2(\text{len}_r) \rceil$ bits) — totaling approximately 19-24 additional bits per instance.

**This was implemented and tested, then reverted.** Testing on real genomes (Arabidopsis, Rice) showed that the complete formula caused unacceptable sensitivity loss (10M+ masked bases lost on Arabidopsis). The overhead per instance was too aggressive: it caused genuine 2-3 copy families to be rejected.

The current implementation uses the original simpler formula:

$$C_{\text{instance}} = L_{\text{int}}(a_i) + L_{\text{int}}(m_i + 1) + m_i \cdot \log_2(3) + \text{position\_encoding}$$

The `mdl_instance_cost_full()` function accepts `consensus_length` and `num_families` parameters for API compatibility but does not use them in the cost computation. This preserves the option to re-introduce overhead terms in future versions with better calibration.

**Lesson**: Per-instance overhead in the MDL formula must be validated on real genomes before committing. Theoretical completeness does not guarantee practical utility.

### 9.3 Spatial Co-occurrence for Fragment Assembly

The original plan proposed k-mer Jaccard for detecting fragment pairs. This is fundamentally wrong for fragment assembly: non-overlapping fragments from different parts of the same TE share zero k-mers. The correct approach is **spatial co-occurrence** — counting how often instances from two families appear near each other in the genome. This was identified during plan review and implemented from the start.

### 9.4 No align_refine During Discovery

Calling `align_refine_family` (which re-collects instances via multi-k-mer seeding and rebuilds the consensus) on families produced by seed-and-extend was found to be **harmful**. The multi-k-mer seeding approach truncates the consensus because it relies on exact k-mer matches to anchor alignments, while the seed-and-extend consensus was built from approximate simultaneous alignment. Coverage dropped from ~48% to ~37% when this was applied. The align-refine functions are used only in the refinement pipeline (merge, split) where they are appropriate.

### 9.5 MINTHRESH=2 vs. RepeatScout's MINTHRESH=3

mdl-repeat defaults to MINTHRESH=2 (minimum l-mer frequency to seed). Since MDL replaces the arbitrary threshold as the final arbiter, allowing seeds with frequency 2 lets MDL evaluate borderline families rather than discarding them a priori. This improves sensitivity for low-copy families that genuinely compress the genome.

### 9.6 Canonical K-mers

Both the discovery hash table and the refinement k-mer table use symmetric/canonical representations:
- Discovery: `max(hash_forward, hash_revcomp)` ensures forward and RC l-mers collide
- Refinement: `min(packed_kmer, revcomp(packed_kmer))` as the canonical form

This halves the table size and ensures that a repeat family is discovered regardless of which strand the seed l-mer comes from.

---

## 10. Parameters Reference

### 10.1 Required Arguments

| Argument | Description |
|----------|-------------|
| `-sequence <file>` | Input FASTA genome file |
| `-output <file>` | Output repeat library (FASTA) |

### 10.2 Discovery Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `-freq <file>` | — | Pre-computed l-mer frequency table |
| `-freq-output <file>` | — | Write l-mer frequency table for reuse |
| `-l <int>` | auto | L-mer length (auto: $\lceil 1 + \log_4(N) \rceil$) |
| `-L <int>` | 10000 | Max extension distance per side (bp) |
| `-minthresh <int>` | 2 | Min l-mer frequency to seed |
| `-goodlength <int>` | 30 | Min consensus length pre-filter |
| `-maxgap <int>` | 5 | Max DP band offset |
| `-match <int>` | 1 | DP match score |
| `-mismatch <int>` | -1 | DP mismatch score |
| `-gap <int>` | -5 | DP gap penalty |
| `-cappenalty <int>` | -20 | Cap on penalty for exiting alignment |
| `-minimprovement <int>` | 3 | Min total score improvement per step |
| `-stopafter <int>` | 100 | Stop after N no-progress positions |
| `-maxentropy <float>` | -0.70 | Shannon entropy filter threshold |
| `-tandemdist <int>` | 500 | Min distance between counted same-strand l-mers |
| `-maxoccurrences <int>` | 10000 | Max occurrences per seed |
| `-maxrepeats <int>` | 100000 | Max families to discover |

### 10.3 Refinement Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `-threads <int>` | 1 | Number of threads for refinement |
| `-mdl-mode <mode>` | exact | MDL position encoding: `none`, `exact`, or `upper` |

### 10.4 Output Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `-instances <file>` | — | Output instance BED file |
| `-stats <file>` | — | Output family statistics TSV |
| `-v` / `-vv` | off | Verbosity level (1 or 2) |

### 10.5 Internal Constants

| Constant | Value | Location | Description |
|----------|-------|----------|-------------|
| `PADLENGTH` | 11000 | types.h | Genome padding (must be >= L) |
| `HASH_SIZE` | 16000057 | discover.c | Discovery hash table size (prime) |
| `KMER_MAX_POSITIONS` | 50000 | kmer.h | Max positions stored per k-mer |
| `ALIGN_MAXOFFSET` | 12 | align.h | Refinement DP band width (±12) |
| `REFINE_MIN_JACCARD` | 0.15 | refine.h | K-mer screening threshold for merge |
| `REFINE_MIN_IDENTITY` | 0.80 | refine.h | Merge identity threshold |
| `REFINE_MIN_COVERAGE` | 0.80 | refine.h | Merge coverage threshold |
| `REFINE_MIN_ALIGNED` | 80 | refine.h | Merge min aligned bases |
| `REFINE_BIMODALITY_THRESH` | 0.20 | refine.h | Split bimodality threshold |
| `REFINE_MAX_DP_CELLS` | 10M | refine.h | Max DP matrix size for alignment |

---

## References

1. Price, A.L., Jones, N.C., & Pevzner, P.A. (2005). De novo identification of repeat families in large genomes. *Bioinformatics*, 21(suppl_1), i351-i358.
2. Rissanen, J. (1978). Modeling by shortest data description. *Automatica*, 14(5), 465-471.
3. Otsu, N. (1979). A threshold selection method from gray-level histograms. *IEEE Transactions on Systems, Man, and Cybernetics*, 9(1), 62-66.
4. Grumbach, S., & Tahi, F. (1994). A new challenge for compression algorithms: genetic sequences. *Information Processing & Management*, 30(6), 875-886.
