# mdl-repeat: Complete Technical Reference for Paper Writing

**Purpose**: This document provides all algorithmic details, mathematical formulations, data structures, thresholds, and complexity analysis needed to write a methods paper for mdl-repeat. Every formula and constant is traceable to exact source code locations.

---

## Table of Contents

1. [Problem Formulation and MDL Framework](#1-problem-formulation-and-mdl-framework)
2. [Discovery Engine: N-Sequence Simultaneous Banded DP](#2-discovery-engine)
3. [K-mer Table and Position Index](#3-k-mer-table-and-position-index)
4. [Seeded Banded Alignment (Refinement)](#4-seeded-banded-alignment)
5. [Refinement Pipeline](#5-refinement-pipeline)
6. [MDL Scoring and Library Selection](#6-mdl-scoring-and-library-selection)
7. [Large Genome Support](#7-large-genome-support)
8. [Data Structures](#8-data-structures)
9. [Complete Parameter Reference](#9-complete-parameter-reference)
10. [Complexity Analysis](#10-complexity-analysis)

---

## 1. Problem Formulation and MDL Framework

### 1.1 Two-Part Code

Given a genome $G$ of length $N$, a repeat library $L$ is a set of $R$ consensus sequences. The Minimum Description Length (MDL) principle selects the library that minimizes:

$$DL(G) = DL(L) + DL(G \mid L)$$

where $DL(L)$ is the cost to describe the library, and $DL(G \mid L)$ is the cost to describe the genome given the library (Rissanen, 1978). No external thresholds (minimum copy count, minimum consensus length) are needed.

### 1.2 Rissanen's Universal Integer Code

To encode a positive integer $n$:

$$L_{\text{int}}(n) = \log_2(c_0) + \log_2(n) + \log_2(\log_2(n)) + \log_2(\log_2(\log_2(n))) + \cdots$$

Only non-negative terms are summed. The normalizing constant $c_0 \approx 2.865064$ ensures the Kraft inequality holds ($\log_2(c_0) = 1.5179605508986484$).

**Implementation** (`mdl.c:15-34`): Iterative loop accumulates `log2(val)` while `val > 0.0`, then adds `LOG2_C0`.

Reference values:

| $n$ | $L_{\text{int}}(n)$ (bits) |
|-----|---------------------------|
| 1   | 1.518                     |
| 2   | 2.518                     |
| 3   | 3.768                     |
| 10  | 7.364                     |
| 100 | 12.344                    |
| 1000| 17.830                    |

### 1.3 Library Encoding Cost

$$DL(L) = L_{\text{int}}(R) + \sum_{r=1}^{R} \left[ L_{\text{int}}(\text{len}_r) + 2 \cdot \text{len}_r \right]$$

Per-family model cost:

$$\text{model\_cost}_r = L_{\text{int}}(\text{len}_r) + 2 \cdot \text{len}_r$$

The $2 \cdot \text{len}_r$ term encodes the consensus at 2 bits/base (uniform base distribution, intentionally conservative). (`mdl.c:123-126`)

### 1.4 Per-Instance Encoding Cost

A repeat instance with aligned length $a_i$ and $m_i$ edit operations:

$$C_{\text{instance}}(a_i, m_i) = L_{\text{int}}(a_i) + L_{\text{int}}(m_i + 1) + m_i \cdot \log_2 3 + P(a_i, m_i)$$

where:
- $L_{\text{int}}(a_i)$: aligned region length
- $L_{\text{int}}(m_i + 1)$: edit count (+1 avoids $L_{\text{int}}(0)$)
- $m_i \cdot \log_2 3$: edit type specification (substitution to 1 of 3 alternative bases)
- $P(a_i, m_i)$: position encoding for where edits occur

**Position encoding modes** (selectable via `-mdl-mode`):

| Mode | Formula | Description |
|------|---------|-------------|
| `exact` (default) | $\log_2 \binom{a_i}{m_i} = \frac{\ln\Gamma(a_i+1) - \ln\Gamma(m_i+1) - \ln\Gamma(a_i - m_i + 1)}{\ln 2}$ | Exact binomial coefficient via lgamma |
| `upper` | $m_i \cdot \log_2(a_i)$ | Upper bound (each edit position independently) |
| `none` | $0$ | Omitted (backward compatibility) |

(`mdl.c:66-107`)

**Clamping**: $a_i$ clamped to $[1, \infty)$; $m_i$ clamped to $[0, a_i]$. (`mdl.c:76-78`)

### 1.5 MDL Selection Criterion

Per-family:

$$\text{mdl\_score}_r = \sum_{i=1}^{n_r} \left( 2 \cdot a_i - C_{\text{instance}}(a_i, m_i) \right) - \text{model\_cost}_r$$

A family is **accepted** iff $\text{mdl\_score}_r > 0$ and $n_r \geq 2$. (`mdl.c:128-156`)

### 1.6 Design Decision: Omitted Per-Instance Overhead

A theoretically complete encoding would include per-instance overhead: type bit (1 bit), family identifier ($\lceil\log_2 R\rceil$ bits), strand bit (1 bit), and consensus start pointer ($\lceil\log_2(\text{len}_r)\rceil$ bits) — totaling ~19–24 additional bits per instance. Testing on real genomes showed 3–6% sensitivity loss. The current formula omits these terms, prioritizing sensitivity over theoretical completeness. The `mdl_instance_cost_full()` function accepts `consensus_length` and `num_families` parameters for API extensibility but does not use them. (`mdl.c:54-64`)

**Consequence**: The iterative R convergence loop (§6.1) becomes effectively inert since per-instance cost does not depend on R.

---

## 2. Discovery Engine

### 2.1 Overview

The discovery engine (`discover.c`, ~2000 lines) implements the RepeatScout (Price et al. 2005) seed-and-extend paradigm with a thread-safe `DiscoverContext` encapsulation. All mutable state lives in a heap-allocated struct, enabling parallel chunked discovery.

### 2.2 Genome Doubling

Before discovery, the genome is doubled to capture both strands:

```
[PADLENGTH N's] [forward genome] [PADLENGTH N's] [reverse complement genome]
```

- `PADLENGTH = 11000` (`types.h:22`)
- Forward copy: positions $[\text{PADLENGTH}, N + \text{PADLENGTH})$
- RC copy: positions $[N + 2 \cdot \text{PADLENGTH}, 2N + 2 \cdot \text{PADLENGTH})$
- Total length: $2N + 2 \cdot \text{PADLENGTH}$

This approximately doubles seed occurrences, providing stronger statistical evidence for simultaneous DP extension. (`discover.c:1575-1597`)

### 2.3 L-mer Frequency Counting

**L-mer length**: Auto-computed as $l = \lceil 1 + \log_4 N \rceil$, max 31. (`discover.c:1758-1759`, `main.c:36-39`)

**Hash function** (`discover.c:126-141`):

$$h_{\text{fwd}} = \left(\sum_{x=0}^{l-1} 4 \cdot h + (\text{base}_x \bmod 4)\right) \bmod H$$
$$h_{\text{rc}} = \left(\sum_{x=0}^{l-1} 4 \cdot h + (3 - \text{base}_{l-1-x}) \bmod 4\right) \bmod H$$
$$h = \max(h_{\text{fwd}}, h_{\text{rc}})$$

where $H = 16{,}000{,}057$ (prime). Symmetric hashing ensures forward and RC l-mers map to the same bucket.

**TANDEMDIST filtering** (`discover.c:332-378`): An l-mer occurrence on a given strand is only counted if its distance from the previous same-strand occurrence exceeds `TANDEMDIST` (default 500 bp). Per-strand tracking via `lastplusocc` (forward) and `lastminusocc` (reverse complement). Prevents inflated counts from tandem arrays.

**Two-pass position deduplication** (`discover.c:459-506`): After counting, a second pass removes tandem duplicates within `TANDEMDIST` on the same strand using signed position markers.

### 2.4 N-Sequence Simultaneous Banded DP Extension

Given a seed l-mer with $N$ occurrences in the genome, **all $N$ sequences are extended simultaneously** using banded dynamic programming.

**DP state**: `score[2][N][2 \cdot \text{MAXOFFSET} + 1]`
- Dimension 1: current/previous row (rolling buffer)
- Dimension 2: one entry per occurrence ($N \leq$ `MAXN` = 10000)
- Dimension 3: band offsets from $-\text{MAXOFFSET}$ to $+\text{MAXOFFSET}$

**Right extension recurrence** (`discover.c:594-655`):

For consensus position $y$, sequence $n$, offset $o$, and candidate base $a$:

$$S_y(n, o) = \max \begin{cases}
S_{y-1}(n, o+1) + \text{GAP} & \text{(gap in sequence, } o < \text{MAXOFFSET}) \\[4pt]
S_{y-1}(n, o) + \delta(a, g_n(o,y)) & \text{(diagonal match/mismatch)} \\[4pt]
\max_{o' \in [-M, o)} \left[ S_{y-1}(n, o') + (o - o') \cdot \text{GAP} + \delta'(a, g_n, o', o, y) \right] & \text{(multi-gap in consensus)}
\end{cases}$$

where:
- $\delta(a, b) = \text{MATCH}$ if $a = b$, else $\text{MISMATCH}$ (default: +1 / -1)
- $g_n(o, y)$: genome base at position computed from `pos[n]`, offset $o$, and consensus position $y$
- $\delta'$: returns MATCH if $a$ matches **any** base in the gap range $[o', o]$, else MISMATCH
- $M = \text{MAXOFFSET}$ (default 5 for discovery, `discover.c:1607`)

**Left extension** (`discover.c:662-723`): Symmetric but direction-reversed. Uses `(w+1) % 2` for row indexing. Asymmetry: left extension does NOT checkpoint `score_of_besty` or `savebestscore`.

**Strand handling**: For reverse-complement occurrences (`rev[n] = 1`), genome positions are mapped as `pos[n] - (offset + y - L - l) - 1` with complement base lookup.

**Boundary check** (`discover.c:601-608`): Before each DP cell computation, verify the genome position stays within the sequence boundary $[\text{bStart}, \text{bEnd})$ computed via `get_boundaries()`.

### 2.5 Consensus Base Selection

At each extension position $y$, for each candidate base $a \in \{A, C, G, T\}$:

$$\text{total}(a) = \sum_{n=0}^{N-1} \max\left(0, \max_{o} S_y(n, o) - \text{CAPPENALTY}\right)$$

where the per-sequence best score is floored at `bestbestscore[n] + CAPPENALTY` (default CAPPENALTY = -20).

**The consensus base is chosen as** $\arg\max_a \text{total}(a)$ — pure argmax over summed best scores, not weighted majority vote. (`discover.c:862-886`)

### 2.6 Stopping Criteria

Three interacting parameters control extension termination:

**MINIMPROVEMENT** (`discover.c:760-773`): Checkpoint updated only if:
$$\text{totalbestscore}_y \geq \text{besttotalbestscore} + (y - y^*) \cdot \text{MINIMPROVEMENT}$$
where $y^*$ is the last checkpoint position. Default MINIMPROVEMENT = 3.

**WHEN_TO_STOP** (`discover.c:889`): Extension halts if $y - y^* \geq \text{WHEN\_TO\_STOP}$. Default = 100.

**CAPPENALTY** (`discover.c:867`): Per-sequence score floor: $\max(0, \text{bestbestscore}[n] + \text{CAPPENALTY})$. Default = -20.

### 2.7 Shannon Entropy Filter

After extension, low-complexity sequences are rejected:

$$H = \sum_{b \in \{A,C,G,T\}} \frac{c_b}{l} \cdot \ln\left(\frac{c_b}{l}\right)$$

where $c_b$ is the count of base $b$ in the l-mer. **Reject if** $H > \text{MAXENTROPY}$ (default -0.70, natural log, negative value). (`discover.c:188-202`)

### 2.8 Masking

After a family is discovered, its occurrences are masked via 1-vs-1 banded DP alignment (not the simultaneous N-sequence DP). The masking DP (`compute_maskscore_right/left`, `discover.c:985-1066`) uses the same recurrence as §2.4 but with a single sequence against the consensus. Extension continues until `WHEN_TO_STOP` consecutive non-improving positions. Masked positions are set to DNA_N. (`discover.c:1073-1216`)

### 2.9 Instance Collection from Extension

For each occurrence $n$ that participated in extension (`discover.c:1501-1565`):

**Forward strand** ($\text{rev}[n] = 0$):
- $\text{genome\_start} = \text{pos}[n] + w^* - L$
- $\text{genome\_end} = \text{pos}[n] + y^* - L + 1$

**Reverse strand** ($\text{rev}[n] = 1$):
- $\text{genome\_start} = \text{pos}[n] - (y^* - L - l) - 1$
- $\text{genome\_end} = \text{pos}[n] - (w^* - L - l)$

where $y^*$ and $w^*$ are the best right and left extension positions. Edit count is computed by direct base-by-base comparison with the consensus. Only instances mapping to the forward copy of the doubled genome are reported (position < $N + \text{PADLENGTH}$).

### 2.10 DiscoverContext: Thread-Safe Encapsulation

All mutable state lives in a heap-allocated `DiscoverContext` struct (`discover.c:62-110`), containing:

| Section | Fields | Purpose |
|---------|--------|---------|
| Parameters | `l, L, MAXOFFSET, MAXN, MAXR, MATCH, MISMATCH, GAP, CAPPENALTY, MINIMPROVEMENT, WHEN_TO_STOP, MAXENTROPY, GOODLENGTH, MINTHRESH, TANDEMDIST` | Algorithm parameters |
| Genome data | `sequence, sequence_owned, removed, length, orig_length, disc_boundaries, disc_num_sequences` | Doubled genome and masking array |
| Workspace | `master, masters[], masterstart[], masterend[], pos[], rev[], upperBoundI[], N` | Extension buffers |
| DP arrays | `score[2][MAXN][2M+1], score_of_besty[MAXN][2M+1], maskscore[2][2M+1], bestbestscore[], savebestscore[]` | DP state |
| Loop state | `besty, bestw, R, totalbestscore, besttotalbestscore` | Extension progress |

Public API: `CandidateList *discover_families(const Genome *genome, const DiscoverParams *params)` — read-only inputs, no global state, multiple calls safe to run in parallel. (`discover.c:1724-1971`)

---

## 3. K-mer Table and Position Index

### 3.1 Canonical K-mers

K-mers are packed into 64-bit integers (2 bits/base, max $k = 31$):

$$\text{packed}(s) = \sum_{i=0}^{k-1} s_i \cdot 4^{k-1-i}$$

where $s_i \in \{0,1,2,3\}$ is the numeric encoding. (`kmer.c:91-100`)

Canonical form: $\text{canonical}(s) = \min(\text{packed}(s), \text{packed}(\overline{s}^R))$ where $\overline{s}^R$ is the reverse complement. (`kmer.c:112-116`)

### 3.2 Hash Table

**Hash function**: Fibonacci hashing for excellent distribution:

$$h(k) = (k \times 11400714819323198485) \bmod T$$

where $11400714819323198485 = \lfloor 2^{64} / \varphi \rfloor$ ($\varphi$ = golden ratio). (`kmer.c:120-125`)

**Table size**: $T = \text{next\_prime}(\max(16{,}000{,}057, N/4))$. (`kmer.c:232-234`)

**Collision resolution**: Chaining with linked list.

### 3.3 Parallel Counting with Striped Locks

- $\text{NUM\_STRIPES} = 4096$ mutex locks (`kmer.c:131`)
- Lock assignment: $\text{stripe} = h(k) \bmod \text{NUM\_STRIPES}$ (`kmer.c:162-163`)
- Per-thread memory pools: `POOL_BLOCK_SIZE = 4096` KmerEntry per block (`kmer.c:10`)
- TANDEMDIST filtering during counting: per-strand tracking with `last_plus_occ` / `last_minus_occ` (`kmer.c:174-185`)

### 3.4 Position Index: Two-Phase Parallel Construction

**Phase 1 — Count** (`kmer.c:480-499`): Workers scan assigned genome chunks, atomically increment position counts.

**Phase 2 — Allocate** (`kmer.c:560-572`): For each k-mer, allocate position array sized at $\min(\text{num\_positions}, \text{KMER\_MAX\_POSITIONS})$ where $\text{KMER\_MAX\_POSITIONS} = 50{,}000$ (`kmer.h:24`).

**Phase 3 — Fill** (`kmer.c:574-597`): Sequential fill in genome-coordinate order for determinism. Lowest-coordinate positions are always kept when truncating at capacity. (`kmer.c:662-665`)

**Position encoding**: Positive values for forward strand, negative for reverse complement (sign-encoded strand).

### 3.5 Trimming

`kmer_trim` removes k-mers with frequency below `MINTHRESH` (default 2). (`kmer.c:399-428`)

---

## 4. Seeded Banded Alignment (Refinement)

The refinement pipeline uses a different alignment strategy from the discovery engine: multi-k-mer seeded banded alignment anchored on exact k-mer matches.

### 4.1 Consensus K-mer Set

For a family consensus of length $L_c$, extract $L_c - k + 1$ k-mers (capped at `MAX_CONS_KMERS` = 10000). Deduplicate canonical k-mers in a hash set. (`align.c:27-93`)

### 4.2 Seed Genome Scan

For each unique consensus k-mer, look up all genome positions from the pre-built KmerTable position index. Complexity: $O(\sum_j f_j)$ where $f_j$ is the frequency of consensus k-mer $j$ — avoids $O(N)$ genome scan. (`align.c:103-149`)

Returns up to `MAX_SEED_HITS = 50000` seed hits, each recording `(genome_pos, cons_pos, strand)`.

### 4.3 Seed Hit Clustering

Hits are sorted by `(strand, genome_pos)` and clustered using distance threshold:

$$d_{\text{merge}} = \lfloor 1.5 \cdot L_c \rfloor$$

Within each cluster, the **anchor** is the hit with consensus position closest to $L_c / 2$ (midpoint). (`align.c:171-210`)

### 4.4 Banded DP Alignment

For each anchor, a banded DP alignment extends bidirectionally.

**Band width**: $W = 2 \cdot \text{MAXOFFSET} + 1$ (default MAXOFFSET = 12, so $W = 25$). (`align.c:266`)

**DP recurrence** (`align.c:297-330`): For consensus position $i$, band offset $o$:

$$D_i(o) = \max \begin{cases}
D_{i-1}(o) + \delta(\text{cons}_i, \text{genome}(i, o)) & \text{(diagonal)} \\
D_{i-1}(o+1) + g & \text{(gap in genome)} \\
D_{i-1}(o-1) + g & \text{(gap in consensus)} \\
S^* + \text{CAPPENALTY} & \text{(score floor)}
\end{cases}$$

where:
- $\delta(a,b) = +1$ if match, $-1$ if mismatch (ALIGN_MATCH / ALIGN_MISMATCH)
- $g = -5$ (g_align_gap, configurable via `-refine-gap`)
- $S^*$: best score seen so far in this direction
- CAPPENALTY = -20

**Genome position mapping** (`align.c:235-242`):
- Forward: $\text{gp} = \text{anchor\_genome} + (i - \text{anchor\_cons}) + o$
- Reverse: $\text{gp} = \text{anchor\_genome} - (i - \text{anchor\_cons}) + o$

**Stall detection**: Extension stops after `ALIGN_WHEN_TO_STOP = 100` consecutive non-improving positions. (`align.c:339-340`)

**Divergence computation** (`align.c:423-445`): Substitution-only counting (gaps not tracked):

$$\text{divergence} = \frac{\text{edits}}{\text{compared\_positions}}$$

**Rejection**: Instance rejected if divergence $> 0.30$ (`g_align_max_divergence`). (`align.c:445`)

### 4.5 Instance Deduplication

New instances overlapping >50% of a shorter existing instance are rejected:

$$\frac{\text{overlap}}{\min(\text{len}_{\text{new}}, \text{len}_{\text{existing}})} > 0.5 \implies \text{reject}$$

(`align.c:525-545`)

### 4.6 Consensus Rebuild: Score-Weighted Majority Voting

For each consensus position $p$, the base is chosen by weighted vote:

$$\text{consensus}[p] = \arg\max_{b \in \{A,C,G,T\}} \sum_{i : p \in [\text{cons\_start}_i, \text{cons\_end}_i)} w_i \cdot \mathbb{1}[\text{genome}_i(p) = b]$$

where $w_i = \max(1, \text{score}_i)$ is the alignment score. (`align.c:579-637`)

### 4.7 Consensus Extension

Bidirectional extension by majority vote on flanking context:

**Dynamic support thresholds** (`align.c:707-708`):
- Extended > 1000 bp: require $\geq 3$ supporting instances
- Extended > 5000 bp: require $\geq 5$ supporting instances
- Always require $\geq 2$ instances and $\geq 50\%$ agreement for best base

**Extension caps** based on instance count (`align.c:791-796`):
- 2–3 instances: max 500 bp
- 4–9 instances: max 3000 bp
- $\geq 10$ instances: max 10000 bp

**Source eligibility**: Only instances within `EXTENSION_SLACK = 15` bp of the consensus edge can contribute to extension. (`align.c:677-679`)

### 4.8 Iterative Refinement

Loop: `collect_instances → extend → rebuild`. Max `ALIGN_MAX_ITERATIONS = 10` iterations. Convergence: no extension AND (<1% consensus change or zero change). (`align.c:808-829`)

### 4.9 Parallelization

Work-stealing via `__atomic_fetch_add` on a shared family counter. Each thread grabs the next unprocessed family. (`align.c:835-952`)

---

## 5. Refinement Pipeline

### 5.1 Merge: 80-80-80 Rule

**Stage 1: K-mer profile screening** (`refine.c:26-57`)
- 8-mer bitset profiles: `PROFILE_BITS = 65536` ($4^8$), stored as 1024 `uint64_t` words
- Jaccard index via popcount: $J = |A \cap B| / |A \cup B|$
- Skip pairs with $J < 0.15$ (`REFINE_MIN_JACCARD`)

**Stage 2: Semi-global alignment** (`refine.c:79-183`)

DP scoring: match $+2$, mismatch $-3$, gap $-2$ (`REFINE_MATCH/MISMATCH/GAP`).

Free leading/trailing gaps on query (shorter sequence) only.

DP cell limit: $(\text{qlen}+1) \times (\text{tlen}+1) \leq 10^7$ (`g_refine_max_dp_cells`). (`refine.c:87-92`)

Both orientations tested (forward and reverse complement of target).

**Merge criteria** (`refine.c:543-569`):

| Criterion | Identity | Coverage | Aligned bp | Additional |
|-----------|----------|----------|------------|------------|
| **Strict** | $\geq 80\%$ | $\geq 80\%$ | $\geq 80$ bp | — |
| **Relaxed** | $\geq 70\%$ | $\geq 70\%$ | — | $\geq 50\%$ instance overlap |
| **DP limit fallback** | — | — | — | Jaccard $\geq 0.80$ + $\geq 50\%$ instance overlap |

**Stage 3: Transitive merges** via union-find with path compression (`refine.c:294-311`). Representative: family with most instances.

**Stage 4: Post-merge re-refinement** via `align_refine_family()` for families that absorbed others. (`refine.c:683-705`)

### 5.2 Split: Otsu's Method on Divergence Distributions

**Histogram**: 100 bins over divergence range $[d_{\min}, d_{\max}]$. (`refine.c:753-762`)

**Otsu's threshold** (`refine.c:787-802`): For each candidate threshold $t$, compute inter-class variance:

$$\sigma^2_B(t) = \frac{w_0 \cdot w_1 \cdot (\mu_0 - \mu_1)^2}{n^2}$$

where $w_0, w_1$ are class sizes, $\mu_0, \mu_1$ are class means. Select $t^* = \arg\max_t \sigma^2_B(t)$.

**Bimodality score**: $B = \sigma^2_B(t^*) / \sigma^2_{\text{total}}$. (`refine.c:804`)

**Acceptance criteria** (`refine.c:944-1061`):

| Check | Threshold | Source |
|-------|-----------|--------|
| Minimum instances | $\geq 10$ | `REFINE_MIN_SPLIT_INSTANCES` |
| Bimodality score | $B \geq 0.20$ | `REFINE_BIMODALITY_THRESH` |
| Valley depth (if $0.20 \leq B < 0.40$) | valley height $< \text{lower\_mode\_height}/2$ | `refine.c:807-825` |
| Minimum cluster size | $\geq 3$ per sub-family | `REFINE_MIN_CLUSTER_SIZE` |
| Mean divergence gap | $\bar{d}_{\text{hi}} - \bar{d}_{\text{lo}} \geq 0.05$ | `REFINE_MIN_DIV_GAP` |
| MDL validation | $\text{mdl}_{\text{lo}} + \text{mdl}_{\text{hi}} > \text{mdl}_{\text{orig}}$ | `refine.c:1053-1061` |

**Consensus rebuild**: Weighted majority voting from reassigned instances (`refine.c:835-888`).

### 5.3 Fragment Assembly: Spatial Co-occurrence

Transposable elements may be discovered as multiple fragments. Non-overlapping fragments share zero k-mers, so **spatial co-occurrence** in the genome is used instead of k-mer Jaccard.

**Distance threshold** (`refine.c:1714-1721`):

$$D = \min(\max(2 \times \text{median\_consensus\_length}, 500), 5000)$$

**Sweep-line algorithm** (`refine.c:1749-1811`):
- Build sorted array of `(start, end, family_id, instance_idx, strand)` for all instances
- Subsampling: families with >10K instances reduced to 10K
- For each pair of instances from different families within distance $D$: increment co-occurrence counter and track orientation consistency

**Pair filtering** (`refine.c:1830-1879`):

| Check | Threshold |
|-------|-----------|
| Co-occurrence count | $\geq 3$ |
| Nesting guard | Neither family $\geq 50\%$ contained in the other |
| Size ratio | $\min(\text{len}_A, \text{len}_B) / \max(\text{len}_A, \text{len}_B) \geq 0.10$ |
| Orientation consistency | $\geq 80\%$ same direction |

**Gap analysis** (`refine.c:1936-1996`):
- Compute gap distribution between co-occurring instance pairs
- Median gap and median absolute deviation (MAD)
- Reject if MAD > 100 (inconsistent gap) or median gap $\notin [-20, \min(\text{median} + 2\text{MAD}, 500)]$
- Negative gaps up to -20 bp allowed (target site duplication overlap)

**Assembly** (`refine.c:1998-2062`):
- Positive gap: `consensus_A + N-padding + consensus_B`
- Overlap ($\leq 0$): `consensus_A[0..len_A-overlap) + consensus_B`

**MDL validation**: Accept only if $\text{mdl}_{\text{assembled}} > \text{mdl}_A + \text{mdl}_B$. (`refine.c:2064-2109`)

**Chain length cap**: Union-find with maximum chain length 10 to prevent runaway merging. (`refine.c:1903-1925`)

### 5.4 Prune: Exclusive Coverage

Removes marginal families whose unique genome contribution does not justify their model cost.

**Algorithm** (`refine.c:1487-1605`):
1. Sort accepted families by MDL score ascending (weakest first)
2. Build genome coverage array: `cov[pos]` = number of accepted families covering position `pos`
3. For each family (weakest first):
   - For each instance, count exclusive bases ($\text{cov}[\text{pos}] = 1$)
   - **Instance-level filter**: skip instance if exclusive bases $< \text{aligned\_length} / 4$ (< 25%) (`refine.c:1565`)
   - Compute exclusive savings: $\sum_i (2 \cdot a_i^{\text{excl}} - C_{\text{instance}}(a_i^{\text{excl}}, m_i^{\text{excl}}))$
   - **Prune** if no instances pass the 25% filter, or if exclusive savings $< \text{model\_cost}$
   - Decrement coverage counts for pruned family

---

## 6. MDL Scoring and Library Selection

### 6.1 Iterative R Convergence

(`mdl.c:172-201`)

1. Initialize $R = \max(n_{\text{families}}, 2)$
2. Repeat (max 3 iterations):
   a. Score all families with current $R$
   b. Count accepted: families with score $> 0$ and instances $\geq 2$
   c. $R_{\text{new}} = \max(\text{accepted}, 2)$
   d. If $R_{\text{new}} = R$: converge. Else $R \leftarrow R_{\text{new}}$.

**Note**: Because per-instance cost does not currently depend on $R$ (§1.6), this loop converges in 1 iteration.

### 6.2 Greedy Selection

(`mdl.c:203-223`)

Sort families by MDL score descending. Accept in order while $\text{mdl\_score} > 0$ and $n_{\text{instances}} \geq 2$.

### 6.3 Coverage Bitmap

(`mdl.c:225-260`)

Bitmap of $\lceil N/8 \rceil$ bytes. For each accepted family's instances, mark genome positions:

```c
BIT_SET(pos): covered[pos >> 3] |= (1 << (pos & 7))
BIT_GET(pos): covered[pos >> 3] & (1 << (pos & 7))
```

Used for reporting only, not selection decisions.

### 6.4 Final Description Length

(`mdl.c:262-268`)

$$DL_{\text{library}} = L_{\text{int}}(R_{\text{final}}) + \sum_{r \in \text{accepted}} \text{model\_cost}_r$$

$$DL_{\text{total}} = 2N - \text{total\_savings} + DL_{\text{library}}$$

$$\text{compression\_ratio} = \frac{DL_{\text{total}}}{2N}$$

### 6.5 Post-Prune Recovery Pass

(`main.c:998-1045`)

After pruning reduces $R$, re-score all rejected families with the reduced $R_{\text{final}} = \max(\text{accepted\_after\_prune}, 2)$. Families that now have positive MDL scores are re-accepted.

---

## 7. Large Genome Support

### 7.1 Genome Sampling

**Activation**: genome length $> \text{sample\_size}$ (default 1 Gb). (`main.c:765`)

**Algorithm** (`main.c:181-274`):

1. **Tile partitioning**: Divide each sequence into non-overlapping tiles of size $w$ (default 1 Mb).
   $$\text{tiles\_per\_seq}[i] = \lfloor \text{seq\_len}_i / w \rfloor$$

2. **Target window count**: $n_w = \lfloor \text{sample\_size} / w \rfloor$, clamped to $[1, \text{total\_tiles}]$.

3. **Partial Fisher-Yates shuffle** (`main.c:221-227`):
   ```
   srand(seed)
   For i = 0 to n_w - 1:
       j = i + rand() % (total_tiles - i)
       swap(tiles[i], tiles[j])
   ```
   Default seed = 42.

4. **Sort selected tiles** by genomic coordinate for sequential memory access (`main.c:230`).

5. **Create sampled genome** via `genome_create_chunk()` with selected tiles as segments.

**Coordinate remapping** (`main.c:280-320`): After discovery, for each instance:
1. Compute raw position: $\text{raw} = \text{position} - \text{PADLENGTH}$
2. Binary search to find containing segment $s$: $\text{raw} \geq \text{seg\_start}[s]$
3. Remap: $\text{original} = \text{segments}[s].\text{raw\_start} + (\text{raw} - \text{seg\_start}[s])$
4. Update position: $\text{position} = \text{original} + \text{PADLENGTH}$
5. Update `seq_index` to original genome's sequence index

**Memory management**: Original genome sequence buffer freed after sampling; metadata (boundaries, sequence_ids) retained for BED output. (`main.c:783-788`)

### 7.2 Chunked Discovery

**Activation**: genome length $> \text{chunk\_size}$ (default 200 Mb). (`main.c:854`)

**Sequence segmentation** (`main.c:377-436`):

1. **Split threshold**: $\tau = 1.8 \times \text{chunk\_size}$

2. If sequence length $> \tau$:
   - Number of parts: $p = 2^k$ where $k = \lceil \log_2(\text{seq\_len} / \tau) \rceil$
   - Part base size: $\text{part\_base} = \text{seq\_len} / p$

3. **Overlap**: Adjacent segments overlap by $L$ bases (default 10000) at each boundary. Clamped to sequence limits. (`main.c:411-416`)

**LPT bin packing** (`main.c:438-487`):

1. Sort segments by length descending (Longest Processing Time first)
2. Number of bins: $b = \lceil \text{total\_seg\_size} / \text{chunk\_size} \rceil$, bounded by $[2, n_{\text{segments}}]$
3. For each segment (largest first): assign to bin with smallest current total
   $$\text{bin}^* = \arg\min_b \text{bin\_size}[b]$$

This minimizes $\max_b(\text{bin\_size}[b])$, optimal for parallel wall-clock time.

**Parallel execution** (`main.c:501-550`):
- Batch size: $\min(n_{\text{bins}}, n_{\text{threads}})$
- Each worker creates a chunk genome via `genome_create_chunk()` (copies sequence data, shares sequence_id strings)
- Each chunk independently calls `discover_families()` with its own `DiscoverContext`
- **Per-chunk l-mer length**: $l_{\text{chunk}} = \lceil 1 + \log_4(N_{\text{chunk}}) \rceil$ — shorter than full-genome $l$, increasing seed sensitivity

**Result concatenation** (`main.c:570-578`): All chunk results merged into single `CandidateList` with renumbered family IDs.

### 7.3 Composition

1. **Sampling** (if triggered): reduce 10+ Gb genome to ~1 Gb
2. **Chunked discovery** (if triggered): split ~1 Gb into ~200 Mb parallel chunks
3. **Refinement**: runs on the full (sampled) genome, not per-chunk
4. **Output**: coordinates remapped to original genome if sampling was used

Small genomes ($< \text{chunk\_size}$): no splitting, single-threaded.

---

## 8. Data Structures

### 8.1 Core Types (`types.h`)

| Type | Definition | Purpose |
|------|-----------|---------|
| `gpos_t` | `int64_t` | Genome position (supports >2 Gb genomes) |
| `glen_t` | `int64_t` | Genome/sequence length |
| `freq_t` | `int32_t` | K-mer frequency |
| `uid_t` | `int32_t` | Family ID |

DNA encoding: A=0, C=1, G=2, T=3, N=99.

Complement: $\overline{c} = 3 - c$ (for $c \neq 99$).

### 8.2 Genome (`genome.h`)

```
Genome {
    char    *sequence       // Numeric DNA encoding, padded at start
    glen_t   length         // raw_length + PADLENGTH
    glen_t   raw_length     // Actual bases (no padding)
    gpos_t  *boundaries     // Cumulative sequence boundaries (pre-padding)
    int      num_sequences  // Number of FASTA records
    char   **sequence_ids   // Sequence identifiers
}
```

**Boundary convention** for $N$ sequences:
- `boundaries[i]` = end of sequence $i$ (= start of sequence $i+1$) in raw coordinates
- `boundaries[N-1]` = raw_length + 1 (sentinel)
- `boundaries[N]` = 0 (terminator)
- Padded coordinate: add PADLENGTH to raw

### 8.3 KmerEntry (`kmer.h`)

```
KmerEntry {
    uint64_t  kmer            // Canonical k-mer, packed 2 bits/base
    freq_t    frequency       // Occurrence count
    gpos_t    last_plus_occ   // Last forward-strand position (TANDEMDIST)
    gpos_t    last_minus_occ  // Last reverse-strand position (TANDEMDIST)
    gpos_t   *positions       // Position array (sign = strand)
    int32_t   num_positions   // Count of stored positions
    int32_t   cap_positions   // Allocated capacity (max 50000)
    KmerEntry *next           // Hash chain
}
```

### 8.4 Instance (`candidates.h`)

```
Instance {
    gpos_t position          // Start in padded genome coordinates
    glen_t aligned_length    // Alignment length
    int    cons_start        // Consensus start (0-based)
    int    cons_end          // Consensus end (exclusive)
    int    num_edits         // Edit count from alignment
    float  divergence        // Divergence (0.0–1.0)
    int    score             // Alignment score
    int8_t strand            // +1 forward, -1 reverse
    int    seq_index         // Which FASTA sequence
}
```

### 8.5 CandidateFamily (`candidates.h`)

```
CandidateFamily {
    uid_t    id                // Family ID
    char    *consensus         // Numeric bases (0/1/2/3)
    int      consensus_length
    int      component_id
    int      topology          // TOPO_LINEAR=0, TOPO_COMPLEX=1, TOPO_CYCLIC=2
    freq_t   estimated_copies
    Instance *instances
    int      num_instances
    int      cap_instances
    double   mdl_score         // >0 if accepted
    double   model_cost        // L_int(len) + 2*len
}
```

### 8.6 MDLResult (`mdl.h`)

```
MDLResult {
    double total_model_cost     // Sum of accepted family model costs
    double total_savings        // Sum of accepted family savings
    double dl_library           // L_int(R) + total_model_cost
    double dl_genome_given_lib  // 2*N - total_savings
    double dl_total             // dl_library + dl_genome_given_lib
    double compression_ratio    // dl_total / (2*N)
    int    num_accepted         // Number of accepted families
    int64_t bases_covered       // Genome bases covered by instances
}
```

### 8.7 SeqSegment (`main.c`)

```
SeqSegment {
    int    seq_index      // Index into original genome's sequences
    gpos_t raw_start      // Start in original genome raw coords
    gpos_t raw_end        // End in original genome raw coords
    glen_t seg_length     // = raw_end - raw_start
}
```

---

## 9. Complete Parameter Reference

### 9.1 CLI Parameters

#### Required

| Parameter | Description |
|-----------|-------------|
| `-sequence <file>` | Input FASTA genome |
| `-output <file>` | Output repeat library (FASTA) |

#### Discovery

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-l #` | auto: $\lceil 1 + \log_4 N \rceil$ | L-mer length |
| `-L #` | 10000 | Max extension distance per side (bp) |
| `-minthresh #` | 2 | Min l-mer frequency to seed |
| `-goodlength #` | 30 | Min consensus length pre-filter |
| `-maxgap #` | 5 | Discovery DP band offset (MAXOFFSET) |
| `-match #` | 1 | DP match score |
| `-mismatch #` | -1 | DP mismatch score |
| `-gap #` | -5 | DP gap penalty |
| `-cappenalty #` | -20 | Cap on penalty for exiting alignment |
| `-minimprovement #` | 3 | Min total score improvement per step |
| `-stopafter #` | 100 | Stop after N no-progress positions (WHEN_TO_STOP) |
| `-maxentropy #` | -0.70 | Shannon entropy filter (natural log, negative) |
| `-tandemdist #` | 500 | Min distance between same-strand l-mers |
| `-maxoccurrences #` | 10000 | Max occurrences per seed (MAXN) |
| `-maxrepeats #` | 100000 | Max families to discover (MAXR) |
| `-freq <file>` | — | Pre-computed l-mer frequency table |
| `-freq-output <file>` | — | Write l-mer frequency table |

#### Refinement

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-threads #` | 1 | Thread count |
| `-mdl-mode <mode>` | exact | Position encoding: `none`, `exact`, `upper` |
| `-max-divergence #` | 0.30 | Max substitution rate for instance acceptance |
| `-refine-gap #` | -5 | Refinement gap penalty (-3 recommended for high-indel) |
| `-refine-maxoffset #` | 12 | Refinement DP band half-width (max 32) |
| `-max-dp-cells #` | 10000000 | Max DP cells for merge alignment (~40 MB) |

#### Large Genome

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-chunk-size #` | 200 | Chunk size in Mb (min 10) |
| `-sample-size #` | 1000 | Genome sampling threshold in Mb (min 100) |
| `-window-size #` | 1000 | Sampling tile size in kb (range [100, 10000]) |
| `-seed #` | 42 | Random seed for sampling |
| `-sample-output <file>` | — | Write sampled genome FASTA |

#### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-instances <file>` | — | Instance BED6 output |
| `-stats <file>` | — | Per-family TSV statistics |
| `-v` / `-vv` | off | Verbosity (1 or 2) |

### 9.2 Internal Constants

| Constant | Value | Source | Description |
|----------|-------|--------|-------------|
| PADLENGTH | 11000 | `types.h:22` | Genome padding ($\geq L$) |
| DEFAULT_TANDEMDIST | 500 | `types.h:23` | Tandem distance |
| DEFAULT_MAXN | 10000 | `types.h:24` | Max instances per candidate |
| HASH_SIZE | 16000057 | `discover.c` | Discovery hash table size (prime) |
| LOG2_C0 | 1.5179605508986484 | `mdl.c:16` | Rissanen constant |
| KMER_MAX_POSITIONS | 50000 | `kmer.h:24` | Position array cap |
| NUM_STRIPES | 4096 | `kmer.c:131` | Parallel counting lock stripes |
| POOL_BLOCK_SIZE | 4096 | `kmer.c:10` | KmerEntry pool block |
| REFINE_SCREEN_K | 8 | `refine.h:9` | K-mer screening k |
| PROFILE_BITS | 65536 | `refine.c:18` | 4^8 profile bitset |
| REFINE_MIN_JACCARD | 0.15 | `refine.h:10` | Merge pre-screen |
| REFINE_MATCH | +2 | `refine.h:23` | Merge alignment match |
| REFINE_MISMATCH | -3 | `refine.h:24` | Merge alignment mismatch |
| REFINE_GAP | -2 | `refine.h:25` | Merge alignment gap |
| REFINE_MIN_IDENTITY | 0.80 | `refine.h:13` | 80% identity |
| REFINE_MIN_COVERAGE | 0.80 | `refine.h:14` | 80% coverage |
| REFINE_MIN_ALIGNED | 80 | `refine.h:15` | 80 bp minimum |
| REFINE_RELAXED_IDENTITY | 0.70 | `refine.h:19` | Relaxed merge identity |
| REFINE_RELAXED_COVERAGE | 0.70 | `refine.h:20` | Relaxed merge coverage |
| REFINE_OVERLAP_RELAX | 0.50 | `refine.h:18` | Instance overlap for relaxed merge |
| REFINE_BIMODALITY_THRESH | 0.20 | `refine.h:35` | Split bimodality minimum |
| REFINE_MIN_SPLIT_INSTANCES | 10 | `refine.h:32` | Min instances for split |
| REFINE_MIN_CLUSTER_SIZE | 3 | `refine.h:33` | Min per sub-family |
| REFINE_MIN_DIV_GAP | 0.05 | `refine.h:34` | Min divergence gap |
| REFINE_DIV_BINS | 100 | `refine.h:36` | Otsu histogram bins |
| ALIGN_MATCH | +1 | `align.h:10` | Refinement match |
| ALIGN_MISMATCH | -1 | `align.h:11` | Refinement mismatch |
| ALIGN_CAPPENALTY | -20 | `align.h:12` | Score floor |
| ALIGN_WHEN_TO_STOP | 100 | `align.h:13` | Stall threshold |
| ALIGN_MAX_ITERATIONS | 10 | `align.h:14` | Convergence limit |
| ALIGN_MAXOFFSET_LIMIT | 32 | `align.h:16` | Max `-refine-maxoffset` |
| ALIGN_MAX_EXTENSION | 10000 | `align.h:17` | Max extension bp |
| EXTENSION_SLACK | 15 | `align.h:15` | Min bases from edge |
| MAX_SEED_HITS | 50000 | `align.c:18` | Seed hit cap |
| MAX_CONS_KMERS | 10000 | `align.c:19` | Consensus k-mer cap |
| DISCOVER_SPLIT_THRESHOLD | 200000000 | `main.c:34` | Default chunk size (200 Mb) |
| Fragment assembly D | [500, 5000] | `refine.c:1719-1721` | Proximity distance |
| Fragment co-occurrence | $\geq 3$ | `refine.c:1834` | Min co-occurrence count |
| Fragment nesting guard | 0.50 | `refine.c:1845` | Containment threshold |
| Fragment size ratio | $\geq 0.10$ | `refine.c:1859` | Minimum relative size |
| Fragment orientation | $\geq 0.80$ | `refine.c:1863` | Same-direction fraction |
| Fragment MAD | $\leq 100$ | `refine.c:1996` | Gap consistency |
| Fragment gap range | $[-20, \text{median}+2\text{MAD}]$ | `refine.c:1993-1996` | Allowed gap |
| Fragment chain cap | 10 | `refine.c:1925` | Max union-find chain |
| Prune exclusive filter | 25% | `refine.c:1565` | Min exclusive bases per instance |

---

## 10. Complexity Analysis

### 10.1 Discovery

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| L-mer counting | $O(N)$ | Single genome scan |
| Per-family extension | $O(\text{WHEN\_TO\_STOP} \cdot N_{\text{occ}} \cdot (2\text{MAXOFFSET}+1)^2)$ | Banded DP |
| Masking | $O(N_{\text{occ}} \cdot L_c)$ | 1-vs-1 DP per occurrence |
| Total discovery | $O(R \cdot (N + N_{\text{occ}} \cdot L_c))$ | $R$ families, amortized masking |

### 10.2 K-mer Table

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Counting | $O(N)$ per thread, $O(N/T)$ with $T$ threads | Striped locks |
| Position index | $O(N)$ count + $O(N)$ fill | Two-phase |
| Lookup | $O(1)$ expected | Chaining, low load factor |
| Trimming | $O(T_{\text{size}})$ | One pass over table |

### 10.3 Seeded Alignment

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Seed scan | $O(\sum f_j)$ | Per-family, via position index |
| Clustering | $O(H \log H)$ | Sort $H$ seed hits |
| Banded DP | $O(L_c \cdot W)$ per anchor | $W = 2\text{MAXOFFSET}+1 = 25$ |
| Instances | $O(A \cdot L_c \cdot W)$ per family | $A$ anchors |
| Rebuild | $O(n_{\text{inst}} \cdot L_c)$ | Per iteration |

### 10.4 Refinement

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| K-mer Jaccard | $O(R^2 \cdot P)$ | $P = 1024$ words popcount |
| Merge alignment | $O(R^2 \cdot L_c^2)$ worst case | DP cell limit prevents worst case |
| Otsu split | $O(R \cdot n_{\text{inst}} + R \cdot B)$ | $B = 100$ bins |
| Fragment assembly | $O(I \log I + I \cdot D_{\text{avg}})$ | $I$ = total instances, sweep-line |
| Prune | $O(R \cdot n_{\text{inst}} + N)$ | Coverage array |

### 10.5 Large Genome

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Sampling | $O(n_w)$ | Fisher-Yates partial shuffle |
| Chunk splitting | $O(S)$ | $S$ = number of sequences |
| LPT packing | $O(n_{\text{seg}} \cdot n_{\text{bins}})$ | Simple loop |
| Parallel discovery | $O(\text{largest\_chunk} / T)$ | LPT-balanced |
| Coordinate remapping | $O(I)$ | Per instance |

### 10.6 Profiling Data (human3M, March 2026)

| Function | % Runtime | Role |
|----------|-----------|------|
| `align_collect_instances` | 59.4% | Seeded alignment (main bottleneck) |
| `semiglobal_align` | 10.8% | Merge stage alignment |
| `cmp_seed_hits` | 9.6% | Seed hit sorting |
| `candidates_extract` | 1.6% | Graph component extraction |

---

## Output Formats

### FASTA (`-output`)

```
>R=<id> length=<consensus_length> copies=<num_instances> mdl=<mdl_score>
<consensus in ASCII, 80 chars/line>
```

Only families with `mdl_score > 0` are written. (`output.c:7-35`)

### BED6 (`-instances`)

```
<chr>  <local_start>  <local_end>  R=<id>  <score>  <strand>
```

- Coordinates: chromosome-local (raw_position - chr_offset - PADLENGTH)
- Score: $\lfloor 1000 \times (1 - \text{divergence}) \rfloor$, clamped [0, 1000]
- Only instances of accepted families. (`output.c:37-88`)

### TSV (`-stats`)

```
family_id  consensus_length  num_instances  divergence_mean  mdl_score  model_cost  topology
```

All families written (no filtering). Topology: "linear", "cyclic", or "complex". (`output.c:90-123`)

---

## References

1. Rissanen, J. (1978). Modeling by shortest data description. *Automatica*, 14(5), 465-471.
2. Price, A.L., Jones, N.C., & Pevzner, P.A. (2005). De novo identification of repeat families in large genomes. *Bioinformatics*, 21(suppl_1), i351-i358.
3. Wicker, T. et al. (2007). A unified classification system for eukaryotic transposable elements. *Nature Reviews Genetics*, 8(12), 973-982.
4. Otsu, N. (1979). A threshold selection method from gray-level histograms. *IEEE Trans. SMC*, 9(1), 62-66.
5. Grumbach, S., & Tahi, F. (1994). A new challenge for compression algorithms: genetic sequences. *Information Processing & Management*, 30(6), 875-886.
