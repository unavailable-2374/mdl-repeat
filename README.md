# mdl-repeat

**MDL-guided de novo repeat element discovery.**

mdl-repeat identifies repetitive element families in genomic DNA sequences without prior knowledge of a repeat library. It uses the **Minimum Description Length (MDL) principle** — an information-theoretic criterion that automatically determines how many families to report, the minimum consensus length, and the minimum copy number that justifies including a family in the library. No manual threshold tuning is required.

## Features

- **MDL-based model selection** — a family is accepted if and only if encoding its instances as references to a consensus saves more bits than encoding the consensus itself
- **Seed-and-extend discovery** — high-frequency l-mer seeds are extended into consensus sequences via N-sequence simultaneous banded DP alignment
- **Refinement pipeline** — merge redundant families (80-80-80 rule), split bimodal divergence (Otsu's method), assemble TE fragments (spatial co-occurrence), prune marginal families (exclusive coverage test)
- **Three MDL encoding modes** — `exact` (lgamma binomial, default), `upper` (conservative bound), `none`
- **Multi-threading** — parallel k-mer position index construction, alignment refinement, and chunked discovery
- **Large genome support** — genome sampling (tile-based Fisher-Yates) and chunked parallel discovery with LPT bin packing for multi-gigabase genomes
- **Frequency table reuse** — save/load l-mer frequency tables to skip recomputation across runs

## Requirements

- GCC (C11 standard)
- GNU Make
- Linux/POSIX
- No external libraries at build time (only `libm` and pthreads)

**Optional runtime dependency:**
- `rmblastn` (RMBlast) — for BLAST-based instance recruitment of short
  families (consensus length < 500 bp). When present in PATH or at
  `/home/shuoc/tool/miniconda3/envs/PGTA/bin/rmblastn`, the repeat-aware
  `rmblastn` task is used with a small word seed (`-word_size 7`),
  RepeatMasker-style default gap costs (`-gapopen 12 -gapextend 2`),
  `-perc_identity 70`, and `-qcov_hsp_perc 60` for bounded short-element
  recruitment. Falls back to the
  built-in k-mer + banded DP path if rmblastn is not found. Install via:
  `conda install -c bioconda rmblast`.

## Build

```bash
make                # Build bin/mdl-repeat (with -march=native -flto)
make PORTABLE=1     # Build without -march=native (for portable binaries)
make clean          # Remove build artifacts
```

## Quick Start

```bash
# Discover repeat families from a genome
./bin/mdl-repeat -sequence genome.fa -output families.fa -v

# Save l-mer frequency table for reuse across runs
./bin/mdl-repeat -sequence genome.fa -output families.fa -freq-output genome.freq
./bin/mdl-repeat -sequence genome.fa -output families.fa -freq genome.freq

# Full output with instance locations and statistics
./bin/mdl-repeat -sequence genome.fa \
    -output families.fa \
    -instances instances.bed \
    -stats stats.tsv \
    -threads 4 -vv

# Large genome (auto-samples if > 1 Gb, auto-chunks if > 200 Mb)
./bin/mdl-repeat -sequence human.fa -output families.fa \
    -threads 8 -chunk-size 200 -sample-size 1000 -v

# Very large genome with custom sampling
./bin/mdl-repeat -sequence axolotl.fa -output families.fa \
    -threads 16 -sample-size 2000 -window-size 500 -seed 123 \
    -sample-output sampled.fa -v
```

## Usage

```
mdl-repeat -sequence <file> -output <file> [options]

Required:
  -sequence <file>   Input FASTA (single or multi-sequence)
  -output <file>     Output repeat library (FASTA)

Discovery:
  -freq <file>       Pre-computed l-mer frequency table
  -freq-output <file> Write l-mer frequency table for reuse
  -l #               L-mer length (default: auto = ceil(1+log4(N)))
  -L #               Max extension distance per side (default: 10000)
  -minthresh #       Min l-mer frequency to seed (default: 2)
  -goodlength #      Min consensus length pre-filter (default: 30)
  -maxgap #          Max DP band offset (default: 5)
  -match #           Match score (default: 1)
  -mismatch #        Mismatch score (default: -1)
  -gap #             Gap penalty (default: -5)
  -cappenalty #      Cap penalty (default: -20)
  -minimprovement #  Min improvement per step (default: 3)
  -stopafter #       Stop after N no-progress positions (default: 100)
  -maxentropy #      Entropy filter threshold (default: -0.70)
  -tandemdist #      Min distance between same-strand l-mers (default: 500)
  -maxoccurrences #  Max occurrences per seed (default: 10000)
  -maxrepeats #      Max families to discover (default: 100000)

Refinement:
  -threads #              Number of threads (default: 1)
  -mdl-mode <mode>        MDL position encoding: none|exact|upper (default: exact)
  -max-divergence #       Max substitution rate (default: 0.30)
  -refine-gap #           Refinement gap penalty (default: -5, try -3 for high-indel)
  -refine-maxoffset #     Refinement DP band half-width (default: 12, max: 32)
  -max-dp-cells #         Max DP cells for merge alignment (default: 10000000)
  -coalesce-factor #      Gap tolerance for tandem-instance coalescing in bases
                          (default: 20.0; 0 = disabled)

Large genome:
  -chunk-size #           Chunk size in Mb for parallel discovery (default: 200)
  -sample-size #          Sampling threshold in Mb (default: 1000)
  -window-size #          Sampling tile size in kb (default: 1000)
  -seed #                 Random seed for sampling (default: 42)
  -sample-output <file>   Write sampled genome FASTA

Output:
  -instances <file>  Instance locations (BED6)
  -stats <file>      Per-family statistics (TSV)
  -v / -vv           Verbosity level
```

## Pipeline

```
Input FASTA
    |
    v
1. Load genome
    |
    v
1b. Sample genome (if > sample_size, default 1 Gb)
    - Tile-based Fisher-Yates random sampling
    |
    v
2. Discover consensus families
   - If > chunk_size (default 200 Mb): chunked parallel discovery
     with LPT bin packing and overlapping segments
   - Otherwise: single-pass seed-and-extend discovery
    |
    v
3. Build k-mer table and position index
    |
    v
4. Compact (remove dead/short families)
    |
    v
5. Merge redundant families (80-80-80 rule + union-find)
    |
    v
6. Split bimodal families (Otsu's method, MDL-validated)
    |
    v
7. Assemble TE fragments (spatial co-occurrence, nesting guard)
    |
    v
8. MDL scoring and library selection (iterative R convergence)
    |
    v
9. Prune marginal families (exclusive coverage vs. model cost)
    |
    v
10. Recovery pass (re-score rejected families with reduced R)
    |
    v
11. Output (remap coordinates if sampled)
    → FASTA library + BED instances + TSV statistics
```

## How It Works

### MDL Model Selection

The Minimum Description Length principle minimizes the total cost of describing the genome:

```
DL(Genome) = DL(Library) + DL(Genome | Library)
```

- **Library cost**: each consensus sequence is encoded at 2 bits/base plus a length header using Rissanen's universal integer code
- **Instance encoding**: each occurrence is encoded as edits from its consensus — cheaper than literal 2 bits/base when divergence is low
- **Selection criterion**: a family is accepted if the total savings across all its instances exceeds the cost of storing its consensus

The number of accepted families R and per-instance costs are circularly dependent (costs depend on R, but R depends on costs). This is resolved by iterating the scoring up to 3 times until R converges.

### Discovery Engine

The discovery engine identifies repeat families through seed-and-extend:

1. **L-mer frequency counting** — count all l-mers in the genome using symmetric hashing (forward and reverse complement share the same bucket), with tandem distance filtering to avoid inflated counts from tandem arrays
2. **Seed selection** — greedily select the highest-frequency unmasked l-mer as the next seed
3. **N-sequence simultaneous extension** — align all N occurrences of the seed simultaneously using banded dynamic programming, extending the consensus one position at a time in both directions; at each position, the consensus base is chosen by score-weighted majority vote across all occurrences
4. **Masking** — mask the discovered family's occurrences in the genome to prevent rediscovery, then repeat from step 2

The default l-mer length is `ceil(1 + log4(N))` where N is the genome length.

### Refinement Pipeline

- **Merge**: Detects redundant families using k-mer Jaccard pre-screening followed by semi-global alignment. Families sharing >= 80% identity over >= 80% coverage with >= 80 bp aligned are merged. Transitive merges are resolved via union-find.
- **Split**: Detects families containing mixed subfamilies by analyzing divergence distributions with Otsu's thresholding method. Splits are accepted only when the combined MDL score of the two sub-families exceeds the original.
- **Fragment assembly**: Detects adjacent TE fragments that were discovered independently. Uses a spatial co-occurrence sweep-line to find families whose instances consistently appear near each other in the genome. A nesting guard prevents merging nested elements. Assembly is accepted only if the assembled MDL score exceeds the sum of the parts.
- **Prune**: Removes accepted families whose exclusive genome coverage (positions not covered by any other family) does not justify their model cost.
- **Recovery**: After pruning reduces R, the per-instance overhead decreases, which may cause previously rejected families to become viable. A recovery pass re-scores all rejected families with the final R.

### Large Genome Support

For multi-gigabase genomes, two scaling mechanisms activate automatically:

- **Genome sampling** (genome > 1 Gb): The genome is reduced to a representative sample via tile-based random selection (1 Mb non-overlapping tiles, Fisher-Yates shuffle). Instance coordinates are remapped to the original genome after discovery. Controlled via `-sample-size`, `-window-size`, `-seed`.
- **Chunked discovery** (genome > 200 Mb): Long sequences are segmented with overlap at boundaries, packed into balanced bins via LPT scheduling, and processed in parallel. Each chunk computes its own l-mer length for increased sensitivity. Controlled via `-chunk-size`.

All position and length variables use 64-bit integers (`gpos_t`/`glen_t`) to avoid overflow on large genomes.

## Output Formats

**FASTA** (`-output`): Consensus sequences for accepted repeat families.

```
>R=0 length=312 copies=15 mdl=1204.3
ACGTACGT...
```

**BED6** (`-instances`): Genomic coordinates of every instance.

```
chr1    10432    10744    R=0    850    +
```

Column 5 is `1000 * (1 - divergence)` on a 0-1000 scale.

**TSV** (`-stats`): Per-family summary (consensus length, copy count, mean divergence, MDL score, model cost, topology).

## Architecture

```
src/
├── main.c              Pipeline driver, CLI parsing, chunked discovery, genome sampling
├── types.h             Core types (gpos_t/glen_t as int64_t, DNA encoding, constants)
├── genome.c/.h         FASTA loading, padding, boundary tracking
├── kmer.c/.h           Canonical k-mer counting, striped-lock parallel hash table, position index
├── discover.c/.h       Seed-and-extend discovery engine (thread-safe DiscoverContext)
├── align.c/.h          Multi-k-mer seeding, banded DP alignment, consensus rebuild/extend
├── candidates.c/.h     CandidateFamily/Instance structures, memory management
├── refine.c/.h         Merge, split, fragment assembly, prune
├── mdl.c/.h            Rissanen's universal integer code, MDL scoring, library selection
├── output.c/.h         FASTA, BED6, TSV writers
└── cmd_line_opts.c/.h  Argument parsing
```

~8,300 lines of C. No external dependencies.

## Testing

```bash
# Generate synthetic test genomes
python3 tests/gen_synthetic.py -o tests/data

# Run full test suite
bash tests/run_tests.sh
```

Test suite includes:
- Multi-family detection (testA)
- Tandem repeat handling (testB)
- Detection limit / sensitivity (testC)
- Nested TE detection (testD)
- L_int unit tests

## References

- Rissanen, J. (1978). Modeling by shortest data description. *Automatica*, 14(5), 465-471.
- Price, A.L., Jones, N.C., & Pevzner, P.A. (2005). De novo identification of repeat families in large genomes. *Bioinformatics*, 21(suppl_1), i351-i358.
- Wicker, T. et al. (2007). A unified classification system for eukaryotic transposable elements. *Nature Reviews Genetics*, 8(12), 973-982.
- Otsu, N. (1979). A threshold selection method from gray-level histograms. *IEEE Trans. SMC*, 9(1), 62-66.

## License

GPL-2.0
