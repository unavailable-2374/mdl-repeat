#ifndef MDL_REFINE_H
#define MDL_REFINE_H

#include "candidates.h"
#include "genome.h"
#include "kmer.h"

/* K-mer screening parameters */
#define REFINE_SCREEN_K        8
#define REFINE_MIN_JACCARD     0.15f

/* 80-80-80 rule (Wicker et al. 2007) */
#define REFINE_MIN_IDENTITY    0.80f
#define REFINE_MIN_COVERAGE    0.80f
#define REFINE_MIN_ALIGNED     80

/* Instance overlap relaxation */
#define REFINE_OVERLAP_RELAX   0.50f
#define REFINE_RELAXED_IDENTITY 0.70f
#define REFINE_RELAXED_COVERAGE 0.70f

/* Consensus-vs-consensus alignment scoring */
#define REFINE_MATCH    2
#define REFINE_MISMATCH -3
#define REFINE_GAP      -2

/* Max DP matrix cells to prevent memory/time blowup on long consensus pairs */
#define REFINE_MAX_DP_CELLS  (10 * 1000 * 1000)  /* 10M cells ~40MB */

/* Splitting parameters */
#define REFINE_MIN_SPLIT_INSTANCES  10   /* minimum instances to attempt split */
#define REFINE_MIN_CLUSTER_SIZE      3   /* minimum instances per sub-family */
#define REFINE_MIN_DIV_GAP        0.05f  /* min mean divergence difference for split */
#define REFINE_BIMODALITY_THRESH  0.20f  /* min inter-class/total variance ratio */
#define REFINE_DIV_BINS            100   /* histogram bins for Otsu's method */

/*
 * Merge redundant families based on consensus similarity (80-80-80 rule).
 * Handles 5' truncation via semi-global alignment and checks both strands.
 * Removes absorbed families and compacts the candidate list.
 * Returns number of merges performed.
 */
int refine_merge_families(CandidateList *cl, const Genome *genome,
                          const KmerTable *kt, int k, int verbose,
                          int num_threads);

/*
 * Split families with bimodal divergence distribution into subfamilies.
 * Uses Otsu's method to detect bimodality, accepts split only if MDL improves.
 * Called before MDL selection (operates on all families).
 * num_families is the current R estimate for MDL scoring.
 * Returns number of splits performed.
 */
int refine_split_families(CandidateList *cl, const Genome *genome,
                          const KmerTable *kt, int k,
                          glen_t genome_len, int verbose,
                          int num_families);

/*
 * Prune marginal families after MDL selection.
 * For each accepted family (weakest first), checks if its exclusive genome
 * coverage justifies its model cost. Removes families where it does not.
 * num_families is the current R for MDL per-instance overhead.
 * Returns number of families pruned.
 */
int refine_prune_families(CandidateList *cl, glen_t genome_len, int verbose,
                          int num_families);

/*
 * Assemble fragment families that are adjacent parts of the same TE.
 * Uses spatial co-occurrence (NOT k-mer Jaccard) to find candidates.
 * Guards against merging nested elements.
 * Validates assemblies via MDL: assembled score must exceed sum of parts.
 * Returns number of assemblies performed.
 */
int refine_assemble_fragments(CandidateList *cl, const Genome *genome,
                              const KmerTable *kt, int k,
                              glen_t genome_len, int verbose,
                              int num_threads);

#endif /* MDL_REFINE_H */
