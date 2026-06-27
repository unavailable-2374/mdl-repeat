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

/*
 * Consensus-vs-consensus alignment scoring (used ONLY by the merge stage
 * to decide whether two candidate families represent the same biological
 * repeat — this is a binary "same family?" decision).
 *
 * SCOPE BOUNDARY (M2#7): These constants do NOT feed into MDL.  The
 * per-instance edit count m_i that mdl.c consumes comes from
 * align_collect_instances() in align.c, which uses ALIGN_MATCH /
 * ALIGN_MISMATCH / g_align_gap (1 / -1 / -5).  The two scoring systems
 * are deliberately decoupled:
 *
 *   - align.c  (1 / -1 / -5): instance-vs-consensus banded DP, produces
 *                              num_edits — the substrate for MDL scoring.
 *   - refine.c (2 / -3 / -2): consensus-vs-consensus semi-global DP,
 *                              produces only an identity/coverage % for
 *                              the 80-80-80 same-family threshold.
 *                              More aggressive match preference because
 *                              the decision is binary and the alignment
 *                              is short relative to consensi.
 */
#define REFINE_MATCH    2
#define REFINE_MISMATCH -3
#define REFINE_GAP      -2

/* Runtime-configurable: max DP matrix cells (default 10M ~40MB).
 * Increase to 50M for long TE families (ERV, Helitron). */
extern int64_t g_refine_max_dp_cells;

/* Optional split-stage audit TSV path.  NULL disables audit output. */
extern const char *g_refine_split_audit_path;

/* Splitting parameters */
#define REFINE_MIN_SPLIT_INSTANCES  3    /* minimum instances (chr4 push experiment 2026-05-02; was 5, was 10) */
#define REFINE_MIN_CLUSTER_SIZE      3   /* minimum instances per sub-family */
#define REFINE_MIN_DIV_GAP        0.03f  /* min mean divergence difference for split (lowered from 0.05 — chr4 90×80 experiment 2026-05-02) */
#define REFINE_BIMODALITY_THRESH  0.20f  /* min inter-class/total variance ratio */
#define REFINE_DIV_BINS            100   /* histogram bins for Otsu's method */

/* Fragment-assembly co-occurrence FRACTION guard (#3 chimera defense): on top
 * of the absolute >=3 co-occurrence floor, require that this fraction of the
 * smaller family's instances co-occur with the other family.  Blocks fusing
 * unrelated-but-adjacent high-copy TEs (3-of-thousands is incidental); real
 * fragments of one element co-occur in most copies. */
#define REFINE_ASSEMBLE_MIN_COOCC_FRAC 0.20f

/* Fragment-assembly divergence-quality guard (#4 chimera defense): the MDL
 * "joined > sum-of-parts" gate is volume-driven — savings scale with copy
 * count, so a chimeric join of two unrelated high-copy families passes with a
 * large positive score while its "copies" fit the concatenated consensus
 * poorly. A genuine assembly of two pieces of ONE element keeps per-copy
 * divergence in the real-TE band (<= max-divergence 0.30); a chimera spikes to
 * 0.4-0.8. Reject the join when the mean divergence of the merged instances
 * exceeds this cap. Set just above the 0.30 refine max-divergence band. */
#define REFINE_ASSEMBLE_MAX_DIV 0.35f

/* Chimeric/over-extended long-family filter (post-selection cleanup).
 * Long consensi (assembly chains, discovery over-extensions to the L-cap,
 * coalesced tandems) whose copies fit the consensus poorly are chimeras: a
 * genuine long element is conserved across copies (low mean divergence), while
 * a chimera's "copies" only match a fragment, scoring high mean divergence (the
 * unmatched bulk counts as edits — divergence here is the per-copy coverage
 * signal in another guise). Length is only the SCOPE (high divergence is normal
 * and harmless in short divergent families); the DECISION is divergence, not
 * length. Empirically (TAIR10) this drops exactly the 40 kb+ chimeras and
 * 20 kb-cap over-extensions with zero family-level recall change. */
#define REFINE_CHIMERA_MIN_LEN 8000
#define REFINE_CHIMERA_MAX_DIV 0.40f

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
 * Uses Otsu's method to detect bimodality.  The current relaxed gate accepts
 * non-negative combined split MDL for positive-score originals; quality
 * classification is handled downstream.
 * Called before MDL selection (operates on all families).
 * num_families is the current R estimate for MDL scoring.
 * Returns number of splits performed.
 */
int refine_split_families(CandidateList *cl, const Genome *genome,
                          const KmerTable *kt, int k,
                          glen_t genome_len, int verbose,
                          int num_families, int num_threads);

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

/*
 * Coalesce same-family tandem-array instances into single longer instances.
 *
 * Background: RM-style ground-truth annotations merge directly-adjacent
 * repeat copies into one long interval (e.g. a 6×1.6 kb tandem array is
 * reported as one 9.6 kb truth interval).  mdl-repeat naturally outputs
 * each copy as a separate instance, which produces a fragment-versus-truth
 * mismatch that drives apparent recall down.
 *
 * This pass walks each accepted family's instance list and merges
 * consecutive instances on the same strand whose gap is less than
 * coalesce_factor × consensus_length (default 1.5).  num_edits and
 * divergence are summed/averaged proportionally; aligned_length grows.
 *
 * Net effect on reporting: long tandem arrays become single long
 * instances that match the truth annotation convention.  No effect on
 * MDL scoring (already done before this pass).
 *
 * Returns: number of instance pairs coalesced.
 */
int refine_coalesce_tandem_instances(CandidateList *cl,
                                     float coalesce_factor,
                                     int verbose);

/*
 * Drop accepted families that are long AND high-divergence — chimeric joins,
 * discovery over-extensions to the L-cap, and coalesced tandems whose copies
 * do not actually match the long consensus. Runs after MDL selection / prune
 * (operates on accepted families; sets accept_state = REJECTED so output skips
 * them). See REFINE_CHIMERA_MIN_LEN / REFINE_CHIMERA_MAX_DIV. Recall-neutral.
 * Returns the number of families dropped.
 */
int refine_drop_chimeric_long(CandidateList *cl, int verbose);

#endif /* MDL_REFINE_H */
