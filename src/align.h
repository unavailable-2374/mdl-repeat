#ifndef MDL_ALIGN_H
#define MDL_ALIGN_H

#include "types.h"
#include "candidates.h"
#include "genome.h"
#include "kmer.h"

/* Alignment parameters (same as RepeatScout defaults) */
#define ALIGN_MATCH           1
#define ALIGN_MISMATCH       -1
#define ALIGN_GAP            -5
#define ALIGN_MAXOFFSET      12
#define ALIGN_CAPPENALTY    -20
#define ALIGN_WHEN_TO_STOP  100
#define ALIGN_MAX_DIVERGENCE  0.30f
#define ALIGN_MAX_ITERATIONS 10
#define EXTENSION_SLACK      15  /* bases from edge to still contribute to extension */

/* Seed hit from multi-k-mer scanning */
typedef struct {
    gpos_t genome_pos;       /* position in genome (padded coords) */
    int    cons_pos;         /* position in consensus */
    int8_t strand;           /* +1 forward, -1 reverse */
} SeedHit;

/* Banded alignment result for one instance */
typedef struct {
    gpos_t genome_start;     /* start in genome (padded coords) */
    gpos_t genome_end;       /* end in genome */
    int    cons_start;       /* start in consensus (for truncation) */
    int    cons_end;         /* end in consensus (exclusive) */
    int    num_edits;        /* actual edit count */
    float  divergence;       /* num_edits / aligned_length */
    int    score;            /* alignment score */
    int8_t strand;
    int    seq_index;
} AlignedInstance;

/* --- API --- */

/*
 * Multi-k-mer seeded instance collection + banded alignment.
 * Replaces the single-anchor collect_instances() from MVP.
 * Modifies fam->instances in place.
 */
int align_collect_instances(CandidateFamily *fam, const Genome *genome,
                            const KmerTable *kt, int k);

/*
 * Rebuild consensus from aligned instances via majority voting.
 * Returns number of positions changed. Updates fam->consensus.
 */
int align_rebuild_consensus(CandidateFamily *fam, const Genome *genome);

/*
 * Extend consensus left/right using flanking genome context from instances.
 * Returns total bases extended (left + right).
 */
int align_extend_consensus(CandidateFamily *fam, const Genome *genome);

/*
 * Full refinement loop: seed -> align -> rebuild -> iterate.
 * Calls align_collect_instances + align_rebuild_consensus in a loop.
 * max_iter iterations (default ALIGN_MAX_ITERATIONS).
 * Returns final number of instances.
 */
int align_refine_family(CandidateFamily *fam, const Genome *genome,
                        const KmerTable *kt, int k, int max_iter);

/*
 * Refine all families in a candidate list.
 * Removes families that lose all instances after refinement.
 * Uses num_threads parallel workers (1 = sequential).
 */
void align_refine_all(CandidateList *cl, const Genome *genome,
                      const KmerTable *kt, int k, int num_threads,
                      int verbose);

#endif /* MDL_ALIGN_H */
