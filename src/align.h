#ifndef MDL_ALIGN_H
#define MDL_ALIGN_H

#include "types.h"
#include "candidates.h"
#include "genome.h"
#include "kmer.h"

/* Alignment parameters (same as RepeatScout defaults) */
#define ALIGN_MATCH           1
#define ALIGN_MISMATCH       -1
#define ALIGN_CAPPENALTY    -20
#define ALIGN_WHEN_TO_STOP  100
#define ALIGN_MAX_ITERATIONS 10
#define EXTENSION_SLACK      15  /* bases from edge to still contribute to extension */

/* Maximum allowed MAXOFFSET for stack array sizing (dp[2][2*MAX+1]) */
#define ALIGN_MAXOFFSET_LIMIT  32

/* Runtime-configurable refinement parameters (set via CLI in main.c) */
extern int   g_align_gap;            /* default: -5 */
extern int   g_align_maxoffset;      /* default: 12, max: ALIGN_MAXOFFSET_LIMIT */
extern float g_align_max_divergence; /* default: 0.30 */
extern int   g_align_max_instances;  /* default: 10000 (per-family recruit cap) */
extern int   g_align_max_seed_hits;  /* default: 50000 (raw seed-hit cap)       */

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
 * For families with consensus_length < 500 bp, tries RMBlast first
 * (more sensitive for short divergent elements); falls back to banded DP
 * if rmblastn is not in PATH.
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

/*
 * RMBlast-based batch instance recruitment for short families
 * (consensus_length < RMBLAST_SHORT_THRESHOLD = 500 bp).
 *
 * Runs ONE rmblastn job with all short-family consensuses as a single
 * multi-FASTA query file (query ID = family index).  Parses BLAST output
 * and distributes hits to the appropriate CandidateFamily.  Much faster
 * than per-family BLAST calls.
 *
 * Called once from main.c after align_refine_all, as a supplementary
 * recruitment step.  Families that already have enough instances are
 * skipped (their instances are merged with BLAST hits).
 *
 * Returns total instances added across all families, or -1 if rmblastn
 * is unavailable (caller should log and skip).
 */
int align_blast_recruit_short_families(CandidateList *cl,
                                        const Genome *genome, int k,
                                        int num_threads, int verbose);

#endif /* MDL_ALIGN_H */
