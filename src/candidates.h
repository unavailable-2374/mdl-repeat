#ifndef MDL_CANDIDATES_H
#define MDL_CANDIDATES_H

#include "types.h"
#include "genome.h"

/* Topology classification */
#define TOPO_LINEAR  0
#define TOPO_COMPLEX 1
#define TOPO_CYCLIC  2

/* A single instance (occurrence) of a candidate in the genome */
typedef struct {
    gpos_t position;          /* start position in padded genome coordinates */
    glen_t aligned_length;    /* approximate alignment length */
    int    cons_start;        /* start in consensus (0-based) */
    int    cons_end;          /* end in consensus (exclusive) */
    int    num_edits;         /* actual edit count from alignment */
    float  divergence;        /* estimated divergence (0.0 - 1.0) */
    int    score;             /* alignment score */
    int8_t strand;            /* +1 forward, -1 reverse */
    int    seq_index;         /* which FASTA sequence it's in */
} Instance;

/* A candidate repeat family */
typedef struct {
    uid_t    id;
    char    *consensus;       /* numeric bases (0/1/2/3) */
    int      consensus_length;
    int      component_id;
    int      topology;        /* TOPO_LINEAR, TOPO_COMPLEX, TOPO_CYCLIC */
    freq_t   estimated_copies; /* estimated copy number */

    Instance *instances;
    int       num_instances;
    int       cap_instances;

    /* MDL fields (filled by mdl_score_family) */
    double   mdl_score;
    double   model_cost;
} CandidateFamily;

typedef struct {
    CandidateFamily *families;
    int              num_families;
    int              cap_families;
} CandidateList;

/*
 * Free candidate list and all associated memory.
 */
void candidates_free(CandidateList *cl);

/*
 * Print candidate statistics.
 */
void candidates_print_stats(const CandidateList *cl);

#endif /* MDL_CANDIDATES_H */
