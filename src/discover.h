#ifndef MDL_DISCOVER_H
#define MDL_DISCOVER_H

#include "types.h"
#include "genome.h"
#include "candidates.h"

/*
 * RepeatScout-style seed-and-extend discovery engine.
 *
 * Faithfully replicates RepeatScout v1.0.7 N-sequence simultaneous
 * banded DP extension + mask-and-extend, with bug fixes:
 *   - Sequence boundary check (removes hardcoded bStart=0; bEnd=50000000)
 *   - Dynamic value[] allocation (fixes overflow when l > 19)
 *
 * MDL replaces GOODLENGTH/MINTHRESH as final arbiter.
 */

/* Discovery parameters */
typedef struct {
    int   l;              /* l-mer length (0 = auto: ceil(1+log4(N))) */
    int   L;              /* max extension distance per side */
    int   MAXOFFSET;      /* max DP band offset */
    int   MAXN;           /* max occurrences per l-mer seed */
    int   MAXR;           /* max families to report */
    int   MATCH;          /* alignment match score */
    int   MISMATCH;       /* alignment mismatch score */
    int   GAP;            /* alignment gap score */
    int   CAPPENALTY;     /* cap on penalty for exiting alignment */
    int   MINIMPROVEMENT; /* min totalbestscore improvement per step */
    int   WHEN_TO_STOP;   /* stop after N no-progress positions */
    float MAXENTROPY;     /* Shannon entropy filter (natural log, negative) */
    int   GOODLENGTH;     /* min consensus length pre-filter */
    int   MINTHRESH;      /* min l-mer frequency to seed */
    int   TANDEMDIST;     /* min distance between counted same-strand l-mers */
    int   VERBOSE;        /* verbosity level */
    char *freq_file;      /* .freq file path (NULL = count internally) */
    char *freq_output;    /* output .freq file path (NULL = don't write) */
} DiscoverParams;

/* Initialize parameters with defaults matching plan */
void discover_default_params(DiscoverParams *params);

/*
 * Main discovery function.
 *
 * Runs the seed-and-extend loop:
 *   1. Build/read l-mer frequency table
 *   2. For each high-frequency seed: extend right, extend left
 *   3. Mask found families to prevent rediscovery
 *   4. Collect instances from extension occurrences
 *
 * Returns a CandidateList with consensus sequences and instances.
 * Caller must free with candidates_free().
 * Returns NULL on failure.
 */
CandidateList *discover_families(const Genome *genome,
                                 const DiscoverParams *params);

#endif /* MDL_DISCOVER_H */
