#ifndef MDL_DISCOVER_INTERNAL_H
#define MDL_DISCOVER_INTERNAL_H

/*
 * Private header — types and helpers shared between the discover.c
 * core and split-out modules (currently: discover_mask.c).
 *
 * This header MUST NOT be included outside the discover_*.c family.
 * It is part of the M3#8 split that pulled mask_headptr() and its DP
 * helpers into discover_mask.c.  Future cuts (extension subsystem,
 * seed-table subsystem) can keep adding declarations here.
 */

#include "types.h"
#include "discover.h"
#include "candidates.h"

/* ================================================================
 * Constants
 * ================================================================ */

#define HASH_SIZE_MIN   16000057   /* prime, same as RepeatScout — minimum */
#define SMALLHASH_SIZE  5003       /* prime, for mask consensus lookup */
#define SMALLL          6          /* length of small hash l-mers */
#define EXTRAMASK       1          /* extend mask region by l-1 on each side */

/*
 * Mask claim accounting.  Currently uses binary semantics
 * (MAX_FAMILY_CLAIMS = 1) — this matches the original RepeatScout
 * design and was empirically the best on TAIR10 chr4 (a value of 2
 * caused over-fragmentation of long families: bp F1 0.632 → 0.609).
 * The constant is retained so that future experiments with adjusted
 * mask logic can flip it without repeatedly editing every check site.
 *
 * CLAIM_PERMANENT is a sentinel for build_all_pos's TANDEMDIST-pruning,
 * which marks positions as never-available regardless of family claims.
 */
#define MAX_FAMILY_CLAIMS  1
#define CLAIM_PERMANENT    255

/* ================================================================
 * Internal hash table types (RepeatScout convention)
 * ================================================================ */

struct posllist {
    gpos_t this;
    struct posllist *next;
};

struct llist {
    int freq;
    gpos_t lastocc;
    gpos_t lastplusocc;       /* last forward-strand occurrence position */
    gpos_t lastminusocc;      /* last reverse-complement occurrence position */
    struct llist *next;
    struct posllist *pos;
};

struct repeatllist {
    char *value;              /* dynamically allocated, length = l */
    int   repeatpos;
    struct repeatllist *next;
};

/* ================================================================
 * DiscoverContext — all mutable state for one discover_families() call
 * ================================================================ */

typedef struct DiscoverContext {
    /* Parameters (from DiscoverParams) */
    gpos_t hash_size;              /* dynamic: max(HASH_SIZE_MIN, 4*genome_len/l) */
    int l, L, MAXOFFSET, MAXN, MAXR;
    int MATCH, MISMATCH, GAP, CAPPENALTY;
    int MINIMPROVEMENT, WHEN_TO_STOP;
    float MAXENTROPY;
    int GOODLENGTH, MINTHRESH, TANDEMDIST, VERBOSE;

    /* Genome data */
    char    *sequence;
    char    *sequence_owned;
    char    *removed;
    glen_t   length;
    glen_t   orig_length;
    gpos_t  *disc_boundaries;
    int      disc_num_sequences;

    /* Extension workspace */
    char    *master;
    char   **masters;
    int     *masters_allocated;
    int     *masterstart;
    int     *masterend;

    gpos_t  *pos;
    char    *rev;
    int     *upperBoundI;
    int      N;

    int   ***score;
    int    **score_of_besty;
    int    **maskscore;

    int      totalbestscore;
    int      besttotalbestscore;
    int     *bestbestscore;
    int     *savebestscore;

    int     *best_left_offset;
    int     *save_left_offset;
    int     *best_right_offset;
    int     *save_right_offset;

    int      besty, bestw;
    int      nrepeatocc, nactiverepeatocc;
    int      bestnrepeatocc, bestnactiverepeatocc;
    int      R;
    int      prevbestfreq, prevbesthash;
} DiscoverContext;

/* ================================================================
 * Cross-module helpers — defined in discover.c, used elsewhere
 * ================================================================ */

char compl_base(char c);
int  hash_function(DiscoverContext *C, const char *lmer);
int  smallhash_function(const char *lmer);
int  lmermatch(DiscoverContext *C, const char *l1, const char *l2);
int  lmermatchrc(DiscoverContext *C, const char *l1, const char *l2);
int  lmermatcheither(DiscoverContext *C, const char *l1, const char *l2);

/* ================================================================
 * Mask subsystem (defined in discover_mask.c)
 * ================================================================ */

void mask_headptr(DiscoverContext *C, struct llist **headptr);

#endif /* MDL_DISCOVER_INTERNAL_H */
