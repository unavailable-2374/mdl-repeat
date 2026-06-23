#ifndef MDL_CANDIDATES_H
#define MDL_CANDIDATES_H

#include <stdio.h>
#include <string.h>

#include "types.h"
#include "genome.h"

/* Topology classification */
#define TOPO_LINEAR  0
#define TOPO_COMPLEX 1
#define TOPO_CYCLIC  2

/* Final library-selection state.  Keep zero as "not yet selected" so
 * memset-initialized families from discovery/refine remain valid. */
typedef enum {
    CAND_ACCEPT_UNSCORED = 0,
    CAND_ACCEPT_REJECTED,
    CAND_ACCEPT_EXCLUSIVE,
    CAND_ACCEPT_STANDALONE,
    CAND_ACCEPT_PRUNED
} CandidateAcceptState;

/* Non-destructive quality tier.  Early phases only annotate; filtering by
 * tier must remain opt-in until benchmarked. */
typedef enum {
    CAND_TIER_UNSET = 0,
    CAND_TIER_CORE,
    CAND_TIER_SUPPORTED,
    CAND_TIER_RESCUE,
    CAND_TIER_WARN,
    CAND_TIER_REJECT
} CandidateQualityTier;

enum {
    CAND_QF_NONE                 = 0,
    CAND_QF_SHORT_CONSENSUS      = 1u << 0,
    CAND_QF_LOW_COPY             = 1u << 1,
    CAND_QF_HIGH_DIVERGENCE      = 1u << 2,
    CAND_QF_STANDALONE_FALLBACK  = 1u << 3,
    CAND_QF_NO_EXCLUSIVE_BASES   = 1u << 4,
    CAND_QF_LOW_EXCLUSIVE_FRAC   = 1u << 5,
    CAND_QF_NONPOSITIVE_EXCL_MDL = 1u << 6,
    CAND_QF_PRUNED_REDUNDANT     = 1u << 7,
    CAND_QF_RESCUE_DISCOVERY     = 1u << 8
};

typedef struct {
    double model_cost;
    double standalone_savings;
    double standalone_score;
    double exclusive_savings;
    double exclusive_score;
    double report_score;
    int64_t exclusive_bases;
    int exclusive_instances;
    CandidateAcceptState accept_state;
    CandidateQualityTier quality_tier;
    uint32_t quality_flags;
} CandidateMdlState;

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
    mdl_uid_t    id;
    char    *consensus;       /* numeric bases (0/1/2/3) */
    int      consensus_length;
    int      component_id;
    int      topology;        /* TOPO_LINEAR, TOPO_COMPLEX, TOPO_CYCLIC */
    freq_t   estimated_copies; /* estimated copy number */

    Instance *instances;
    int       num_instances;
    int       cap_instances;

    /* Legacy MDL aliases.  During the score-state migration, keep these
     * synchronized with mdl.report_score/model_cost for existing callers. */
    double   mdl_score;
    double   model_cost;
    uint32_t discovery_flags;
    CandidateMdlState mdl;
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

static inline int candidate_is_accepted(const CandidateFamily *f)
{
    if (!f) return 0;
    switch (f->mdl.accept_state) {
    case CAND_ACCEPT_EXCLUSIVE:
    case CAND_ACCEPT_STANDALONE:
        return 1;
    case CAND_ACCEPT_REJECTED:
    case CAND_ACCEPT_PRUNED:
        return 0;
    case CAND_ACCEPT_UNSCORED:
    default:
        /* Backward compatibility for unit tests and pre-selection helpers
         * that still mark accepted records by a positive legacy score. */
        return f->mdl_score > 0.0;
    }
}

static inline double candidate_report_score(const CandidateFamily *f)
{
    if (!f) return 0.0;
    return (f->mdl.report_score != 0.0) ? f->mdl.report_score : f->mdl_score;
}

static inline double candidate_model_cost(const CandidateFamily *f)
{
    if (!f) return 0.0;
    return (f->mdl.model_cost != 0.0) ? f->mdl.model_cost : f->model_cost;
}

static inline const char *candidate_accept_state_name(CandidateAcceptState s)
{
    switch (s) {
    case CAND_ACCEPT_UNSCORED:   return "unscored";
    case CAND_ACCEPT_REJECTED:   return "rejected";
    case CAND_ACCEPT_EXCLUSIVE:  return "exclusive";
    case CAND_ACCEPT_STANDALONE: return "standalone";
    case CAND_ACCEPT_PRUNED:     return "pruned";
    default:                     return "unknown";
    }
}

static inline const char *candidate_quality_tier_name(CandidateQualityTier t)
{
    switch (t) {
    case CAND_TIER_UNSET:     return "unset";
    case CAND_TIER_CORE:      return "core";
    case CAND_TIER_SUPPORTED: return "supported";
    case CAND_TIER_RESCUE:    return "rescue";
    case CAND_TIER_WARN:      return "warn";
    case CAND_TIER_REJECT:    return "reject";
    default:                  return "unknown";
    }
}

static inline const char *candidate_quality_flags_string(uint32_t flags,
                                                         char *buf,
                                                         size_t buf_size)
{
    if (!buf || buf_size == 0) return "";
    buf[0] = '\0';
    if (flags == CAND_QF_NONE) {
        snprintf(buf, buf_size, "none");
        return buf;
    }

    #define CAND_APPEND_FLAG(mask, name) do {                              \
        if (flags & (mask)) {                                              \
            size_t _len = strlen(buf);                                     \
            if (_len + 1 < buf_size && buf[0] != '\0') {                   \
                strncat(buf, "|", buf_size - _len - 1);                   \
                _len = strlen(buf);                                       \
            }                                                              \
            if (_len + 1 < buf_size)                                      \
                strncat(buf, (name), buf_size - _len - 1);                \
        }                                                                  \
    } while (0)

    CAND_APPEND_FLAG(CAND_QF_SHORT_CONSENSUS, "short_consensus");
    CAND_APPEND_FLAG(CAND_QF_LOW_COPY, "low_copy");
    CAND_APPEND_FLAG(CAND_QF_HIGH_DIVERGENCE, "high_divergence");
    CAND_APPEND_FLAG(CAND_QF_STANDALONE_FALLBACK, "standalone_fallback");
    CAND_APPEND_FLAG(CAND_QF_NO_EXCLUSIVE_BASES, "no_exclusive_bases");
    CAND_APPEND_FLAG(CAND_QF_LOW_EXCLUSIVE_FRAC, "low_exclusive_fraction");
    CAND_APPEND_FLAG(CAND_QF_NONPOSITIVE_EXCL_MDL, "nonpositive_exclusive_mdl");
    CAND_APPEND_FLAG(CAND_QF_PRUNED_REDUNDANT, "pruned_redundant");
    CAND_APPEND_FLAG(CAND_QF_RESCUE_DISCOVERY, "rescue_discovery");

    #undef CAND_APPEND_FLAG
    return buf;
}

#endif /* MDL_CANDIDATES_H */
