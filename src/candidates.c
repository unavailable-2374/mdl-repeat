/*
 * candidates.c — Candidate family utility functions
 *
 * Provides candidates_free() and candidates_print_stats().
 * The actual family discovery is in discover.c.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "candidates.h"

void candidates_free(CandidateList *cl)
{
    if (!cl) return;
    for (int i = 0; i < cl->num_families; i++) {
        free(cl->families[i].consensus);
        free(cl->families[i].instances);
    }
    free(cl->families);
    free(cl);
}

void candidates_print_stats(const CandidateList *cl)
{
    if (!cl) return;

    int max_len = 0, min_len = INT32_MAX;
    int max_copies = 0;
    int64_t total_instance_bases = 0;
    float min_div = 1.0f, max_div = 0.0f;

    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *f = &cl->families[i];
        if (f->consensus_length > max_len) max_len = f->consensus_length;
        if (f->consensus_length < min_len) min_len = f->consensus_length;
        if (f->num_instances > max_copies) max_copies = f->num_instances;

        for (int j = 0; j < f->num_instances; j++) {
            total_instance_bases += f->instances[j].aligned_length;
            if (f->instances[j].divergence < min_div)
                min_div = f->instances[j].divergence;
            if (f->instances[j].divergence > max_div)
                max_div = f->instances[j].divergence;
        }
    }

    fprintf(stderr, "Candidate statistics:\n");
    fprintf(stderr, "  Families:             %d\n", cl->num_families);
    fprintf(stderr, "  Consensus length:     min=%d, max=%d\n",
            cl->num_families > 0 ? min_len : 0, max_len);
    fprintf(stderr, "  Max copies:           %d\n", max_copies);
    fprintf(stderr, "  Total instance bases: %" PRId64 "\n", total_instance_bases);
    if (cl->num_families > 0 && total_instance_bases > 0)
        fprintf(stderr, "  Divergence range:     %.1f%% - %.1f%%\n",
                min_div * 100, max_div * 100);
}
