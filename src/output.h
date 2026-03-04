#ifndef MDL_OUTPUT_H
#define MDL_OUTPUT_H

#include "types.h"
#include "candidates.h"
#include "mdl.h"
#include "genome.h"

/*
 * Write accepted repeat families to FASTA.
 * Header: >R=<id> length=<len> copies=<n> mdl=<score>
 * Only writes families with mdl_score > 0.
 */
int output_fasta(const char *filename, const CandidateList *cl);

/*
 * Write repeat instances to BED6 format.
 * Columns: chr start end name score strand
 * Only writes instances of accepted families.
 */
int output_bed(const char *filename, const CandidateList *cl,
               const Genome *genome);

/*
 * Write per-family statistics to TSV.
 * Columns: id, consensus_length, num_instances, divergence_mean,
 *          mdl_score, model_cost, topology
 */
int output_stats(const char *filename, const CandidateList *cl);

#endif /* MDL_OUTPUT_H */
