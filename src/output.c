#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "output.h"

int output_fasta(const char *filename, const CandidateList *cl)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "output_fasta: could not open '%s'\n", filename);
        return -1;
    }

    int written = 0;
    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *f = &cl->families[i];
        if (f->mdl_score <= 0) continue;

        fprintf(fp, ">R=%d length=%d copies=%d mdl=%.1f\n",
                f->id, f->consensus_length, f->num_instances, f->mdl_score);

        for (int j = 0; j < f->consensus_length; j++) {
            fputc(num_to_char(f->consensus[j]), fp);
            if ((j + 1) % 80 == 0) fputc('\n', fp);
        }
        if (f->consensus_length % 80 != 0) fputc('\n', fp);

        written++;
    }

    fclose(fp);
    fprintf(stderr, "Wrote %d families to %s\n", written, filename);
    return written;
}

int output_bed(const char *filename, const CandidateList *cl,
               const Genome *genome)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "output_bed: could not open '%s'\n", filename);
        return -1;
    }

    int written = 0;
    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *f = &cl->families[i];
        if (f->mdl_score <= 0) continue;

        for (int j = 0; j < f->num_instances; j++) {
            Instance *inst = &f->instances[j];
            gpos_t raw_start = inst->position - PADLENGTH;
            gpos_t raw_end = raw_start + (glen_t)inst->aligned_length;

            /* Get chromosome name */
            const char *chrom = "unknown";
            if (inst->seq_index >= 0 && inst->seq_index < genome->num_sequences &&
                genome->sequence_ids && genome->sequence_ids[inst->seq_index]) {
                chrom = genome->sequence_ids[inst->seq_index];
            }

            /* Adjust to local chromosome coordinates */
            gpos_t chr_offset = 0;
            if (inst->seq_index > 0 && inst->seq_index < genome->num_sequences) {
                chr_offset = genome->boundaries[inst->seq_index - 1];
            }

            gpos_t local_start = raw_start - chr_offset;
            gpos_t local_end = raw_end - chr_offset;

            int score = (int)(1000 * (1.0 - inst->divergence));
            if (score < 0) score = 0;
            if (score > 1000) score = 1000;

            fprintf(fp, "%s\t%" PRId64 "\t%" PRId64 "\tR=%d\t%d\t%c\n",
                    chrom,
                    (int64_t)local_start, (int64_t)local_end,
                    f->id, score,
                    inst->strand > 0 ? '+' : '-');
            written++;
        }
    }

    fclose(fp);
    fprintf(stderr, "Wrote %d instances to %s\n", written, filename);
    return written;
}

int output_stats(const char *filename, const CandidateList *cl)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "output_stats: could not open '%s'\n", filename);
        return -1;
    }

    /* Header */
    fprintf(fp, "family_id\tconsensus_length\tnum_instances\t"
                "divergence_mean\tmdl_score\tmodel_cost\ttopology\n");

    int written = 0;
    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *f = &cl->families[i];

        float avg_div = 0;
        for (int j = 0; j < f->num_instances; j++)
            avg_div += f->instances[j].divergence;
        if (f->num_instances > 0) avg_div /= f->num_instances;

        const char *topo = (f->topology == TOPO_LINEAR) ? "linear" :
                           (f->topology == TOPO_CYCLIC) ? "cyclic" : "complex";

        fprintf(fp, "%d\t%d\t%d\t%.4f\t%.1f\t%.1f\t%s\n",
                f->id, f->consensus_length, f->num_instances,
                avg_div, f->mdl_score, f->model_cost, topo);
        written++;
    }

    fclose(fp);
    fprintf(stderr, "Wrote %d family stats to %s\n", written, filename);
    return written;
}
