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
        if (!candidate_is_accepted(f)) continue;

        /* Mean instance divergence — surfaces a per-family identity signal
         * downstream tools (Refiner_mdl) need without re-parsing the stats TSV. */
        float avg_div = 0.0f;
        for (int j = 0; j < f->num_instances; j++)
            avg_div += f->instances[j].divergence;
        if (f->num_instances > 0) avg_div /= (float)f->num_instances;

        const char *topo = (f->topology == TOPO_LINEAR) ? "linear" :
                           (f->topology == TOPO_CYCLIC) ? "cyclic" : "complex";
        char qflags[256];
        candidate_quality_flags_string(f->mdl.quality_flags, qflags,
                                       sizeof(qflags));

        /* Header schema (v6.1+): backward-compatible — older parsers that
         * only match `>R=N length=L copies=C mdl=M` still succeed; new
         * fields `div=` and `topo=` are appended. */
        fprintf(fp, ">R=%d length=%d copies=%d mdl=%.1f div=%.3f topo=%s "
                    "accept=%s tier=%s flags=0x%08x qflags=%s\n",
                f->id, f->consensus_length, f->num_instances,
                candidate_report_score(f), avg_div, topo,
                candidate_accept_state_name(f->mdl.accept_state),
                candidate_quality_tier_name(f->mdl.quality_tier),
                (unsigned)f->mdl.quality_flags, qflags);

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
        if (!candidate_is_accepted(f)) continue;

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
                "divergence_mean\tmdl_score\tmodel_cost\ttopology\t"
                "standalone_score\texclusive_score\texclusive_bases\t"
                "exclusive_instances\tacceptance\tquality_tier\tquality_flags\t"
                "quality_notes\n");

    int written = 0;
    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *f = &cl->families[i];

        float avg_div = 0;
        for (int j = 0; j < f->num_instances; j++)
            avg_div += f->instances[j].divergence;
        if (f->num_instances > 0) avg_div /= f->num_instances;

        const char *topo = (f->topology == TOPO_LINEAR) ? "linear" :
                           (f->topology == TOPO_CYCLIC) ? "cyclic" : "complex";
        char qflags[256];
        candidate_quality_flags_string(f->mdl.quality_flags, qflags,
                                       sizeof(qflags));

        fprintf(fp, "%d\t%d\t%d\t%.4f\t%.1f\t%.1f\t%s\t"
                    "%.1f\t%.1f\t%" PRId64 "\t%d\t%s\t%s\t0x%08x\t%s\n",
                f->id, f->consensus_length, f->num_instances,
                avg_div, candidate_report_score(f), candidate_model_cost(f), topo,
                f->mdl.standalone_score, f->mdl.exclusive_score,
                f->mdl.exclusive_bases, f->mdl.exclusive_instances,
                candidate_accept_state_name(f->mdl.accept_state),
                candidate_quality_tier_name(f->mdl.quality_tier),
                (unsigned)f->mdl.quality_flags, qflags);
        written++;
    }

    fclose(fp);
    fprintf(stderr, "Wrote %d family stats to %s\n", written, filename);
    return written;
}
