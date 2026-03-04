#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "genome.h"

Genome *genome_load(const char *filename)
{
    FILE *fp;
    glen_t file_size;
    Genome *g;
    glen_t i;
    int seq;
    int j;
    char c;
    char id_buffer[1024];
    int boundaries_cap = 100;

    /* Open file and get size estimate */
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "genome_load: could not open '%s'\n", filename);
        return NULL;
    }
    fseek(fp, 0, SEEK_END);
    file_size = (glen_t)ftell(fp);
    fclose(fp);

    /* Allocate genome struct */
    g = calloc(1, sizeof(Genome));
    if (g == NULL) {
        fprintf(stderr, "genome_load: out of memory for Genome struct\n");
        return NULL;
    }

    /* Allocate sequence buffer: file_size is an overestimate of bases */
    g->sequence = malloc((size_t)(file_size + PADLENGTH + 1));
    if (g->sequence == NULL) {
        fprintf(stderr, "genome_load: could not allocate %" PRId64 " bytes for sequence\n",
                (int64_t)(file_size + PADLENGTH + 1));
        free(g);
        return NULL;
    }

    /* Allocate boundaries */
    g->boundaries = malloc((size_t)boundaries_cap * sizeof(gpos_t));
    if (g->boundaries == NULL) {
        fprintf(stderr, "genome_load: could not allocate boundaries\n");
        free(g->sequence);
        free(g);
        return NULL;
    }

    /* Pad the beginning with N */
    for (i = 0; i < PADLENGTH; i++)
        g->sequence[i] = DNA_N;

    /* Re-open for reading */
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "genome_load: could not reopen '%s'\n", filename);
        genome_free(g);
        return NULL;
    }

    i = 0;   /* bases read (not counting padding) */
    seq = 0; /* number of FASTA records seen */

    while (!feof(fp)) {
        c = (char)getc(fp);
        if (c == EOF) continue;
        if (c == '\n') continue;

        if (c == '>') {
            /* Parse sequence ID (up to first whitespace) */
            j = 0;
            while ((c = (char)getc(fp)) != '\n' && c != EOF) {
                if (c == ' ' || c == '\t') {
                    /* Skip rest of header line */
                    while ((c = (char)getc(fp)) != '\n' && c != EOF)
                        ;
                    break;
                }
                if (j < (int)sizeof(id_buffer) - 1)
                    id_buffer[j++] = c;
            }
            id_buffer[j] = '\0';

            /* Store sequence ID */
            char **tmp_ids = realloc(g->sequence_ids,
                                     (size_t)(seq + 2) * sizeof(char *));
            if (tmp_ids == NULL) {
                fprintf(stderr, "genome_load: out of memory for sequence IDs\n");
                fclose(fp);
                genome_free(g);
                return NULL;
            }
            g->sequence_ids = tmp_ids;
            g->sequence_ids[seq] = malloc(strlen(id_buffer) + 1);
            if (g->sequence_ids[seq] == NULL) {
                fprintf(stderr, "genome_load: out of memory for ID string\n");
                fclose(fp);
                genome_free(g);
                return NULL;
            }
            strcpy(g->sequence_ids[seq], id_buffer);
            g->sequence_ids[seq + 1] = NULL;

            /* Store boundary (start of this sequence in base coords) */
            if (seq > 0) {
                if (seq >= boundaries_cap) {
                    boundaries_cap *= 2;
                    gpos_t *tmp_b = realloc(g->boundaries,
                                            (size_t)boundaries_cap * sizeof(gpos_t));
                    if (tmp_b == NULL) {
                        fprintf(stderr, "genome_load: out of memory for boundaries\n");
                        fclose(fp);
                        genome_free(g);
                        return NULL;
                    }
                    g->boundaries = tmp_b;
                }
                g->boundaries[seq - 1] = i; /* pre-padding position */
            }
            seq++;
        } else {
            /* Read sequence data */
            if (c > 64) {
                g->sequence[i + PADLENGTH] = char_to_num(c);
                i++;
            }
            while ((c = (char)getc(fp)) != '\n' && !feof(fp)) {
                if (c > 64) {
                    g->sequence[i + PADLENGTH] = char_to_num(c);
                    i++;
                }
            }
        }
    }

    fclose(fp);

    /* Ensure at least 1 sequence recorded */
    if (seq == 0)
        seq = 1;

    /* Store final boundary sentinel */
    if (seq + 1 >= boundaries_cap) {
        boundaries_cap = seq + 2;
        gpos_t *tmp_b = realloc(g->boundaries,
                                (size_t)boundaries_cap * sizeof(gpos_t));
        if (tmp_b == NULL) {
            fprintf(stderr, "genome_load: out of memory for final boundary\n");
            genome_free(g);
            return NULL;
        }
        g->boundaries = tmp_b;
    }
    g->boundaries[seq - 1] = i + 1;  /* end sentinel */
    g->boundaries[seq] = 0;          /* zero terminator */

    g->num_sequences = seq;
    g->raw_length = i;
    g->length = i + PADLENGTH;

    return g;
}

void genome_free(Genome *g)
{
    if (g == NULL) return;
    free(g->sequence);
    free(g->boundaries);
    if (g->sequence_ids) {
        for (int i = 0; i < g->num_sequences; i++)
            free(g->sequence_ids[i]);
        free(g->sequence_ids);
    }
    free(g);
}

void genome_print_stats(const Genome *g)
{
    if (g == NULL) {
        fprintf(stderr, "genome_print_stats: NULL genome\n");
        return;
    }
    fprintf(stderr, "Genome statistics:\n");
    fprintf(stderr, "  Total length (with padding): %" PRId64 "\n", (int64_t)g->length);
    fprintf(stderr, "  Raw sequence bases:          %" PRId64 "\n", (int64_t)g->raw_length);
    fprintf(stderr, "  Number of sequences:         %d\n", g->num_sequences);
    for (int i = 0; i < g->num_sequences; i++) {
        gpos_t start = (i == 0) ? 0 : g->boundaries[i - 1];
        gpos_t end = g->boundaries[i];
        fprintf(stderr, "  [%d] %s: %" PRId64 " bp\n",
                i,
                (g->sequence_ids && g->sequence_ids[i]) ? g->sequence_ids[i] : "unnamed",
                (int64_t)(end - start - (i == g->num_sequences - 1 ? 1 : 0)));
    }
}

int genome_check_boundary(const Genome *g, gpos_t pos, glen_t len)
{
    /* pos is in padded coordinates */
    gpos_t raw_pos = pos - PADLENGTH;
    gpos_t raw_end = raw_pos + len;

    if (raw_pos < 0 || raw_end > g->raw_length)
        return 0;

    /* Find which sequence raw_pos belongs to */
    int seq_start = 0;
    for (int i = 0; i < g->num_sequences; i++) {
        gpos_t seq_end = g->boundaries[i];
        if (i < g->num_sequences - 1) {
            /* boundaries[i] stores start of sequence i+1 */
            if (raw_pos >= seq_start && raw_pos < seq_end) {
                /* Check if end is also within this sequence */
                return (raw_end <= seq_end) ? 1 : 0;
            }
            seq_start = seq_end;
        } else {
            /* Last sequence */
            return (raw_pos >= seq_start && raw_end <= seq_end) ? 1 : 0;
        }
    }
    return 0;
}

int genome_get_seq_index(const Genome *g, gpos_t pos)
{
    gpos_t raw_pos = pos - PADLENGTH;
    if (raw_pos < 0)
        return -1;

    for (int i = 0; i < g->num_sequences; i++) {
        if (i == 0) {
            if (raw_pos < g->boundaries[0])
                return 0;
        } else {
            if (raw_pos >= g->boundaries[i - 1] && raw_pos < g->boundaries[i])
                return i;
        }
    }
    return -1;
}
