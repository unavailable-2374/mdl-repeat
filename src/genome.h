#ifndef MDL_GENOME_H
#define MDL_GENOME_H

#include "types.h"

typedef struct {
    char    *sequence;       /* numeric encoding (0/1/2/3/99), padded at start */
    glen_t   length;         /* total length including PADLENGTH padding */
    glen_t   raw_length;     /* actual sequence bases (no padding) */
    gpos_t  *boundaries;     /* boundary positions for each sequence (pre-padding) */
    int      num_sequences;  /* number of FASTA records */
    char   **sequence_ids;   /* identifier for each FASTA record */
} Genome;

/*
 * Load a multi-FASTA file into a Genome struct.
 * Sequence is padded with DNA_N at the front (PADLENGTH bases).
 * Returns NULL on failure.
 */
Genome *genome_load(const char *filename);

/*
 * Free all memory associated with a Genome.
 */
void genome_free(Genome *g);

/*
 * Print summary statistics for a loaded genome.
 */
void genome_print_stats(const Genome *g);

/*
 * Check if a position crosses a sequence boundary.
 * Returns 1 if the range [pos, pos+len) stays within one sequence, 0 otherwise.
 */
int genome_check_boundary(const Genome *g, gpos_t pos, glen_t len);

/*
 * Get the sequence index that contains position pos.
 * Returns -1 if pos is in the padding region.
 */
int genome_get_seq_index(const Genome *g, gpos_t pos);

#endif /* MDL_GENOME_H */
