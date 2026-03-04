#ifndef MDL_KMER_H
#define MDL_KMER_H

#include "types.h"
#include "genome.h"

/*
 * K-mer hash table entry.
 * After kmer_build_positions(): stores all occurrence positions for
 * efficient index-based lookup (avoids per-family genome scans).
 */
typedef struct kmer_entry {
    uint64_t kmer;              /* canonical k-mer, packed 2 bits/base */
    freq_t   frequency;
    gpos_t   last_plus_occ;     /* last forward-strand occurrence position */
    gpos_t   last_minus_occ;    /* last reverse-strand occurrence position */
    gpos_t  *positions;         /* all positions (sign encodes strand:
                                   positive = forward, negative = RC of canonical) */
    int32_t  num_positions;     /* number of stored positions */
    int32_t  cap_positions;     /* allocated capacity */
    struct kmer_entry *next;    /* hash chain */
} KmerEntry;

#define KMER_MAX_POSITIONS  50000  /* max positions stored per k-mer */

typedef struct {
    KmerEntry **buckets;
    int         k;
    size_t      table_size;     /* prime, scaled to genome size */
    int64_t     num_entries;    /* total distinct k-mers stored */
    int64_t     num_frequent;   /* entries with freq >= minthresh */
} KmerTable;

/* --- Packed k-mer operations --- */

/* Pack numeric sequence[0..k-1] into a uint64_t (2 bits/base).
 * Returns UINT64_MAX if any base is DNA_N. Max k=31. */
uint64_t kmer_pack(const char *seq, int k);

/* Reverse complement of a packed k-mer */
uint64_t kmer_revcomp(uint64_t kmer, int k);

/* Canonical form = min(kmer, revcomp) */
uint64_t kmer_canonical(uint64_t kmer, int k);

/* --- Hash table operations --- */

/* Count all k-mers in genome with TANDEMDIST filtering.
 * Returns a new KmerTable, or NULL on failure. */
KmerTable *kmer_count(const Genome *g, int k, int tandemdist);

/* Remove entries with frequency < min_freq. Returns number removed. */
int64_t kmer_trim(KmerTable *kt, freq_t min_freq);

/* Build position index: second genome pass to collect ALL occurrence
 * positions for surviving (post-trim) k-mers.
 * Must be called after kmer_trim() and before candidates_extract().
 * Replaces per-family O(N) genome scans with O(freq) lookups.
 * num_threads > 1 enables parallel genome scan with atomic position append. */
void kmer_build_positions(KmerTable *kt, const Genome *g, int num_threads);

/* Look up a packed canonical k-mer. Returns entry or NULL. */
KmerEntry *kmer_lookup(const KmerTable *kt, uint64_t canonical_kmer);

/* Get total number of distinct k-mers (including trimmed) */
int64_t kmer_count_entries(const KmerTable *kt);

/* Free the k-mer table */
void kmer_free(KmerTable *kt);

/* Print frequency distribution summary */
void kmer_print_stats(const KmerTable *kt);

#endif /* MDL_KMER_H */
