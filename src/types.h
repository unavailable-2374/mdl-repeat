#ifndef MDL_TYPES_H
#define MDL_TYPES_H

#include <stdint.h>
#include <stddef.h>
#include <inttypes.h>

/* --- Genome position and length types (human genome > INT_MAX) --- */
typedef int64_t gpos_t;    /* genome position */
typedef int64_t glen_t;    /* genome/sequence length */
typedef int32_t freq_t;    /* k-mer frequency (max ~1M for Alu, fits int32) */
typedef int32_t uid_t;     /* unitig/family ID */

/* --- DNA encoding (matches RepeatScout convention) --- */
#define DNA_A   0
#define DNA_C   1
#define DNA_G   2
#define DNA_T   3
#define DNA_N   99

/* --- Constants --- */
#define PADLENGTH       11000   /* must be >= max extension distance */
#define DEFAULT_TANDEMDIST  500
#define DEFAULT_MINTHRESH   3
#define DEFAULT_MAXN        10000   /* max instances per candidate */

/* --- Inline DNA helpers (carried from RepeatScout) --- */
static inline char char_to_num(char c)
{
    switch (c) {
        case 'A': case 'a': return DNA_A;
        case 'C': case 'c': return DNA_C;
        case 'G': case 'g': return DNA_G;
        case 'T': case 't': return DNA_T;
        default:            return DNA_N;
    }
}

static inline char num_to_char(char z)
{
    switch (z) {
        case DNA_A: return 'A';
        case DNA_C: return 'C';
        case DNA_G: return 'G';
        case DNA_T: return 'T';
        default:    return 'N';
    }
}

static inline char dna_complement(char c)
{
    if (c == DNA_N) return DNA_N;
    return 3 - c;
}

#endif /* MDL_TYPES_H */
