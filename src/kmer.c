#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "kmer.h"

/* --- KmerEntry memory pool --- */

#define POOL_BLOCK_SIZE 4096

typedef struct kmer_pool_block {
    KmerEntry entries[POOL_BLOCK_SIZE];
    struct kmer_pool_block *next;
} KmerPoolBlock;

typedef struct {
    KmerPoolBlock *head;       /* linked list of blocks */
    int            next_idx;   /* next free index in current block */
} KmerPool;

static KmerPool g_kmer_pool = { NULL, POOL_BLOCK_SIZE };

static KmerEntry *pool_alloc(void)
{
    if (g_kmer_pool.next_idx >= POOL_BLOCK_SIZE) {
        KmerPoolBlock *blk = malloc(sizeof(KmerPoolBlock));
        if (!blk) return NULL;
        blk->next = g_kmer_pool.head;
        g_kmer_pool.head = blk;
        g_kmer_pool.next_idx = 0;
    }
    return &g_kmer_pool.head->entries[g_kmer_pool.next_idx++];
}

static void pool_free_all(void)
{
    KmerPoolBlock *blk = g_kmer_pool.head;
    while (blk) {
        KmerPoolBlock *next = blk->next;
        free(blk);
        blk = next;
    }
    g_kmer_pool.head = NULL;
    g_kmer_pool.next_idx = POOL_BLOCK_SIZE;
}

/* --- Primes for table sizing --- */

/* Find the smallest prime >= n (simple trial division, called once) */
static size_t next_prime(size_t n)
{
    if (n <= 2) return 2;
    if (n % 2 == 0) n++;
    for (;; n += 2) {
        int is_prime = 1;
        for (size_t d = 3; d * d <= n; d += 2) {
            if (n % d == 0) { is_prime = 0; break; }
        }
        if (is_prime) return n;
    }
}

/* --- Packed k-mer operations --- */

uint64_t kmer_pack(const char *seq, int k)
{
    uint64_t packed = 0;
    for (int i = 0; i < k; i++) {
        if (seq[i] == DNA_N)
            return UINT64_MAX;
        packed = (packed << 2) | ((uint64_t)(seq[i] & 3));
    }
    return packed;
}

uint64_t kmer_revcomp(uint64_t kmer, int k)
{
    uint64_t rc = 0;
    for (int i = 0; i < k; i++) {
        rc = (rc << 2) | (3ULL - (kmer & 3ULL));
        kmer >>= 2;
    }
    return rc;
}

uint64_t kmer_canonical(uint64_t kmer, int k)
{
    uint64_t rc = kmer_revcomp(kmer, k);
    return (kmer <= rc) ? kmer : rc;
}

/* --- Hash function --- */

static inline size_t kmer_hash(uint64_t kmer, size_t table_size)
{
    /* Fibonacci hashing: good mixing for sequential k-mers */
    uint64_t h = kmer * 11400714819323198485ULL; /* 2^64 / phi */
    return (size_t)(h % table_size);
}

/* --- Hash table operations --- */

KmerTable *kmer_count(const Genome *g, int k, int tandemdist)
{
    if (k > 31) {
        fprintf(stderr, "kmer_count: k=%d exceeds maximum 31 (64-bit packing limit)\n", k);
        return NULL;
    }

    KmerTable *kt = calloc(1, sizeof(KmerTable));
    if (kt == NULL) {
        fprintf(stderr, "kmer_count: out of memory for KmerTable\n");
        return NULL;
    }

    kt->k = k;
    /* Table size: genome_length/4 rounded to next prime.
     * Larger than RepeatScout's fixed 16M — avoids chaining for large genomes. */
    kt->table_size = next_prime((size_t)(g->length / 4));
    if (kt->table_size < 16000057)
        kt->table_size = 16000057; /* minimum for small genomes */

    kt->buckets = calloc(kt->table_size, sizeof(KmerEntry *));
    if (kt->buckets == NULL) {
        fprintf(stderr, "kmer_count: could not allocate %" PRId64 " buckets\n",
                (int64_t)kt->table_size);
        free(kt);
        return NULL;
    }

    /* Single-pass counting with TANDEMDIST filtering */
    glen_t end = g->length - k + 1;

    for (glen_t i = PADLENGTH; i < end; i++) {
        /* Progress indicator */
        if (i % 500000 == 0)
            fprintf(stderr, "  kmer counting: %" PRId64 " / %" PRId64 " (%.0f%%)\r",
                    (int64_t)(i - PADLENGTH), (int64_t)(end - PADLENGTH),
                    100.0 * (double)(i - PADLENGTH) / (double)(end - PADLENGTH));

        /* Pack the k-mer at position i */
        uint64_t fwd = kmer_pack(g->sequence + i, k);
        if (fwd == UINT64_MAX)
            continue; /* contains N */

        uint64_t rc = kmer_revcomp(fwd, k);
        uint64_t canon = (fwd <= rc) ? fwd : rc;
        int is_reverse = (fwd > rc);  /* this occurrence is reverse-strand */

        size_t h = kmer_hash(canon, kt->table_size);

        /* Search for existing entry */
        KmerEntry *entry = kt->buckets[h];
        while (entry != NULL) {
            if (entry->kmer == canon)
                break;
            entry = entry->next;
        }

        if (entry != NULL) {
            /* Existing k-mer: apply TANDEMDIST filtering */
            if (is_reverse) {
                if ((gpos_t)i - entry->last_minus_occ >= tandemdist)
                    entry->frequency++;
                entry->last_minus_occ = (gpos_t)i;
            } else {
                if ((gpos_t)i - entry->last_plus_occ >= tandemdist)
                    entry->frequency++;
                entry->last_plus_occ = (gpos_t)i;
            }
        } else {
            /* New k-mer (pool-allocated) */
            entry = pool_alloc();
            if (entry == NULL) {
                fprintf(stderr, "kmer_count: out of memory at position %" PRId64 "\n",
                        (int64_t)i);
                kmer_free(kt);
                return NULL;
            }
            entry->kmer = canon;
            entry->frequency = 1;
            if (is_reverse) {
                entry->last_plus_occ = -1000000;
                entry->last_minus_occ = (gpos_t)i;
            } else {
                entry->last_plus_occ = (gpos_t)i;
                entry->last_minus_occ = -1000000;
            }
            entry->positions = NULL;
            entry->num_positions = 0;
            entry->cap_positions = 0;
            entry->next = kt->buckets[h];
            kt->buckets[h] = entry;
            kt->num_entries++;
        }
    }

    fprintf(stderr, "  kmer counting: done. %" PRId64 " distinct k-mers\n",
            kt->num_entries);

    return kt;
}

int64_t kmer_trim(KmerTable *kt, freq_t min_freq)
{
    int64_t removed = 0;

    for (size_t h = 0; h < kt->table_size; h++) {
        KmerEntry *prev = NULL;
        KmerEntry *entry = kt->buckets[h];
        while (entry != NULL) {
            if (entry->frequency < min_freq) {
                KmerEntry *to_unlink = entry;
                entry = entry->next;
                if (prev == NULL)
                    kt->buckets[h] = entry;
                else
                    prev->next = entry;
                /* Free only the positions array; entry itself is pool-managed */
                free(to_unlink->positions);
                to_unlink->positions = NULL;
                kt->num_entries--;
                removed++;
            } else {
                prev = entry;
                entry = entry->next;
            }
        }
    }

    kt->num_frequent = kt->num_entries;
    return removed;
}

KmerEntry *kmer_lookup(const KmerTable *kt, uint64_t canonical_kmer)
{
    size_t h = kmer_hash(canonical_kmer, kt->table_size);
    KmerEntry *entry = kt->buckets[h];
    while (entry != NULL) {
        if (entry->kmer == canonical_kmer)
            return entry;
        entry = entry->next;
    }
    return NULL;
}

int64_t kmer_count_entries(const KmerTable *kt)
{
    return kt->num_entries;
}

void kmer_free(KmerTable *kt)
{
    if (kt == NULL) return;
    if (kt->buckets) {
        /* Free positions arrays (allocated individually) */
        for (size_t h = 0; h < kt->table_size; h++) {
            KmerEntry *entry = kt->buckets[h];
            while (entry != NULL) {
                free(entry->positions);
                entry = entry->next;
            }
        }
        free(kt->buckets);
    }
    /* Free all KmerEntry structs in bulk via pool */
    pool_free_all();
    free(kt);
}

/* --- Parallel position index building --- */

typedef struct {
    const KmerTable *kt;
    const Genome *g;
    glen_t chunk_start;
    glen_t chunk_end;
    int64_t positions_stored;
} PosWorkerArgs;

static void *pos_worker(void *arg)
{
    PosWorkerArgs *a = (PosWorkerArgs *)arg;
    int k = a->kt->k;
    a->positions_stored = 0;

    for (glen_t i = a->chunk_start; i < a->chunk_end; i++) {
        uint64_t fwd = kmer_pack(a->g->sequence + i, k);
        if (fwd == UINT64_MAX) continue;

        uint64_t rc = kmer_revcomp(fwd, k);
        uint64_t canon = (fwd <= rc) ? fwd : rc;
        int is_reverse = (fwd > rc);

        KmerEntry *entry = kmer_lookup(a->kt, canon);
        if (!entry) continue;

        /* Atomic claim a slot; silently drop if over capacity */
        int32_t idx = __atomic_fetch_add(&entry->num_positions, 1,
                                         __ATOMIC_RELAXED);
        if (idx < entry->cap_positions) {
            entry->positions[idx] = is_reverse ? -(gpos_t)i : (gpos_t)i;
            a->positions_stored++;
        }
    }
    return NULL;
}

void kmer_build_positions(KmerTable *kt, const Genome *g, int num_threads)
{
    int k = kt->k;
    glen_t end = g->length - k + 1;

    fprintf(stderr, "  Building k-mer position index...\n");

    /* Pre-allocate position arrays for all surviving entries.
     * This enables lock-free parallel filling via atomic slot claiming. */
    for (size_t h = 0; h < kt->table_size; h++) {
        for (KmerEntry *e = kt->buckets[h]; e != NULL; e = e->next) {
            int cap = (int)e->frequency * 6;
            if (cap < 64) cap = 64;
            if (cap > KMER_MAX_POSITIONS) cap = KMER_MAX_POSITIONS;
            e->positions = malloc((size_t)cap * sizeof(gpos_t));
            if (!e->positions) cap = 0;
            e->cap_positions = cap;
            e->num_positions = 0;
        }
    }

    int64_t total_positions = 0;

    if (num_threads <= 1) {
        /* Sequential path */
        for (glen_t i = PADLENGTH; i < end; i++) {
            if (i % 2000000 == 0)
                fprintf(stderr, "  position index: %" PRId64 " / %" PRId64
                        " (%.0f%%)\r",
                        (int64_t)(i - PADLENGTH), (int64_t)(end - PADLENGTH),
                        100.0 * (double)(i - PADLENGTH) /
                        (double)(end - PADLENGTH));

            uint64_t fwd = kmer_pack(g->sequence + i, k);
            if (fwd == UINT64_MAX) continue;

            uint64_t rc = kmer_revcomp(fwd, k);
            uint64_t canon = (fwd <= rc) ? fwd : rc;
            int is_reverse = (fwd > rc);

            KmerEntry *entry = kmer_lookup(kt, canon);
            if (!entry) continue;

            if (entry->num_positions >= entry->cap_positions) continue;

            entry->positions[entry->num_positions++] =
                is_reverse ? -(gpos_t)i : (gpos_t)i;
            total_positions++;
        }
    } else {
        /* Parallel path: divide genome into chunks */
        glen_t range = end - PADLENGTH;
        PosWorkerArgs *args = malloc((size_t)num_threads * sizeof(PosWorkerArgs));
        pthread_t *threads = malloc((size_t)num_threads * sizeof(pthread_t));
        if (!args || !threads) {
            free(args); free(threads);
            fprintf(stderr, "WARNING: alloc failed, falling back to single-threaded\n");
            kmer_build_positions(kt, g, 1);
            return;
        }

        for (int t = 0; t < num_threads; t++) {
            args[t].kt = kt;
            args[t].g = g;
            args[t].chunk_start = PADLENGTH + range * t / num_threads;
            args[t].chunk_end = PADLENGTH + range * (t + 1) / num_threads;
            args[t].positions_stored = 0;
        }

        int launched = 0;
        for (int t = 0; t < num_threads; t++) {
            if (pthread_create(&threads[t], NULL, pos_worker, &args[t]) != 0)
                break;
            launched++;
        }
        /* If some threads failed, process remaining chunks on main thread */
        for (int t = launched; t < num_threads; t++) {
            args[t].chunk_start = args[t].chunk_start; /* already set */
            pos_worker(&args[t]);
        }
        for (int t = 0; t < launched; t++)
            pthread_join(threads[t], NULL);

        for (int t = 0; t < num_threads; t++)
            total_positions += args[t].positions_stored;

        free(args);
        free(threads);

        /* Clamp num_positions that may have overshot cap_positions */
        for (size_t h = 0; h < kt->table_size; h++) {
            for (KmerEntry *e = kt->buckets[h]; e != NULL; e = e->next) {
                if (e->num_positions > e->cap_positions)
                    e->num_positions = e->cap_positions;
            }
        }
    }

    fprintf(stderr, "  Position index: %" PRId64 " positions stored (%.0f MB)       \n",
            total_positions,
            (double)total_positions * (double)sizeof(gpos_t) / (1024.0 * 1024.0));
}

void kmer_print_stats(const KmerTable *kt)
{
    /* Compute frequency distribution: bins for freq 1,2,3,...,9, 10-99, 100-999, 1000+ */
    int64_t bins[13] = {0}; /* [0]=freq1, [1]=freq2, ..., [8]=freq9, [9]=10-99, [10]=100-999, [11]=1000-9999, [12]=10000+ */
    freq_t max_freq = 0;
    int64_t total_entries = 0;

    for (size_t h = 0; h < kt->table_size; h++) {
        for (KmerEntry *e = kt->buckets[h]; e != NULL; e = e->next) {
            total_entries++;
            if (e->frequency > max_freq)
                max_freq = e->frequency;
            if (e->frequency <= 9)
                bins[e->frequency - 1]++;
            else if (e->frequency <= 99)
                bins[9]++;
            else if (e->frequency <= 999)
                bins[10]++;
            else if (e->frequency <= 9999)
                bins[11]++;
            else
                bins[12]++;
        }
    }

    fprintf(stderr, "K-mer table statistics (k=%d):\n", kt->k);
    fprintf(stderr, "  Total distinct k-mers: %" PRId64 "\n", total_entries);
    fprintf(stderr, "  Max frequency:         %d\n", (int)max_freq);
    fprintf(stderr, "  Frequency distribution:\n");
    for (int i = 0; i < 9; i++) {
        if (bins[i] > 0)
            fprintf(stderr, "    freq=%d: %" PRId64 "\n", i + 1, bins[i]);
    }
    if (bins[9] > 0)  fprintf(stderr, "    freq=10-99:    %" PRId64 "\n", bins[9]);
    if (bins[10] > 0) fprintf(stderr, "    freq=100-999:  %" PRId64 "\n", bins[10]);
    if (bins[11] > 0) fprintf(stderr, "    freq=1000-9999:%" PRId64 "\n", bins[11]);
    if (bins[12] > 0) fprintf(stderr, "    freq=10000+:   %" PRId64 "\n", bins[12]);
}
