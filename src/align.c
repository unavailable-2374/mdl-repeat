#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <pthread.h>
#include <unistd.h>

#include "align.h"

/* Local strdup replacement to avoid _POSIX_C_SOURCE conflicts with mdl_uid_t */
static char *align_strdup(const char *s)
{
    size_t n = strlen(s) + 1;
    char *p = malloc(n);
    if (p) memcpy(p, s, n);
    return p;
}

/* Runtime-configurable refinement parameters (defaults match original #defines) */
int   g_align_gap            = -5;
int   g_align_maxoffset      = 12;
float g_align_max_divergence = 0.30f;

/* ================================================================
 * Internal constants
 * ================================================================ */

#define MAX_SEED_HITS    50000   /* cap per-family seed hits */
#define MAX_CONS_KMERS   10000   /* cap consensus k-mer set entries */
#define MAX_BAND_WIDTH   (2 * ALIGN_MAXOFFSET_LIMIT + 1)
#define ALIGN_MAX_EXTENSION 10000  /* max bases to extend per direction per iteration */
/* J' (v6) reverted: bumping to 20000 caused align_refine to overshoot true element
 * boundaries by 3-5 kb on synthetic ATHILA-like elements, producing chimeric
 * consensus that MDL rejects. K's banded DP + larger cap lets consensus drift
 * into random background. See V6_PHASE3_RESULT.md BIO-N2 failure. Keep at 10000;
 * for very large LTR elements (>20kb), the per-iteration cap × 10 iters = 100kb
 * still allows full-length recovery. */

/* ================================================================
 * Consensus k-mer hash set (small, for one family at a time)
 * ================================================================ */

typedef struct ckmer_entry {
    uint64_t canon;            /* canonical k-mer */
    uint64_t fwd;              /* forward (non-canonical) packed k-mer */
    int      cons_pos;         /* position in consensus */
    struct ckmer_entry *next;
} CKmerEntry;

typedef struct {
    CKmerEntry **buckets;
    size_t       size;         /* prime table size */
    int          count;
} CKmerSet;

static CKmerSet *ckmer_set_create(int expected)
{
    CKmerSet *cs = calloc(1, sizeof(CKmerSet));
    if (!cs) return NULL;
    /* Use ~2x expected entries for low collision */
    size_t sz = (size_t)(expected * 2 + 1);
    if (sz < 127) sz = 127;
    cs->size = sz;
    cs->buckets = calloc(sz, sizeof(CKmerEntry *));
    if (!cs->buckets) { free(cs); return NULL; }
    return cs;
}

static void ckmer_set_free(CKmerSet *cs)
{
    if (!cs) return;
    for (size_t i = 0; i < cs->size; i++) {
        CKmerEntry *e = cs->buckets[i];
        while (e) {
            CKmerEntry *next = e->next;
            free(e);
            e = next;
        }
    }
    free(cs->buckets);
    free(cs);
}

static void ckmer_set_insert(CKmerSet *cs, uint64_t canon, uint64_t fwd,
                              int cons_pos)
{
    size_t h = (size_t)(canon % cs->size);
    /* Check for duplicate canonical k-mer — keep first occurrence */
    for (CKmerEntry *e = cs->buckets[h]; e; e = e->next) {
        if (e->canon == canon) return;
    }
    CKmerEntry *e = malloc(sizeof(CKmerEntry));
    if (!e) return;
    e->canon = canon;
    e->fwd = fwd;
    e->cons_pos = cons_pos;
    e->next = cs->buckets[h];
    cs->buckets[h] = e;
    cs->count++;
}

static CKmerEntry *ckmer_set_lookup(const CKmerSet *cs, uint64_t canon)
{
    size_t h = (size_t)(canon % cs->size);
    for (CKmerEntry *e = cs->buckets[h]; e; e = e->next) {
        if (e->canon == canon) return e;
    }
    return NULL;
}

/* ================================================================
 * Step 1: Multi-k-mer seeding via position index lookup
 * ================================================================ */

/*
 * For each unique consensus k-mer, look up its positions from the
 * KmerTable position index.  O(sum of freqs) instead of O(genome_length).
 */
static int seed_genome_scan(const CandidateFamily *fam, const KmerTable *kt,
                            int k, SeedHit *hits, int max_hits)
{
    int cons_len = fam->consensus_length;
    if (cons_len < k) return 0;

    int n_cons_kmers = cons_len - k + 1;
    if (n_cons_kmers > MAX_CONS_KMERS) n_cons_kmers = MAX_CONS_KMERS;

    /* Use CKmerSet to deduplicate canonical k-mers in the consensus
     * (same k-mer appearing at multiple positions → only process once) */
    CKmerSet *processed = ckmer_set_create(n_cons_kmers);
    if (!processed) return 0;

    int n_hits = 0;

    for (int i = 0; i < n_cons_kmers && n_hits < max_hits; i++) {
        uint64_t fwd = kmer_pack(fam->consensus + i, k);
        if (fwd == UINT64_MAX) continue;
        uint64_t canon = kmer_canonical(fwd, k);

        /* Skip if we've already processed this canonical k-mer */
        if (ckmer_set_lookup(processed, canon)) continue;
        ckmer_set_insert(processed, canon, fwd, i);

        KmerEntry *entry = kmer_lookup(kt, canon);
        if (!entry || entry->num_positions == 0) continue;

        int cons_is_rc = (fwd != canon);

        for (int pi = 0; pi < entry->num_positions && n_hits < max_hits; pi++) {
            gpos_t sp = entry->positions[pi];
            gpos_t pos = (sp < 0) ? -sp : sp;
            int genome_is_rc = (sp < 0);

            int8_t strand = (cons_is_rc == genome_is_rc) ? +1 : -1;

            hits[n_hits].genome_pos = pos;
            hits[n_hits].cons_pos = i;
            hits[n_hits].strand = strand;
            n_hits++;
        }
    }

    ckmer_set_free(processed);
    return n_hits;
}

/* ================================================================
 * Step 2: Seed clustering (cluster_seed_hits)
 * ================================================================ */

/* Sort seed hits by (strand, genome_pos) */
static int cmp_seed_hits(const void *a, const void *b)
{
    const SeedHit *ha = (const SeedHit *)a;
    const SeedHit *hb = (const SeedHit *)b;
    if (ha->strand != hb->strand) return ha->strand - hb->strand;
    if (ha->genome_pos < hb->genome_pos) return -1;
    if (ha->genome_pos > hb->genome_pos) return 1;
    return 0;
}

/*
 * Cluster nearby same-strand hits into candidate instances.
 * For each cluster, pick the hit closest to consensus center as anchor.
 * Returns number of clusters. Fills anchors[] with selected anchor hits.
 */
static int cluster_seed_hits(SeedHit *hits, int n_hits, int cons_len,
                             SeedHit *anchors, int max_anchors)
{
    if (n_hits == 0) return 0;

    qsort(hits, (size_t)n_hits, sizeof(SeedHit), cmp_seed_hits);

    glen_t merge_dist = (glen_t)(cons_len * 1.5);
    int cons_mid = cons_len / 2;
    int n_anchors = 0;

    int cluster_start = 0;
    for (int i = 1; i <= n_hits && n_anchors < max_anchors; i++) {
        int end_cluster = 0;
        if (i == n_hits) {
            end_cluster = 1;
        } else if (hits[i].strand != hits[i - 1].strand) {
            end_cluster = 1;
        } else if (hits[i].genome_pos - hits[i - 1].genome_pos > merge_dist) {
            end_cluster = 1;
        }

        if (end_cluster) {
            /* Pick anchor closest to consensus center */
            int best = cluster_start;
            int best_dist = abs(hits[cluster_start].cons_pos - cons_mid);
            for (int j = cluster_start + 1; j < i; j++) {
                int d = abs(hits[j].cons_pos - cons_mid);
                if (d < best_dist) {
                    best_dist = d;
                    best = j;
                }
            }
            anchors[n_anchors++] = hits[best];
            cluster_start = i;
        }
    }

    return n_anchors;
}

/* ================================================================
 * Step 3: Banded alignment (align_banded)
 * ================================================================ */

/*
 * Access a genome base, with strand handling.
 * For forward strand: genome[pos]
 * For reverse strand: complement of genome[pos]
 */
static inline char genome_base(const Genome *genome, gpos_t pos, int8_t strand)
{
    if (pos < 0 || pos >= genome->length) return DNA_N;
    char b = genome->sequence[pos];
    if (strand < 0) return dna_complement(b);
    return b;
}

/*
 * Compute the genome position for consensus position i, given an anchor.
 * For forward strand: anchor_genome + (i - anchor_cons)
 * For reverse strand: anchor_genome - (i - anchor_cons)
 *   (reverse strand reads genome backwards)
 */
static inline gpos_t genome_pos_for_cons(gpos_t anchor_genome, int anchor_cons,
                                         int cons_i, int8_t strand)
{
    if (strand > 0)
        return anchor_genome + (cons_i - anchor_cons);
    else
        return anchor_genome - (cons_i - anchor_cons);
}

/*
 * Banded DP alignment extending from a seed anchor.
 * Extends right from anchor, then left, using RepeatScout's scoring.
 *
 * Returns an AlignedInstance with boundaries and score.
 * Returns .score = INT_MIN on failure.
 */
static AlignedInstance align_banded(const CandidateFamily *fam,
                                   const Genome *genome,
                                   const SeedHit *anchor)
{
    AlignedInstance result;
    memset(&result, 0, sizeof(result));
    result.score = INT_MIN;
    result.strand = anchor->strand;

    int cons_len = fam->consensus_length;
    int anchor_cons = anchor->cons_pos;
    gpos_t anchor_genome = anchor->genome_pos;

    /* DP arrays: two rows, each of width (2*maxoffset+1).
     * offset index: o ∈ [0, band_width-1], actual offset = o - maxoffset */
    int band_width = 2 * g_align_maxoffset + 1;
    int dp[2][MAX_BAND_WIDTH];
    int prev_row, curr_row;

    /* --- Extend RIGHT from anchor --- */
    int best_right_score = 0;
    int best_right_i = anchor_cons;
    int right_end_cons = anchor_cons; /* exclusive */

    /* Initialize DP at anchor position */
    prev_row = 0;
    for (int o = 0; o < band_width; o++)
        dp[prev_row][o] = ALIGN_CAPPENALTY;
    dp[prev_row][g_align_maxoffset] = 0; /* offset 0 = exact alignment */

    /* Score the anchor position itself */
    {
        gpos_t gp = genome_pos_for_cons(anchor_genome, anchor_cons,
                                        anchor_cons, anchor->strand);
        char gb = genome_base(genome, gp, anchor->strand);
        char cb = fam->consensus[anchor_cons];
        if (gb != DNA_N && cb != DNA_N) {
            int s = (gb == cb) ? ALIGN_MATCH : ALIGN_MISMATCH;
            dp[prev_row][g_align_maxoffset] = s;
            best_right_score = s;
        }
        right_end_cons = anchor_cons + 1;
    }

    int stall = 0;
    for (int i = anchor_cons + 1; i < cons_len; i++) {
        curr_row = 1 - prev_row;

        for (int o = 0; o < band_width; o++)
            dp[curr_row][o] = ALIGN_CAPPENALTY;

        for (int o = 0; o < band_width; o++) {
            int actual_offset = o - g_align_maxoffset;
            gpos_t gp = genome_pos_for_cons(anchor_genome, anchor_cons,
                                            i, anchor->strand);
            gp += actual_offset;

            char gb = genome_base(genome, gp, anchor->strand);
            char cb = fam->consensus[i];

            if (gb == DNA_N || cb == DNA_N) continue;

            int diag = (gb == cb) ? ALIGN_MATCH : ALIGN_MISMATCH;

            /* Match/mismatch from same offset */
            int val = dp[prev_row][o] + diag;

            /* Gap in genome: offset shifts by +1 (previous o+1) */
            if (o + 1 < band_width && dp[prev_row][o + 1] + g_align_gap > val)
                val = dp[prev_row][o + 1] + g_align_gap;

            /* Gap in consensus: offset shifts by -1 (previous o-1) */
            if (o > 0 && dp[prev_row][o - 1] + g_align_gap > val)
                val = dp[prev_row][o - 1] + g_align_gap;

            /* Cap penalty floor */
            if (val < best_right_score + ALIGN_CAPPENALTY)
                val = best_right_score + ALIGN_CAPPENALTY;

            dp[curr_row][o] = val;

            if (val > best_right_score) {
                best_right_score = val;
                best_right_i = i;
                stall = 0;
            }
        }

        stall++;
        if (stall >= ALIGN_WHEN_TO_STOP) break;

        prev_row = curr_row;
    }
    right_end_cons = best_right_i + 1;

    /* --- Extend LEFT from anchor --- */
    int best_left_score = 0;
    int best_left_i = anchor_cons;
    int left_start_cons = anchor_cons;

    /* Re-initialize DP */
    prev_row = 0;
    for (int o = 0; o < band_width; o++)
        dp[prev_row][o] = ALIGN_CAPPENALTY;
    dp[prev_row][g_align_maxoffset] = 0;

    stall = 0;
    for (int i = anchor_cons - 1; i >= 0; i--) {
        curr_row = 1 - prev_row;

        for (int o = 0; o < band_width; o++)
            dp[curr_row][o] = ALIGN_CAPPENALTY;

        for (int o = 0; o < band_width; o++) {
            int actual_offset = o - g_align_maxoffset;
            gpos_t gp = genome_pos_for_cons(anchor_genome, anchor_cons,
                                            i, anchor->strand);
            gp += actual_offset;

            char gb = genome_base(genome, gp, anchor->strand);
            char cb = fam->consensus[i];

            if (gb == DNA_N || cb == DNA_N) continue;

            int diag = (gb == cb) ? ALIGN_MATCH : ALIGN_MISMATCH;

            int val = dp[prev_row][o] + diag;

            if (o + 1 < band_width && dp[prev_row][o + 1] + g_align_gap > val)
                val = dp[prev_row][o + 1] + g_align_gap;

            if (o > 0 && dp[prev_row][o - 1] + g_align_gap > val)
                val = dp[prev_row][o - 1] + g_align_gap;

            if (val < best_left_score + ALIGN_CAPPENALTY)
                val = best_left_score + ALIGN_CAPPENALTY;

            dp[curr_row][o] = val;

            if (val > best_left_score) {
                best_left_score = val;
                best_left_i = i;
                stall = 0;
            }
        }

        stall++;
        if (stall >= ALIGN_WHEN_TO_STOP) break;

        prev_row = curr_row;
    }
    left_start_cons = best_left_i;

    /* Compute aligned region */
    int aligned_cons_len = right_end_cons - left_start_cons;
    if (aligned_cons_len < 10) return result; /* too short */

    int total_score = best_right_score + best_left_score;
    if (total_score <= 0) return result;

    /* Compute genome coordinates */
    gpos_t gstart = genome_pos_for_cons(anchor_genome, anchor_cons,
                                        left_start_cons, anchor->strand);
    gpos_t gend = genome_pos_for_cons(anchor_genome, anchor_cons,
                                      right_end_cons - 1, anchor->strand);
    if (gstart > gend) {
        gpos_t tmp = gstart;
        gstart = gend;
        gend = tmp;
    }
    gend++; /* exclusive */

    /* Count actual edits in the aligned region.
     *
     * DESIGN NOTE (Audit §C3): This counts substitutions only, not gaps.
     * The banded DP above handles gaps for scoring and boundary detection,
     * but the traceback needed to identify gap positions is not performed.
     * Effect on MDL: instance cost is slightly underestimated → savings
     * biased upward → favors acceptance.  Impact is < 3% (TE copies have
     * gap rates ~1/5 to 1/10 of substitution rates).  A full traceback
     * would add complexity for marginal MDL accuracy improvement. */
    int edits = 0;
    int compared = 0;
    for (int i = left_start_cons; i < right_end_cons; i++) {
        gpos_t gp = genome_pos_for_cons(anchor_genome, anchor_cons,
                                        i, anchor->strand);
        char gb = genome_base(genome, gp, anchor->strand);
        char cb = fam->consensus[i];
        if (gb == DNA_N || cb == DNA_N) continue;
        compared++;
        if (gb != cb) edits++;
    }

    float div = (compared > 0) ? (float)edits / (float)compared : 1.0f;
    if (div > g_align_max_divergence) return result;

    /* Fill result */
    result.genome_start = gstart;
    result.genome_end = gend;
    result.cons_start = left_start_cons;
    result.cons_end = right_end_cons;
    result.num_edits = edits;
    result.divergence = div;
    result.score = total_score;
    result.strand = anchor->strand;

    return result;
}

/* ================================================================
 * Step 4a: RMBlast-based instance recruitment for short families
 * ================================================================ */

/*
 * Resolve rmblastn: honor $RMBLASTN_BIN if set, otherwise find it on PATH.
 * Returns a malloc'd path string on success, NULL if not found.
 * Caller must free() the returned string.
 */
static char *find_rmblastn(void)
{
    /* Explicit override (e.g. set by the conda package or the user) */
    const char *env = getenv("RMBLASTN_BIN");
    if (env && *env && access(env, X_OK) == 0)
        return align_strdup(env);
    /* Resolve rmblastn from PATH via 'which' */
    FILE *fp = popen("which rmblastn 2>/dev/null", "r");
    if (!fp) return NULL;
    char buf[512];
    buf[0] = '\0';
    if (fgets(buf, sizeof(buf), fp)) {
        /* Strip trailing newline */
        size_t l = strlen(buf);
        while (l > 0 && (buf[l-1] == '\n' || buf[l-1] == '\r'))
            buf[--l] = '\0';
    }
    pclose(fp);
    if (buf[0] != '\0' && access(buf, X_OK) == 0)
        return align_strdup(buf);
    return NULL;
}

/* Cached rmblastn path result: 0 = unchecked, 1 = found (g_rmblastn_path set),
 * -1 = not available.  Protected by a simple flag (called from parallel
 * threads, but all threads do the same idempotent check; worst case the
 * check runs multiple times which is harmless). */
static volatile int g_rmblastn_checked = 0;
static char *g_rmblastn_path = NULL;

static int rmblastn_available(void)
{
    if (g_rmblastn_checked == 0) {
        char *p = find_rmblastn();
        if (p) {
            g_rmblastn_path = p;
            g_rmblastn_checked = 1;
        } else {
            g_rmblastn_checked = -1;
        }
    }
    return (g_rmblastn_checked == 1);
}

/*
 * Write genome raw sequences (unpadded, ASCII) to a temp FASTA.
 * seq_raw_starts[i] = raw offset where sequence i starts (for coord remapping).
 * Returns 0 on success, -1 on error.
 */
static int write_genome_fasta(const Genome *genome, const char *path,
                               gpos_t *seq_raw_starts, int max_seqs)
{
    FILE *fp = fopen(path, "w");
    if (!fp) return -1;

    for (int i = 0; i < genome->num_sequences && i < max_seqs; i++) {
        gpos_t raw_start = (i == 0) ? 0 : genome->boundaries[i - 1];
        gpos_t raw_end   = genome->boundaries[i];
        /* Last sequence boundary has sentinel +1; strip it */
        if (i == genome->num_sequences - 1 && raw_end > 0)
            raw_end -= 1;

        seq_raw_starts[i] = raw_start;

        const char *seqid = (genome->sequence_ids && genome->sequence_ids[i])
                            ? genome->sequence_ids[i] : "seq";
        fprintf(fp, ">%s\n", seqid);
        for (gpos_t j = raw_start; j < raw_end; j++) {
            fputc(num_to_char(genome->sequence[PADLENGTH + j]), fp);
            if ((j - raw_start + 1) % 80 == 0) fputc('\n', fp);
        }
        if ((raw_end - raw_start) % 80 != 0) fputc('\n', fp);
    }
    fclose(fp);
    return 0;
}

/*
 * RMBlast recruitment is a bounded supplemental pass for short consensuses.
 * It is not a replacement for the k-mer + banded-DP refinement path:
 * parameters are chosen to recover divergent full/near-full short instances
 * without admitting seed-sized local hits.
 *
 * RepeatMasker/RepeatModeler use rmblastn's repeat-aware task with small
 * word seeds.  We follow that direction, but keep identity/coverage gates
 * explicit because this pipeline scores accepted instances with its own MDL
 * model rather than RepeatMasker's SW score model.  Do not enable
 * RepeatMasker's query-domain masklevel filtering here: it intentionally
 * suppresses overlapping hits on the same consensus query, while this pass
 * needs to recover all genomic instances for downstream duplicate/overlap
 * filtering in mdl-repeat.
 */
#define RMBLAST_SHORT_THRESHOLD  500
#define RMBLAST_MIN_PIDENT       70.0f
#define RMBLAST_MIN_QCOV_HSP     60.0f
#define RMBLAST_WORD_SIZE        7
#define RMBLAST_GAPOPEN          12
#define RMBLAST_GAPEXTEND        2
#define RMBLAST_MAX_HSPS         1000000


/* ================================================================
 * Step 4: align_collect_instances — orchestrator
 * ================================================================ */

int align_collect_instances(CandidateFamily *fam, const Genome *genome,
                            const KmerTable *kt, int k)
{
    /* Free old instances */
    free(fam->instances);
    fam->instances = NULL;
    fam->num_instances = 0;
    fam->cap_instances = 0;

    if (fam->consensus_length < k) return 0;

    /* Step 1: Multi-k-mer genome scan */
    SeedHit *hits = malloc((size_t)MAX_SEED_HITS * sizeof(SeedHit));
    if (!hits) return 0;

    int n_hits = seed_genome_scan(fam, kt, k, hits, MAX_SEED_HITS);
    if (n_hits == 0) {
        free(hits);
        return 0;
    }

    /* Step 2: Cluster hits into candidate instances */
    int max_anchors = DEFAULT_MAXN;
    SeedHit *anchors = malloc((size_t)max_anchors * sizeof(SeedHit));
    if (!anchors) {
        free(hits);
        return 0;
    }

    int n_anchors = cluster_seed_hits(hits, n_hits, fam->consensus_length,
                                      anchors, max_anchors);
    free(hits);

    if (n_anchors == 0) {
        free(anchors);
        return 0;
    }

    /* Pre-allocate instance array */
    int cap = (n_anchors < 64) ? 64 : n_anchors;
    if (cap > DEFAULT_MAXN) cap = DEFAULT_MAXN;
    fam->instances = malloc((size_t)cap * sizeof(Instance));
    if (!fam->instances) {
        free(anchors);
        fam->cap_instances = 0;
        return 0;
    }
    fam->cap_instances = cap;

    /* Minimum instance length: an alignment shorter than this is
     * almost certainly a noisy seed-level match, not a real repeat
     * occurrence.  Without this filter the output BED can be polluted
     * with 11-15 bp fragments (observed on testD, where SINE seeds
     * find spurious tiny matches in the genome background).
     *
     * Threshold: max(k + 10, 30).  Below k+10 the alignment can't
     * extend meaningfully past the seed; below 30 bp it has no
     * biological meaning regardless of seed length. */
    int min_instance_len = (k + 10 > 30) ? (k + 10) : 30;

    /* Step 3: Banded alignment for each anchor */
    for (int i = 0; i < n_anchors && fam->num_instances < DEFAULT_MAXN; i++) {
        AlignedInstance ai = align_banded(fam, genome, &anchors[i]);
        if (ai.score == INT_MIN) continue;

        /* Boundary check */
        glen_t alen = ai.genome_end - ai.genome_start;
        if (ai.genome_start < PADLENGTH) continue;
        if (ai.genome_end > genome->length) continue;
        if (!genome_check_boundary(genome, ai.genome_start, alen))
            continue;
        if (alen < min_instance_len) continue;

        /* Check for duplicate: skip if overlapping an existing instance */
        int is_dup = 0;
        for (int j = 0; j < fam->num_instances; j++) {
            Instance *existing = &fam->instances[j];
            gpos_t ex_end = existing->position + existing->aligned_length;
            /* Overlap if they share > 50% of the shorter one */
            gpos_t ov_start = (ai.genome_start > existing->position)
                              ? ai.genome_start : existing->position;
            gpos_t ov_end = (ai.genome_end < ex_end)
                            ? ai.genome_end : ex_end;
            if (ov_end > ov_start) {
                glen_t overlap = ov_end - ov_start;
                glen_t shorter = (alen < existing->aligned_length)
                                 ? alen : existing->aligned_length;
                if (overlap > shorter / 2) {
                    is_dup = 1;
                    break;
                }
            }
        }
        if (is_dup) continue;

        /* Grow array if needed */
        if (fam->num_instances >= fam->cap_instances) {
            int new_cap = fam->cap_instances * 2;
            if (new_cap > DEFAULT_MAXN) new_cap = DEFAULT_MAXN;
            Instance *tmp = realloc(fam->instances,
                                    (size_t)new_cap * sizeof(Instance));
            if (!tmp) break;
            fam->instances = tmp;
            fam->cap_instances = new_cap;
        }

        Instance *inst = &fam->instances[fam->num_instances];
        inst->position = ai.genome_start;
        inst->aligned_length = ai.genome_end - ai.genome_start;
        inst->cons_start = ai.cons_start;
        inst->cons_end = ai.cons_end;
        inst->num_edits = ai.num_edits;
        inst->divergence = ai.divergence;
        inst->score = ai.score;
        inst->strand = ai.strand;
        inst->seq_index = genome_get_seq_index(genome, ai.genome_start);
        fam->num_instances++;
    }

    free(anchors);
    return fam->num_instances;
}

/* ================================================================
 * Step 5: Consensus rebuild via majority voting
 * ================================================================ */

int align_rebuild_consensus(CandidateFamily *fam, const Genome *genome)
{
    if (fam->num_instances < 2) return 0;

    int cons_len = fam->consensus_length;
    int changed = 0;

    for (int p = 0; p < cons_len; p++) {
        int count[4] = {0, 0, 0, 0};

        for (int j = 0; j < fam->num_instances; j++) {
            Instance *inst = &fam->instances[j];

            /* Only count if this instance covers consensus position p */
            if (p < inst->cons_start || p >= inst->cons_end)
                continue;

            /* Weight by alignment score (higher quality → more influence) */
            int weight = inst->score;
            if (weight < 1) weight = 1;

            /* Map consensus position to genome position.
             * For forward strand: genome_pos = position + (p - cons_start)
             * For reverse strand: genome_pos = position + (cons_end - 1 - p) */
            gpos_t gp;
            if (inst->strand > 0) {
                gp = inst->position + (p - inst->cons_start);
            } else {
                gp = inst->position + (inst->cons_end - 1 - p);
            }

            if (gp < 0 || gp >= genome->length) continue;

            char base = genome->sequence[gp];
            if (inst->strand < 0)
                base = dna_complement(base);

            if (base >= 0 && base <= 3)
                count[(int)base] += weight;
        }

        /* Find majority base */
        int total = count[0] + count[1] + count[2] + count[3];
        if (total == 0) continue;

        int best_base = 0;
        for (int b = 1; b < 4; b++) {
            if (count[b] > count[best_base])
                best_base = b;
        }

        if (fam->consensus[p] != (char)best_base) {
            fam->consensus[p] = (char)best_base;
            changed++;
        }
    }

    return changed;
}

/* ================================================================
 * Step 5b: Consensus extension via flanking context (banded DP)
 * ================================================================ */

/*
 * Compute the genome base contributed by instance `inst` at extension
 * column `y` (0-based, distance from current consensus edge) under
 * a given DP `offset` for direction +1/-1.  Returns DNA_N if out of
 * bounds.  Reverse-strand complementing is applied here.
 *
 * Coordinate mapping mirrors the existing linear formula in
 * align_rebuild_consensus, augmented with the offset adjustment that
 * the banded DP uses to absorb small indels:
 *   direction = +1 (right): consensus position p = consensus_length + y
 *   direction = -1 (left):  consensus position p = -(y + 1)
 *   forward strand: gp = position + (p - cons_start) + offset
 *   reverse strand: gp = position + (cons_end - 1 - p) - offset
 */
static inline char ext_genome_base(const Instance *inst, const Genome *genome,
                                   int direction, int cons_len,
                                   int y, int offset)
{
    int p;
    if (direction > 0) p = cons_len + y;
    else               p = -(y + 1);

    gpos_t gp;
    if (inst->strand > 0) {
        gp = inst->position + (p - inst->cons_start) + offset;
    } else {
        gp = inst->position + (inst->cons_end - 1 - p) - offset;
    }

    if (gp < 0 || gp >= genome->length) return DNA_N;
    char base = genome->sequence[gp];
    if (inst->strand < 0) base = dna_complement(base);
    if (base < 0 || base > 3) return DNA_N;
    return base;
}

/*
 * Mirror of discover.c::compute_score_right adapted for refine-stage
 * extension.  Computes the best DP score for instance n at column y,
 * offset `offset`, assuming consensus base `a` is appended.
 *
 * The scoring follows RepeatScout's three-case structure:
 *   A) gap in sequence  (oldoffset = offset+1)
 *   B) diagonal         (oldoffset = offset)
 *   C) multiple-gap     (oldoffset < offset; takes match if any pos matches)
 *
 * Returns 0 when the underlying genome position is out of bounds for
 * that instance — same convention as discover.c (boundary → score 0).
 */
typedef struct {
    int   N;                    /* number of contributing instances */
    int   maxoffset;
    int   band_width;
    int   direction;            /* +1 right, -1 left */
    int   cons_len;             /* consensus_length at start of extension */
    int **score[2];             /* score[parity][n] -> band_width ints */
    int  *bestbestscore;        /* per-instance best score so far */
    Instance     *insts;        /* contributing instance pointers, length N */
    const Genome *genome;
} ExtDP;

static int ext_compute_score(const ExtDP *D, int y, int n, int offset, char a)
{
    int oldoffset, tempscore, ans;
    int prev = (y - 1) & 1;
    int M = D->maxoffset;

    /* Boundary check: if reading column y at this offset is out of
     * the genome for this instance, abort this DP cell with score 0. */
    char gb = ext_genome_base(&D->insts[n], D->genome, D->direction,
                              D->cons_len, y, offset);
    if (gb == DNA_N) return 0;

    ans = -1000000000;

    /* Case A: gap in sequence (oldoffset = offset+1) */
    if (offset < M) {
        oldoffset = offset + 1;
        tempscore = D->score[prev][n][oldoffset + M] + g_align_gap;
        if (tempscore > ans) ans = tempscore;
    }

    /* Case B: diagonal (oldoffset = offset) */
    oldoffset = offset;
    tempscore = D->score[prev][n][oldoffset + M];
    tempscore += (a == gb) ? ALIGN_MATCH : ALIGN_MISMATCH;
    if (tempscore > ans) ans = tempscore;

    /* Case C: multiple gaps in consensus (oldoffset < offset).
     * Match if any position in the spanned range matches `a`. */
    for (oldoffset = -M; oldoffset < offset; oldoffset++) {
        int ismatch = 0;
        for (int x = oldoffset; x <= offset; x++) {
            char gbx = ext_genome_base(&D->insts[n], D->genome, D->direction,
                                       D->cons_len, y, x);
            if (gbx != DNA_N && a == gbx) { ismatch = 1; break; }
        }
        tempscore = D->score[prev][n][oldoffset + M];
        tempscore += (offset - oldoffset) * g_align_gap;
        tempscore += ismatch ? ALIGN_MATCH : ALIGN_MISMATCH;
        if (tempscore > ans) ans = tempscore;
    }

    return ans;
}

/*
 * Extend consensus in one direction using banded dynamic programming.
 *   direction: +1 = right, -1 = left.
 * Indels of width <= MAXOFFSET in any single instance are absorbed by
 * the band, so a single insertion no longer causes that instance to
 * vote for the wrong column for the rest of the extension.
 *
 * Returns number of bases extended.
 */
static int extend_direction(CandidateFamily *fam, const Genome *genome,
                            int direction, int max_ext)
{
    int extended = 0;

    if (max_ext > ALIGN_MAX_EXTENSION) max_ext = ALIGN_MAX_EXTENSION;
    if (max_ext < 1) return 0;

    int cons_len_at_start = fam->consensus_length;
    int M = g_align_maxoffset;
    if (M > ALIGN_MAXOFFSET_LIMIT) M = ALIGN_MAXOFFSET_LIMIT;
    int band_width = 2 * M + 1;

    /* Collect the instances contributing to this extension (those whose
     * end reaches near the edge being extended).  Same SLACK gating as
     * the column-vote version. */
    int max_insts = fam->num_instances;
    if (max_insts == 0) return 0;
    Instance *contrib = malloc((size_t)max_insts * sizeof(Instance));
    if (!contrib) return 0;
    int N = 0;
    for (int j = 0; j < fam->num_instances; j++) {
        Instance *inst = &fam->instances[j];
        if (direction > 0) {
            if (inst->cons_end < cons_len_at_start - EXTENSION_SLACK) continue;
        } else {
            if (inst->cons_start > EXTENSION_SLACK) continue;
        }
        contrib[N++] = *inst;
    }
    if (N < 2) { free(contrib); return 0; }

    /* Allocate score arrays: score[2][N][band_width] as flat block. */
    int *score_block = malloc((size_t)2 * N * band_width * sizeof(int));
    int **row0 = malloc((size_t)N * sizeof(int *));
    int **row1 = malloc((size_t)N * sizeof(int *));
    int *bestbest = malloc((size_t)N * sizeof(int));
    char *ext_buf = malloc((size_t)max_ext);
    if (!score_block || !row0 || !row1 || !bestbest || !ext_buf) {
        free(score_block); free(row0); free(row1);
        free(bestbest); free(ext_buf); free(contrib);
        return 0;
    }
    for (int n = 0; n < N; n++) {
        row0[n] = score_block + (0 * N + n) * band_width;
        row1[n] = score_block + (1 * N + n) * band_width;
    }

    ExtDP D;
    D.N = N;
    D.maxoffset = M;
    D.band_width = band_width;
    D.direction = direction;
    D.cons_len = cons_len_at_start;
    D.score[0] = row0;
    D.score[1] = row1;
    D.bestbestscore = bestbest;
    D.insts = contrib;
    D.genome = genome;

    /* Initialize column "y = -1" (the imaginary anchor inside the
     * existing consensus).  We seed all offsets to 0 so the DP has a
     * neutral starting point — analogous to discover.c initialising
     * from `l * MATCH` after l-mer seed.  Here we use 0 because no
     * pre-extension match has been accumulated; the DP starts fresh. */
    {
        int prev = (-1) & 1;  /* equivalent to ((-1) % 2 + 2) % 2 = 1 */
        for (int n = 0; n < N; n++) {
            for (int o = 0; o < band_width; o++) {
                int actual = o - M;
                int gp_pen = 0;
                if (actual < 0) gp_pen = -actual * g_align_gap;
                if (actual > 0) gp_pen =  actual * g_align_gap;
                D.score[prev][n][o] = gp_pen;  /* offset 0 = 0; off-band = neg */
            }
            bestbest[n] = 0;
        }
    }

    int besty = -1;                /* best column index so far */
    int besttotalbestscore = 0;
    char besta = 0;

    /* Main extension loop */
    int y;
    for (y = 0; y < max_ext; y++) {
        int newtotalbestscore = 0;
        int besta_local = 0;
        int curr = y & 1;

        /* For each candidate base, compute total best score */
        for (char a = 0; a < 4; a++) {
            int total_a = 0;
            for (int n = 0; n < N; n++) {
                int floor_score = bestbest[n] + ALIGN_CAPPENALTY;
                if (floor_score < 0) floor_score = 0;
                int bestscore_n = floor_score;
                for (int off = -M; off <= M; off++) {
                    int s = ext_compute_score(&D, y, n, off, a);
                    if (s > bestscore_n) bestscore_n = s;
                }
                total_a += bestscore_n;
            }
            if (total_a > newtotalbestscore) {
                newtotalbestscore = total_a;
                besta_local = a;
            }
        }
        besta = (char)besta_local;

        /* Commit the chosen base: fill score[curr][n][*] for besta */
        for (int n = 0; n < N; n++) {
            for (int off = -M; off <= M; off++) {
                D.score[curr][n][off + M] =
                    ext_compute_score(&D, y, n, off, besta);
            }
        }

        /* Update bestbestscore and totalbestscore (analogue of
         * compute_totalbestscore_right minus the savebestscore /
         * score_of_besty bookkeeping that's only needed for a
         * subsequent left-extension hand-off). */
        int totalbestscore = 0;
        int nrepeatocc = 0;
        for (int n = 0; n < N; n++) {
            int floor_score = bestbest[n] + ALIGN_CAPPENALTY;
            if (floor_score < 0) floor_score = 0;
            int bs = floor_score;
            for (int off = -M; off <= M; off++) {
                if (D.score[curr][n][off + M] > bs)
                    bs = D.score[curr][n][off + M];
            }
            if (bs > 0) nrepeatocc++;
            if (bs > bestbest[n]) bestbest[n] = bs;
            totalbestscore += bs;
        }

        ext_buf[y] = besta;

        /* Segdup-prevention checks (carried from column-vote version):
         * require minimum number of contributing instances as the
         * extension grows long, so 2-3 identical copies can't drag the
         * consensus along forever. */
        if (nrepeatocc < 2) break;
        if (y > 1000 && nrepeatocc < 3) break;
        if (y > 5000 && nrepeatocc < 5) break;

        /* Column-majority hard stop (mirror of HEAD column-vote).
         * Without this, K's banded DP keeps extending on noise when N is
         * large (e.g. N=27 post-merge): random background alignments in
         * 5-6 of 27 instances are enough to keep nrepeatocc-pos and
         * pick some besta, but no real majority exists. Result: 5 kb
         * overshoot past true element boundary on synthetic ATHILA
         * (V6_PHASE3_RESULT.md BIO-N2 failure).
         *
         * Count instances whose besta vote is strongly supported
         * (per-instance best score with besta > floor_score by ≥ MATCH).
         * If less than 50% of N support, this column is noise; break. */
        int n_majority = 0;
        for (int n = 0; n < N; n++) {
            int floor_score = bestbest[n] + ALIGN_CAPPENALTY;
            if (floor_score < 0) floor_score = 0;
            for (int off = -M; off <= M; off++) {
                if (D.score[curr][n][off + M] >= floor_score + ALIGN_MATCH) {
                    n_majority++;
                    break;
                }
            }
        }
        if (n_majority * 2 < N) break;

        /* Score plateau / WHEN_TO_STOP termination */
        if (totalbestscore > besttotalbestscore) {
            besttotalbestscore = totalbestscore;
            besty = y;
        }
        /* Adaptive WHEN_TO_STOP — long extensions earn a longer quiet
         * window (mirrors discover.c::extend_right). */
        int extended_so_far = besty + 1;
        int adaptive_stop = ALIGN_WHEN_TO_STOP;
        if (extended_so_far / 10 > adaptive_stop)
            adaptive_stop = extended_so_far / 10;
        if (y - besty >= adaptive_stop) break;
    }

    /* Trim to best position: extended = besty + 1 (besty is 0-based) */
    extended = (besty >= 0) ? (besty + 1) : 0;

    free(score_block);
    free(row0);
    free(row1);
    free(bestbest);
    free(contrib);

    if (extended == 0) {
        free(ext_buf);
        return 0;
    }

    /* Reallocate consensus */
    int new_len = fam->consensus_length + extended;
    char *new_cons = realloc(fam->consensus, (size_t)new_len);
    if (!new_cons) {
        free(ext_buf);
        return 0;
    }

    if (direction > 0) {
        /* Append to right */
        memcpy(new_cons + fam->consensus_length, ext_buf, (size_t)extended);
    } else {
        /* Prepend to left: shift existing consensus right, insert reversed */
        memmove(new_cons + extended, new_cons, (size_t)fam->consensus_length);
        for (int i = 0; i < extended; i++)
            new_cons[i] = ext_buf[extended - 1 - i];

        /* Shift all instance coordinates to match new consensus origin */
        for (int j = 0; j < fam->num_instances; j++) {
            fam->instances[j].cons_start += extended;
            fam->instances[j].cons_end += extended;
        }
    }

    fam->consensus = new_cons;
    fam->consensus_length = new_len;

    free(ext_buf);
    return extended;
}

int align_extend_consensus(CandidateFamily *fam, const Genome *genome)
{
    if (fam->num_instances < 2) return 0;

    /* Extension cap based on instance count.
     * More instances = higher confidence = allow longer extension per iteration.
     * Few instances (2-3) get conservative cap to prevent segdup runaway.
     * Many instances (10+) get full extension — the WHEN_TO_STOP and CAPPENALTY
     * mechanisms naturally stop at the TE boundary.
     * This replaces the old consensus-length-based cap that limited growth of
     * initially short graph unitigs (median 222bp → stuck at short lengths). */
    int cap;
    if (fam->num_instances <= 3)
        cap = 500;
    else if (fam->num_instances <= 9)
        cap = 3000;
    else
        cap = ALIGN_MAX_EXTENSION;

    int ext_right = extend_direction(fam, genome, +1, cap);
    int ext_left  = extend_direction(fam, genome, -1, cap);

    return ext_right + ext_left;
}

/* ================================================================
 * Step 6: Iteration loop
 * ================================================================ */

int align_refine_family(CandidateFamily *fam, const Genome *genome,
                        const KmerTable *kt, int k, int max_iter)
{
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        align_collect_instances(fam, genome, kt, k);

        if (fam->num_instances < 2)
            break;

        int extended = align_extend_consensus(fam, genome);
        int changed = align_rebuild_consensus(fam, genome);

        /* Converged: no extension and consensus stable */
        if (extended == 0 && changed == 0)
            break;
        if (extended == 0 && changed < fam->consensus_length / 100 + 1)
            break;
    }

    return fam->num_instances;
}

/* ================================================================
 * Step 7: Batch refinement (parallel)
 * ================================================================ */

/* Shared state for worker threads */
typedef struct {
    CandidateList  *cl;
    const Genome   *genome;
    const KmerTable *kt;
    int             k;
    int             verbose;
    int             next_family;    /* atomic: next family index to process */
    int             done_count;     /* atomic: families completed so far */
    pthread_mutex_t progress_lock;  /* protects stderr progress output */
} RefineWorkerState;

static void *refine_worker(void *arg)
{
    RefineWorkerState *state = (RefineWorkerState *)arg;

    while (1) {
        int idx = __atomic_fetch_add(&state->next_family, 1, __ATOMIC_SEQ_CST);
        if (idx >= state->cl->num_families)
            break;

        CandidateFamily *fam = &state->cl->families[idx];
        align_refine_family(fam, state->genome, state->kt, state->k,
                            ALIGN_MAX_ITERATIONS);

        int done = __atomic_add_fetch(&state->done_count, 1, __ATOMIC_SEQ_CST);
        if (state->verbose && done % 100 == 0) {
            pthread_mutex_lock(&state->progress_lock);
            fprintf(stderr, "  Refined %d / %d families\r",
                    done, state->cl->num_families);
            pthread_mutex_unlock(&state->progress_lock);
        }
    }

    return NULL;
}

void align_refine_all(CandidateList *cl, const Genome *genome,
                      const KmerTable *kt, int k, int num_threads,
                      int verbose)
{
    if (!cl || cl->num_families == 0) return;

    if (num_threads < 1) num_threads = 1;

    /* Cap threads at family count */
    if (num_threads > cl->num_families)
        num_threads = cl->num_families;

    if (num_threads == 1) {
        /* Sequential path (no thread overhead) */
        int refined = 0, lost = 0;
        for (int i = 0; i < cl->num_families; i++) {
            CandidateFamily *fam = &cl->families[i];
            align_refine_family(fam, genome, kt, k, ALIGN_MAX_ITERATIONS);

            if (fam->num_instances >= 2) refined++;
            else lost++;

            if (verbose && (i + 1) % 100 == 0)
                fprintf(stderr, "  Refined %d / %d families\r",
                        i + 1, cl->num_families);
        }
        if (verbose)
            fprintf(stderr, "  Alignment refinement: %d families refined, "
                    "%d lost all instances\n", refined, lost);
        return;
    }

    /* Parallel path */
    if (verbose)
        fprintf(stderr, "  Using %d threads\n", num_threads);

    RefineWorkerState state;
    state.cl = cl;
    state.genome = genome;
    state.kt = kt;
    state.k = k;
    state.verbose = verbose;
    state.next_family = 0;
    state.done_count = 0;
    pthread_mutex_init(&state.progress_lock, NULL);

    pthread_t *threads = malloc((size_t)num_threads * sizeof(pthread_t));
    if (!threads) {
        fprintf(stderr, "WARNING: thread allocation failed, falling back to single-threaded\n");
        align_refine_all(cl, genome, kt, k, 1, verbose);
        return;
    }

    /* Launch workers */
    for (int t = 0; t < num_threads; t++) {
        if (pthread_create(&threads[t], NULL, refine_worker, &state) != 0) {
            fprintf(stderr, "WARNING: pthread_create failed for thread %d, "
                    "continuing with %d threads\n", t, t);
            num_threads = t;
            break;
        }
    }

    /* Wait for all workers */
    for (int t = 0; t < num_threads; t++)
        pthread_join(threads[t], NULL);

    free(threads);
    pthread_mutex_destroy(&state.progress_lock);

    /* Count results */
    int refined = 0, lost = 0;
    for (int i = 0; i < cl->num_families; i++) {
        if (cl->families[i].num_instances >= 2) refined++;
        else lost++;
    }

    if (verbose)
        fprintf(stderr, "  Alignment refinement: %d families refined, "
                "%d lost all instances\n", refined, lost);
}

/* ================================================================
 * Batch RMBlast recruitment for short families
 * ================================================================ */

/*
 * Add an Instance to a CandidateFamily, growing the array as needed.
 * Returns 1 on success, 0 on allocation failure.
 */
static int family_add_instance(CandidateFamily *fam, const Instance *inst)
{
    if (fam->num_instances >= DEFAULT_MAXN) return 0;

    if (!fam->instances) {
        int cap = 64;
        fam->instances = malloc((size_t)cap * sizeof(Instance));
        if (!fam->instances) return 0;
        fam->cap_instances = cap;
        fam->num_instances = 0;
    }

    if (fam->num_instances >= fam->cap_instances) {
        int new_cap = fam->cap_instances * 2;
        if (new_cap > DEFAULT_MAXN) new_cap = DEFAULT_MAXN;
        Instance *tmp = realloc(fam->instances,
                                (size_t)new_cap * sizeof(Instance));
        if (!tmp) return 0;
        fam->instances = tmp;
        fam->cap_instances = new_cap;
    }

    fam->instances[fam->num_instances++] = *inst;
    return 1;
}

/*
 * Check if a new instance overlaps an existing one by > 50% of the shorter.
 */
static int instance_overlaps_existing(const CandidateFamily *fam,
                                       gpos_t genome_start, gpos_t genome_end)
{
    glen_t alen = genome_end - genome_start;
    for (int j = 0; j < fam->num_instances; j++) {
        const Instance *ex = &fam->instances[j];
        gpos_t ex_end = ex->position + ex->aligned_length;
        gpos_t ov_s = (genome_start > ex->position) ? genome_start : ex->position;
        gpos_t ov_e = (genome_end   < ex_end)        ? genome_end   : ex_end;
        if (ov_e > ov_s) {
            glen_t ovl = ov_e - ov_s;
            glen_t shorter = (alen < ex->aligned_length) ? alen : ex->aligned_length;
            if (ovl > shorter / 2) return 1;
        }
    }
    return 0;
}

int align_blast_recruit_short_families(CandidateList *cl,
                                        const Genome *genome, int k,
                                        int verbose)
{
    if (!rmblastn_available()) return -1;

    /* Identify families with consensus_length < RMBLAST_SHORT_THRESHOLD */
    int n_short = 0;
    for (int i = 0; i < cl->num_families; i++) {
        CandidateFamily *f = &cl->families[i];
        if (f->consensus && f->consensus_length < RMBLAST_SHORT_THRESHOLD
                         && f->consensus_length >= k)
            n_short++;
    }
    if (n_short == 0) return 0;

    if (verbose)
        fprintf(stderr, "  RMBlast recruitment: %d short families (<%d bp) — "
                "running batch rmblastn with %s...\n",
                n_short, RMBLAST_SHORT_THRESHOLD, g_rmblastn_path);

    /* Build temp file paths */
    int pid = (int)getpid();
    char query_fa[256], subj_fa[256], out_tab[256];
    snprintf(query_fa, sizeof(query_fa), "/tmp/mdlr_batch_q_%d.fa",  pid);
    snprintf(subj_fa,  sizeof(subj_fa),  "/tmp/mdlr_batch_s_%d.fa",  pid);
    snprintf(out_tab,  sizeof(out_tab),  "/tmp/mdlr_batch_o_%d.tab", pid);

    /* Map from genome sequence names to raw start offsets */
    int nseq = genome->num_sequences;
    if (nseq < 1) return -1;
    gpos_t *seq_raw_starts = malloc((size_t)nseq * sizeof(gpos_t));
    if (!seq_raw_starts) return -1;

    int rc = -1;

    /* Write genome subject FASTA */
    if (write_genome_fasta(genome, subj_fa, seq_raw_starts, nseq) != 0)
        goto cleanup;

    /* Write multi-FASTA query: one record per short family.
     * Header format: >F<family_index>  (e.g., ">F42") */
    {
        FILE *fp = fopen(query_fa, "w");
        if (!fp) goto cleanup;
        for (int i = 0; i < cl->num_families; i++) {
            CandidateFamily *f = &cl->families[i];
            if (!f->consensus || f->consensus_length < k
                               || f->consensus_length >= RMBLAST_SHORT_THRESHOLD)
                continue;
            fprintf(fp, ">F%d\n", i);
            for (int j = 0; j < f->consensus_length; j++) {
                fputc(num_to_char(f->consensus[j]), fp);
                if ((j + 1) % 80 == 0) fputc('\n', fp);
            }
            if (f->consensus_length % 80 != 0) fputc('\n', fp);
        }
        fclose(fp);
    }

    /* Run single rmblastn for all short families */
    {
        char cmd[2048];
        snprintf(cmd, sizeof(cmd),
                 "%s -task rmblastn -query %s -subject %s "
                 "-outfmt \"6 qseqid sseqid sstart send pident length "
                 "qstart qend evalue sstrand\" "
                 "-word_size %d -perc_identity %.1f -qcov_hsp_perc %.1f "
                 "-gapopen %d -gapextend %d "
                 "-max_hsps %d -dust no "
                 "-out %s 2>/dev/null",
                 g_rmblastn_path, query_fa, subj_fa,
                 RMBLAST_WORD_SIZE, (double)RMBLAST_MIN_PIDENT,
                 (double)RMBLAST_MIN_QCOV_HSP,
                 RMBLAST_GAPOPEN, RMBLAST_GAPEXTEND, RMBLAST_MAX_HSPS,
                 out_tab);
        system(cmd);  /* non-zero exit = no hits, not a hard error */
    }

    /* Parse BLAST output and distribute hits to families */
    FILE *fp = fopen(out_tab, "r");
    if (!fp) { rc = 0; goto cleanup; }

    char line[1024];
    int total_added = 0;

    while (fgets(line, sizeof(line), fp)) {
        /* Fields: qseqid sseqid sstart send pident length qstart qend evalue sstrand */
        char qseqid[64], sseqid[256], sstrand_str[16];
        long long sstart, send, qstart, qend, length;
        double pident, evalue;

        if (sscanf(line, "%63s %255s %lld %lld %lf %lld %lld %lld %lf %15s",
                   qseqid, sseqid, &sstart, &send, &pident, &length,
                   &qstart, &qend, &evalue, sstrand_str) != 10)
            continue;

        /* Parse family index from qseqid "F<idx>" */
        if (qseqid[0] != 'F') continue;
        int fam_idx = atoi(qseqid + 1);
        if (fam_idx < 0 || fam_idx >= cl->num_families) continue;

        CandidateFamily *fam = &cl->families[fam_idx];
        if (!fam->consensus) continue;

        /* Filters */
        if (pident < (double)RMBLAST_MIN_PIDENT) continue;
        if (length < (long long)k) continue;

        /* Strand */
        int8_t strand = (strcmp(sstrand_str, "minus") == 0) ? -1 : +1;

        /* Find genome sequence */
        int seq_idx = -1;
        for (int i = 0; i < nseq; i++) {
            const char *seqid = (genome->sequence_ids && genome->sequence_ids[i])
                                ? genome->sequence_ids[i] : "seq";
            if (strcmp(seqid, sseqid) == 0) { seq_idx = i; break; }
        }
        if (seq_idx < 0) continue;

        /* Convert BLAST coords → padded internal coords */
        gpos_t raw_start, raw_end;
        if (strand > 0) {
            raw_start = seq_raw_starts[seq_idx] + (sstart - 1);
            raw_end   = seq_raw_starts[seq_idx] + send;
        } else {
            raw_start = seq_raw_starts[seq_idx] + (send - 1);
            raw_end   = seq_raw_starts[seq_idx] + sstart;
        }
        gpos_t genome_start = raw_start + PADLENGTH;
        gpos_t genome_end   = raw_end   + PADLENGTH;
        glen_t alen         = genome_end - genome_start;

        /* Validate */
        if (genome_start < PADLENGTH) continue;
        if (genome_end > genome->length) continue;
        if (alen < (glen_t)k) continue;
        if (!genome_check_boundary(genome, genome_start, alen)) continue;

        /* Consensus coords: qstart/qend are 1-based in query (always forward) */
        int cons_start = (int)(qstart - 1);
        int cons_end   = (int)qend;
        if (cons_start < 0) cons_start = 0;
        if (cons_end > fam->consensus_length) cons_end = fam->consensus_length;
        if (cons_end <= cons_start) continue;

        float div = (float)((100.0 - pident) / 100.0);
        if (div > g_align_max_divergence) continue;

        /* Skip if overlapping existing instance */
        if (instance_overlaps_existing(fam, genome_start, genome_end)) continue;

        Instance inst;
        inst.position       = genome_start;
        inst.aligned_length = alen;
        inst.cons_start     = cons_start;
        inst.cons_end       = cons_end;
        inst.num_edits      = (int)(div * (float)length + 0.5f);
        inst.divergence     = div;
        inst.score          = (int)(pident / 100.0 * (double)length);
        if (inst.score < 1) inst.score = 1;
        inst.strand         = strand;
        inst.seq_index      = genome_get_seq_index(genome, genome_start);

        if (family_add_instance(fam, &inst))
            total_added++;
    }

    fclose(fp);
    rc = total_added;

    if (verbose)
        fprintf(stderr, "  RMBlast recruitment: added %d instances across %d families\n",
                total_added, n_short);

cleanup:
    free(seq_raw_starts);
    remove(query_fa);
    remove(subj_fa);
    remove(out_tab);
    return rc;
}
