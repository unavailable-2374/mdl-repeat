/*
 * main.c — MDL-guided seed-and-extend de novo repeat discovery
 *
 * Pipeline:
 *   1. Load genome
 *   2. Discover consensus families (RepeatScout-style seed-and-extend)
 *   3. Build k-mer table + position index (for refine pipeline)
 *   3b. Align-refine all families (find all genome instances)
 *   4. Compact (remove dead/short families)
 *   5. Merge redundant families
 *   6. Split bimodal families
 *   7. MDL scoring and library selection
 *   8. Prune marginal families
 *   9. Output (FASTA, BED, TSV)
 *
 * -trace-dir <dir>: dump FASTA+BED after each stage for recall diagnostics.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <errno.h>
#include <sys/stat.h>

#include "types.h"
#include "genome.h"
#include "kmer.h"
#include "candidates.h"
#include "discover.h"
#include "align.h"
#include "mdl.h"
#include "refine.h"
#include "output.h"
#include "cmd_line_opts.h"
#include "rescue.h"
#include "external_qc.h"

/* ================================================================
 * Trace-dump helpers (enabled by -trace-dir <dir>)
 * Writes ALL families (regardless of mdl_score) for diagnostic purposes.
 * ================================================================ */


/*
 * Dump one pipeline stage: writes .fa, .bed, .tsv under trace_dir.
 * stage_num is 1-based (01..08), stage_name is a short label.
 * genome may be NULL at stages where instance coords are not meaningful.
 * mdl_filter: if non-zero, only dump accepted families
 *             (use after mdl_select_library to show only accepted families).
 */
static void trace_dump_filtered(const char *trace_dir, int stage_num,
                                const char *stage_name,
                                const CandidateList *cl, const Genome *genome,
                                int mdl_filter)
{
    char path[4096];
    int n_fam = 0, n_inst = 0;

    /* Build a temporary view — families that pass the filter */
    CandidateList view;
    view.families = cl->families;
    view.num_families = cl->num_families;
    view.cap_families = cl->cap_families;

    if (mdl_filter) {
        /* Count only accepted families for the header */
        for (int i = 0; i < cl->num_families; i++)
            if (candidate_is_accepted(&cl->families[i])) n_fam++;
    }

    snprintf(path, sizeof(path), "%s/%02d_%s.fa", trace_dir, stage_num, stage_name);
    {
        FILE *fp = fopen(path, "w");
        if (fp) {
            int w = 0;
            for (int i = 0; i < cl->num_families; i++) {
                const CandidateFamily *f = &cl->families[i];
                if (!f->consensus || f->consensus_length <= 0) continue;
                if (mdl_filter && !candidate_is_accepted(f)) continue;
                fprintf(fp, ">R=%d length=%d copies=%d mdl=%.1f\n",
                        f->id, f->consensus_length, f->num_instances,
                        candidate_report_score(f));
                for (int j = 0; j < f->consensus_length; j++) {
                    fputc(num_to_char(f->consensus[j]), fp);
                    if ((j + 1) % 80 == 0) fputc('\n', fp);
                }
                if (f->consensus_length % 80 != 0) fputc('\n', fp);
                w++;
            }
            fclose(fp);
            if (!mdl_filter) n_fam = w;
        }
    }

    snprintf(path, sizeof(path), "%s/%02d_%s.bed", trace_dir, stage_num, stage_name);
    {
        FILE *fp = fopen(path, "w");
        if (fp) {
            for (int i = 0; i < cl->num_families; i++) {
                const CandidateFamily *f = &cl->families[i];
                if (!f->instances || f->num_instances <= 0) continue;
                if (mdl_filter && !candidate_is_accepted(f)) continue;
                for (int j = 0; j < f->num_instances; j++) {
                    const Instance *inst = &f->instances[j];
                    gpos_t raw_start = inst->position - PADLENGTH;
                    gpos_t raw_end = raw_start + (glen_t)inst->aligned_length;
                    const char *chrom = "unknown";
                    gpos_t chr_offset = 0;
                    if (genome && inst->seq_index >= 0
                               && inst->seq_index < genome->num_sequences
                               && genome->sequence_ids
                               && genome->sequence_ids[inst->seq_index]) {
                        chrom = genome->sequence_ids[inst->seq_index];
                        if (inst->seq_index > 0)
                            chr_offset = genome->boundaries[inst->seq_index - 1];
                    }
                    gpos_t local_start = raw_start - chr_offset;
                    gpos_t local_end   = raw_end   - chr_offset;
                    int score = (int)(1000 * (1.0 - inst->divergence));
                    if (score < 0) score = 0;
                    if (score > 1000) score = 1000;
                    fprintf(fp, "%s\t%" PRId64 "\t%" PRId64 "\tR=%d\t%d\t%c\n",
                            chrom, (int64_t)local_start, (int64_t)local_end,
                            f->id, score, inst->strand > 0 ? '+' : '-');
                    n_inst++;
                }
            }
            fclose(fp);
        }
    }

    snprintf(path, sizeof(path), "%s/%02d_%s.tsv", trace_dir, stage_num, stage_name);
    {
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "family_id\tconsensus_length\tnum_instances\t"
                        "divergence_mean\tmdl_score\tmodel_cost\ttopology\t"
                        "standalone_score\texclusive_score\texclusive_bases\t"
                        "exclusive_instances\tacceptance\tquality_tier\tquality_flags\t"
                        "quality_notes\n");
            for (int i = 0; i < cl->num_families; i++) {
                const CandidateFamily *f = &cl->families[i];
                if (mdl_filter && !candidate_is_accepted(f)) continue;
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
                        avg_div, candidate_report_score(f),
                        candidate_model_cost(f), topo,
                        f->mdl.standalone_score, f->mdl.exclusive_score,
                        f->mdl.exclusive_bases, f->mdl.exclusive_instances,
                        candidate_accept_state_name(f->mdl.accept_state),
                        candidate_quality_tier_name(f->mdl.quality_tier),
                        (unsigned)f->mdl.quality_flags, qflags);
            }
            fclose(fp);
        }
    }

    fprintf(stderr, "[trace] stage %02d %-12s: %d families, %d instances%s\n",
            stage_num, stage_name, n_fam, n_inst,
            mdl_filter ? " (accepted only)" : "");
    (void)view; /* suppress unused warning */
}

/*
 * Convenience wrapper — no mdl_score filter (stages 01-05).
 */
static void trace_dump(const char *trace_dir, int stage_num,
                       const char *stage_name,
                       const CandidateList *cl, const Genome *genome)
{
    trace_dump_filtered(trace_dir, stage_num, stage_name, cl, genome, 0);
}

#define DISCOVER_SPLIT_THRESHOLD_DEFAULT  200000000LL  /* 200M bases */

static int default_k(glen_t len)
{
    return (int)ceil(1.0 + log((double)len) / log(4.0));
}

/* ----------------------------------------------------------------
 * Chunked discovery helpers (large genomes)
 * ---------------------------------------------------------------- */

/*
 * Segment descriptor: sub-range of a genome sequence.
 * Long sequences (> 1.8 * chunk_size) are recursively halved with overlap
 * at split points to avoid missing repeats at boundaries.
 */
static int cmp_seg_length_desc(const void *a, const void *b)
{
    glen_t la = ((const SeqSegment *)a)->seg_length;
    glen_t lb = ((const SeqSegment *)b)->seg_length;
    if (lb > la) return 1;
    if (lb < la) return -1;
    return 0;
}

/*
 * Create a sub-Genome from selected segments of the original genome.
 * Each segment is a sub-range [raw_start, raw_end) of a sequence.
 * The new Genome has its own sequence buffer and boundaries array,
 * but shares sequence_ids string pointers (do NOT free those individually).
 */
static Genome *genome_create_chunk(const Genome *g, const SeqSegment *segments, int num_segs)
{
    glen_t total_len = 0;
    for (int i = 0; i < num_segs; i++)
        total_len += segments[i].seg_length;

    Genome *chunk = calloc(1, sizeof(Genome));
    if (!chunk) return NULL;

    chunk->raw_length = total_len;
    chunk->length = PADLENGTH + total_len;
    chunk->num_sequences = num_segs;

    chunk->sequence = malloc((size_t)(chunk->length + 1));
    if (!chunk->sequence) { free(chunk); return NULL; }

    for (glen_t p = 0; p < PADLENGTH; p++)
        chunk->sequence[p] = DNA_N;

    chunk->boundaries = malloc((size_t)(num_segs + 1) * sizeof(gpos_t));
    if (!chunk->boundaries) { free(chunk->sequence); free(chunk); return NULL; }

    chunk->sequence_ids = malloc((size_t)(num_segs + 1) * sizeof(char *));
    if (!chunk->sequence_ids) {
        free(chunk->boundaries);
        free(chunk->sequence);
        free(chunk);
        return NULL;
    }

    glen_t offset = 0;
    for (int i = 0; i < num_segs; i++) {
        memcpy(chunk->sequence + PADLENGTH + offset,
               g->sequence + PADLENGTH + segments[i].raw_start,
               (size_t)segments[i].seg_length);

        offset += segments[i].seg_length;
        chunk->boundaries[i] = offset;
        chunk->sequence_ids[i] = g->sequence_ids[segments[i].seq_index];
    }

    /* +1 sentinel on last boundary, matching genome_load convention */
    chunk->boundaries[num_segs - 1] = offset + 1;
    chunk->boundaries[num_segs] = 0;
    chunk->sequence_ids[num_segs] = NULL;

    return chunk;
}

/*
 * Free a chunk Genome. Does NOT free individual sequence_ids strings
 * (they are shared with the original Genome).
 */
static void genome_free_chunk(Genome *chunk)
{
    if (!chunk) return;
    free(chunk->sequence);
    free(chunk->boundaries);
    free(chunk->sequence_ids);  /* free array only, not the strings */
    free(chunk);
}

/*
 * Write a Genome struct to FASTA file.
 * Each segment becomes one FASTA record with its original sequence name
 * plus a suffix indicating the window coordinates.
 */
static int genome_write_fasta(const Genome *g, const char *filename,
                               const SeqSegment *segments, int num_segments)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "genome_write_fasta: could not open '%s'\n", filename);
        return -1;
    }

    gpos_t offset = 0;
    for (int i = 0; i < num_segments; i++) {
        const char *name = (g->sequence_ids && g->sequence_ids[i])
                           ? g->sequence_ids[i] : "unnamed";
        fprintf(fp, ">%s:%"PRId64"-%"PRId64"\n",
                name, (int64_t)segments[i].raw_start, (int64_t)segments[i].raw_end);

        for (glen_t j = 0; j < segments[i].seg_length; j++) {
            fputc(num_to_char(g->sequence[PADLENGTH + offset + j]), fp);
            if ((j + 1) % 80 == 0) fputc('\n', fp);
        }
        if (segments[i].seg_length % 80 != 0) fputc('\n', fp);
        offset += segments[i].seg_length;
    }

    fclose(fp);
    fprintf(stderr, "Wrote sampled genome (%d windows) to %s\n", num_segments, filename);
    return 0;
}

static int cmp_int_asc(const void *a, const void *b)
{
    return *(const int *)a - *(const int *)b;
}

/*
 * Sample fixed-size windows from a large genome.
 * Divides each sequence into non-overlapping tiles of window_size,
 * then randomly selects tiles to build a sampled Genome.
 *
 * Returns the sampled Genome and the segment mapping for coordinate remapping.
 * Caller must free segments with free() and sampled Genome with genome_free_chunk().
 */
static Genome *genome_sample_windows(const Genome *genome,
                                      glen_t sample_size,
                                      glen_t window_size,
                                      SeqSegment **out_segments,
                                      int *out_num_segments,
                                      unsigned int seed,
                                      int verbose)
{
    int n = genome->num_sequences;

    /* Count tiles per sequence and cumulative offsets */
    int *tiles_per_seq = calloc((size_t)n, sizeof(int));
    int *tile_offset = calloc((size_t)n, sizeof(int));
    int total_tiles = 0;

    for (int i = 0; i < n; i++) {
        gpos_t start = (i == 0) ? 0 : genome->boundaries[i - 1];
        gpos_t end = genome->boundaries[i];
        glen_t seq_len = end - start;
        if (i == genome->num_sequences - 1) seq_len -= 1;

        tile_offset[i] = total_tiles;
        tiles_per_seq[i] = (int)(seq_len / window_size);
        total_tiles += tiles_per_seq[i];
    }

    int num_windows = (int)(sample_size / window_size);
    if (num_windows > total_tiles) num_windows = total_tiles;
    if (num_windows < 1) num_windows = 1;

    if (verbose)
        fprintf(stderr, "  Sampling: %d tiles of %.0fkb available, selecting %d\n",
                total_tiles, (double)window_size / 1000, num_windows);

    /* Partial Fisher-Yates shuffle to select num_windows random tiles */
    int *all_tiles = malloc((size_t)total_tiles * sizeof(int));
    if (!all_tiles) { free(tiles_per_seq); free(tile_offset); return NULL; }
    for (int i = 0; i < total_tiles; i++)
        all_tiles[i] = i;

    srand(seed);
    for (int i = 0; i < num_windows; i++) {
        int j = i + (int)((unsigned long)rand() % (unsigned long)(total_tiles - i));
        int tmp = all_tiles[i];
        all_tiles[i] = all_tiles[j];
        all_tiles[j] = tmp;
    }

    /* Sort selected tiles for sequential memory access */
    qsort(all_tiles, (size_t)num_windows, sizeof(int), cmp_int_asc);

    /* Build segments from selected tiles */
    SeqSegment *segments = malloc((size_t)num_windows * sizeof(SeqSegment));
    if (!segments) { free(all_tiles); free(tiles_per_seq); free(tile_offset); return NULL; }

    for (int i = 0; i < num_windows; i++) {
        int tile_id = all_tiles[i];

        /* Find which sequence this tile belongs to */
        int seq_idx = 0;
        for (int s = n - 1; s >= 0; s--) {
            if (tile_id >= tile_offset[s]) {
                seq_idx = s;
                break;
            }
        }

        int tile_within_seq = tile_id - tile_offset[seq_idx];
        gpos_t seq_start = (seq_idx == 0) ? 0 : genome->boundaries[seq_idx - 1];

        segments[i].seq_index = seq_idx;
        segments[i].raw_start = seq_start + (gpos_t)tile_within_seq * window_size;
        segments[i].raw_end = segments[i].raw_start + window_size;
        segments[i].seg_length = window_size;
    }

    free(all_tiles);
    free(tiles_per_seq);
    free(tile_offset);

    /* Build sampled genome */
    Genome *sampled = genome_create_chunk(genome, segments, num_windows);

    if (sampled)
        fprintf(stderr, "Sampled genome: %.1fG → %.1fG bases (%d windows of %.0fkb, %.1f%%)\n",
                (double)genome->raw_length / 1e9,
                (double)sampled->raw_length / 1e9,
                num_windows, (double)window_size / 1000,
                100.0 * sampled->raw_length / genome->raw_length);

    *out_segments = segments;
    *out_num_segments = num_windows;
    return sampled;
}

/*
 * Remap instance coordinates from sampled genome back to original genome.
 * Updates both position and seq_index for each instance.
 */
static void remap_instance_coordinates(CandidateList *candidates,
                                        const SeqSegment *segments,
                                        int num_segments)
{
    /* Build cumulative start positions for segments in sample */
    gpos_t *seg_start = malloc((size_t)(num_segments + 1) * sizeof(gpos_t));
    seg_start[0] = 0;
    for (int s = 0; s < num_segments; s++)
        seg_start[s + 1] = seg_start[s] + segments[s].seg_length;

    int remapped = 0;
    for (int f = 0; f < candidates->num_families; f++) {
        CandidateFamily *fam = &candidates->families[f];
        if (!fam->instances) continue;

        for (int i = 0; i < fam->num_instances; i++) {
            gpos_t raw_pos = fam->instances[i].position - PADLENGTH;

            /* Find which segment this position falls in */
            int seg_idx = 0;
            for (int s = num_segments - 1; s >= 0; s--) {
                if (raw_pos >= seg_start[s]) {
                    seg_idx = s;
                    break;
                }
            }

            /* Remap position to original genome coordinates */
            gpos_t offset_in_seg = raw_pos - seg_start[seg_idx];
            gpos_t original_raw = segments[seg_idx].raw_start + offset_in_seg;
            fam->instances[i].position = original_raw + PADLENGTH;

            /* Update seq_index to original genome's sequence */
            fam->instances[i].seq_index = segments[seg_idx].seq_index;
            remapped++;
        }
    }

    free(seg_start);
    fprintf(stderr, "Remapped %d instance coordinates to original genome\n", remapped);
}

/*
 * Thread worker for parallel chunk discovery.
 */
typedef struct {
    const Genome        *full_genome;
    const SeqSegment    *segments;
    int                  num_segments;
    const DiscoverParams *params;
    CandidateList       *result;
    int                  chunk_id;
    int                  success;
} ChunkWorkerArgs;

static void *chunk_discover_worker(void *arg)
{
    ChunkWorkerArgs *w = (ChunkWorkerArgs *)arg;
    Genome *chunk_genome = genome_create_chunk(w->full_genome,
                                                w->segments, w->num_segments);
    if (!chunk_genome) {
        w->result = NULL;
        w->success = 0;
        return NULL;
    }
    /* Pass num_threads=1 inside chunk workers: chunks already run in
     * parallel at the chunk level, and kmer.c's global pool is not
     * safe under simultaneous parallel kmer_count calls (HASH_PORT_DESIGN §7.2). */
    w->result = discover_families(chunk_genome, w->params, 1);
    /* Remap instance positions from chunk-local to full-genome coordinates.
     * Without this, downstream output_bed produces negative BED start
     * coordinates because chunk positions are concatenated against the chunk
     * genome but final output subtracts full-genome chr_offset.
     * Discovered in V6_PHASE3_RESULT.md (TAIR10 nuclear had 5904 negative-start
     * BED rows on chr2). The sample-segments code path correctly calls
     * remap_instance_coordinates at main.c:1197; the chunked path missed it. */
    if (w->result)
        remap_instance_coordinates(w->result, w->segments, w->num_segments);
    genome_free_chunk(chunk_genome);
    w->success = (w->result != NULL);
    return NULL;
}

/*
 * Chunked discovery for large genomes.
 * Splits genome sequences into bins via first-fit-decreasing,
 * runs discover_families() on each chunk in parallel, and concatenates results.
 */
static CandidateList *discover_chunked(const Genome *genome,
                                       const DiscoverParams *params,
                                       int verbose, int num_threads,
                                       glen_t chunk_size)
{
    int n = genome->num_sequences;
    glen_t split_threshold = (glen_t)(chunk_size * 1.8);
    glen_t overlap = params->L;  /* overlap at each split point side */

    /* Each chunk computes its own l from its size (l = ceil(1+log4(N))).
     * With LPT balancing, all chunks are similar size → same l value.
     * Chunk-based l is more sensitive than full-genome l (one base shorter). */
    DiscoverParams chunk_params = *params;

    /* Build segment array: one segment per sequence, splitting long ones */
    int seg_cap = n + 16;
    int num_segments = 0;
    int num_splits = 0;
    SeqSegment *segments = malloc((size_t)seg_cap * sizeof(SeqSegment));
    if (!segments) return NULL;

    for (int i = 0; i < n; i++) {
        gpos_t seq_start = (i == 0) ? 0 : genome->boundaries[i - 1];
        gpos_t seq_end   = genome->boundaries[i];
        glen_t seq_len = seq_end - seq_start;
        if (i == genome->num_sequences - 1) seq_len -= 1;

        if (seq_len > split_threshold) {
            /* Recursively halve until each part <= split_threshold */
            int num_parts = 2;
            while (seq_len / num_parts > split_threshold)
                num_parts *= 2;

            glen_t part_base = seq_len / num_parts;

            fprintf(stderr, "  Splitting %s (%.1fM bases) into %d parts "
                    "(%.1fM each, %lldkb overlap)\n",
                    (genome->sequence_ids && genome->sequence_ids[i])
                        ? genome->sequence_ids[i] : "unnamed",
                    (double)seq_len / 1e6, num_parts,
                    (double)part_base / 1e6, (long long)overlap / 1000);

            /* Ensure capacity */
            while (num_segments + num_parts > seg_cap) {
                seg_cap *= 2;
                segments = realloc(segments, (size_t)seg_cap * sizeof(SeqSegment));
            }

            for (int j = 0; j < num_parts; j++) {
                gpos_t start = seq_start + j * part_base;
                gpos_t end   = (j == num_parts - 1)
                             ? (seq_start + seq_len)
                             : (seq_start + (j + 1) * part_base);

                /* Add overlap at split points */
                if (j > 0)              start -= overlap;
                if (j < num_parts - 1)  end   += overlap;

                /* Clamp to sequence boundaries */
                if (start < seq_start)          start = seq_start;
                if (end > seq_start + seq_len)  end   = seq_start + seq_len;

                segments[num_segments].seq_index  = i;
                segments[num_segments].raw_start  = start;
                segments[num_segments].raw_end    = end;
                segments[num_segments].seg_length = end - start;
                num_segments++;
            }
            num_splits++;
        } else {
            if (num_segments >= seg_cap) {
                seg_cap *= 2;
                segments = realloc(segments, (size_t)seg_cap * sizeof(SeqSegment));
            }
            segments[num_segments].seq_index  = i;
            segments[num_segments].raw_start  = seq_start;
            segments[num_segments].raw_end    = seq_start + seq_len;
            segments[num_segments].seg_length = seq_len;
            num_segments++;
        }
    }

    /* Sort segments by length descending (for LPT scheduling) */
    qsort(segments, (size_t)num_segments, sizeof(SeqSegment), cmp_seg_length_desc);

    /* Determine number of bins from total segment size (includes overlaps) */
    glen_t total_seg_size = 0;
    for (int i = 0; i < num_segments; i++)
        total_seg_size += segments[i].seg_length;

    int num_bins = (int)((total_seg_size + chunk_size - 1) / chunk_size);
    if (num_bins < 2) num_bins = 2;
    if (num_bins > num_segments) num_bins = num_segments;

    /* LPT (Longest Processing Time first) bin packing for balanced parallel load.
     * Each segment (sorted desc) is assigned to the bin with smallest current total.
     * This minimizes max(bin_size), optimal for parallel wall-clock time. */
    glen_t *bin_sizes = calloc((size_t)num_bins, sizeof(glen_t));
    SeqSegment **bin_segments = calloc((size_t)num_bins, sizeof(SeqSegment *));
    int *bin_counts = calloc((size_t)num_bins, sizeof(int));
    int *bin_caps = calloc((size_t)num_bins, sizeof(int));

    if (!bin_sizes || !bin_segments || !bin_counts || !bin_caps) {
        free(segments);
        free(bin_sizes);
        free(bin_segments);
        free(bin_counts);
        free(bin_caps);
        return NULL;
    }

    for (int b = 0; b < num_bins; b++) {
        bin_caps[b] = 16;
        bin_segments[b] = malloc((size_t)bin_caps[b] * sizeof(SeqSegment));
    }

    for (int i = 0; i < num_segments; i++) {
        /* Find bin with smallest current total (LPT) */
        int min_b = 0;
        for (int b = 1; b < num_bins; b++) {
            if (bin_sizes[b] < bin_sizes[min_b])
                min_b = b;
        }

        if (bin_counts[min_b] >= bin_caps[min_b]) {
            bin_caps[min_b] *= 2;
            bin_segments[min_b] = realloc(bin_segments[min_b],
                                           (size_t)bin_caps[min_b] * sizeof(SeqSegment));
        }
        bin_segments[min_b][bin_counts[min_b]++] = segments[i];
        bin_sizes[min_b] += segments[i].seg_length;
    }

    free(segments);

    int max_parallel = num_bins < num_threads ? num_bins : num_threads;
    fprintf(stderr, "Large genome (%.1fG bases): %d segments (%d sequences split) → "
            "%d chunks (%d parallel)\n",
            (double)genome->raw_length / 1e9, num_segments, num_splits,
            num_bins, max_parallel);

    for (int b = 0; b < num_bins; b++)
        fprintf(stderr, "  Chunk %d/%d: %d segments, %.1fM bases\n",
                b + 1, num_bins, bin_counts[b], (double)bin_sizes[b] / 1e6);

    /* Run discovery on chunks in parallel (batches of max_parallel) */
    int total_families = 0;
    CandidateList **chunk_results = calloc((size_t)num_bins, sizeof(CandidateList *));
    ChunkWorkerArgs *workers = calloc((size_t)num_bins, sizeof(ChunkWorkerArgs));
    pthread_t *threads = malloc((size_t)max_parallel * sizeof(pthread_t));

    if (!chunk_results || !workers || !threads) {
        free(chunk_results); free(workers); free(threads);
        goto cleanup_bins;
    }

    for (int batch_start = 0; batch_start < num_bins; batch_start += max_parallel) {
        int batch_end = batch_start + max_parallel;
        if (batch_end > num_bins) batch_end = num_bins;
        int batch_size = batch_end - batch_start;

        /* Launch batch */
        for (int i = 0; i < batch_size; i++) {
            int b = batch_start + i;
            workers[b].full_genome = genome;
            workers[b].segments = bin_segments[b];
            workers[b].num_segments = bin_counts[b];
            workers[b].params = &chunk_params;
            workers[b].result = NULL;
            workers[b].chunk_id = b;
            workers[b].success = 0;
            pthread_create(&threads[i], NULL, chunk_discover_worker, &workers[b]);
        }

        /* Wait for batch */
        for (int i = 0; i < batch_size; i++)
            pthread_join(threads[i], NULL);

        /* Check results */
        for (int i = 0; i < batch_size; i++) {
            int b = batch_start + i;
            if (!workers[b].success) {
                fprintf(stderr, "ERROR: Discovery failed on chunk %d\n", b + 1);
                for (int j = 0; j < num_bins; j++)
                    if (chunk_results[j]) candidates_free(chunk_results[j]);
                free(chunk_results); free(workers); free(threads);
                goto cleanup_bins;
            }
            chunk_results[b] = workers[b].result;
            total_families += chunk_results[b]->num_families;
            if (verbose)
                fprintf(stderr, "  Chunk %d: %d families discovered\n",
                        b + 1, chunk_results[b]->num_families);
        }
    }

    free(workers);
    free(threads);

    fprintf(stderr, "Total candidates from all chunks: %d\n", total_families);

    /* Concatenate all chunk CandidateLists into one */
    CandidateList *combined = malloc(sizeof(CandidateList));
    if (!combined) goto cleanup_bins;

    combined->num_families = total_families;
    combined->cap_families = total_families > 0 ? total_families : 1;
    combined->families = malloc((size_t)combined->cap_families * sizeof(CandidateFamily));
    if (!combined->families) {
        free(combined);
        combined = NULL;
        goto cleanup_bins;
    }

    int write_idx = 0;
    for (int b = 0; b < num_bins; b++) {
        CandidateList *cl = chunk_results[b];
        if (!cl) continue;
        for (int f = 0; f < cl->num_families; f++) {
            combined->families[write_idx] = cl->families[f];
            combined->families[write_idx].id = (mdl_uid_t)write_idx;
            write_idx++;
        }
        /* Free the chunk CandidateList shell (but NOT the families' data,
         * which is now owned by the combined list) */
        free(cl->families);
        cl->families = NULL;
        cl->num_families = 0;
        candidates_free(cl);
    }

    free(chunk_results);

    /* Clean up bin-packing arrays */
    for (int b = 0; b < num_bins; b++)
        free(bin_segments[b]);
    free(bin_segments);
    free(bin_sizes);
    free(bin_counts);
    free(bin_caps);

    return combined;

cleanup_bins:
    free(chunk_results);
    for (int b = 0; b < num_bins; b++)
        free(bin_segments[b]);
    free(bin_segments);
    free(bin_sizes);
    free(bin_counts);
    free(bin_caps);
    return NULL;
}

static void usage(void)
{
    fprintf(stderr,
        "mdl-repeat: MDL-guided seed-and-extend de novo repeat discovery\n\n"
        "Usage:\n"
        "  mdl-repeat -sequence <file> -output <file> [options]\n\n"
        "Required:\n"
        "  -sequence <file>   Input FASTA\n"
        "  -output <file>     Output repeat library (FASTA)\n\n"
        "Optional (discovery):\n"
        "  -freq <file>       Pre-computed l-mer frequency table (build_lmer_table format)\n"
        "  -freq-output <file> Write l-mer frequency table for reuse\n"
        "  -l #               L-mer length (default: auto = ceil(1+log4(N)))\n"
        "  -L #               Max extension distance per side (default: 10000)\n"
        "  -minthresh #       Min l-mer frequency to seed (default: 2)\n"
        "  -goodlength #      Min consensus length pre-filter (default: 30)\n"
        "  -maxgap #          Max DP band offset (default: 5)\n"
        "  -match #           Match score (default: 1)\n"
        "  -mismatch #        Mismatch score (default: -1)\n"
        "  -gap #             Gap penalty (default: -5)\n"
        "  -cappenalty #      Cap penalty for exiting alignment (default: -20)\n"
        "  -minimprovement #  Min improvement per step (default: 3)\n"
        "  -stopafter #       Stop after N no-progress positions (default: 100)\n"
        "  -maxentropy #      Entropy filter threshold (default: -0.70)\n"
        "  -tandemdist #      Min distance between same-strand l-mers (default: 500)\n"
        "  -maxoccurrences #  Max occurrences per seed (default: 10000)\n"
        "  -maxrepeats #      Max families to discover (default: 100000)\n"
        "  -recall-rescue    Run bounded secondary discovery with shorter l-mer seeds\n"
        "  -recruit-short    Recruit extra instances for short families (<500 bp) via\n"
        "                     rmblastn before selection (needs rmblastn in PATH or\n"
        "                     $RMBLASTN_BIN); rescues divergent SINE/MITE families\n"
        "  -rescue-full-genome Run rescue discovery on the full genome instead of gaps\n"
        "  -rescue-l-delta # Shorten rescue l-mer length by N (default: 1, min l=8)\n"
        "  -rescue-maxrepeats # Max rescue families to append (default: 2000)\n"
        "  -rescue-min-gap # Min uncovered gap length for targeted rescue (default: 200)\n"
        "  -chunk-size #      Chunk size in Mb for large genome discovery (default: 200)\n"
        "  -sample-size #     Sample size in Mb for large genomes (default: 1000)\n"
        "  -sample-output <file> Write sampled genome FASTA for reproducibility\n"
        "  -window-size #     Sampling window size in kb (default: 1000)\n"
        "  -seed #            Random seed for sampling (default: 42)\n\n"
        "Optional (refinement):\n"
        "  -threads #         Number of threads for refinement (default: 1)\n"
        "  -mdl-mode <mode>   MDL position encoding: none|exact|upper (default: exact)\n"
        "  -max-divergence #  Max substitution rate for instance acceptance (default: 0.30)\n"
        "  -refine-gap #      Refinement gap penalty (default: -5; try -3 for high-indel species)\n"
        "  -refine-maxoffset # Refinement DP band half-width (default: 12, max: 32)\n"
        "  -max-dp-cells #    Max DP cells for consensus merge (default: 10000000)\n"
        "  -coalesce-factor # Gap tolerance for tandem-instance coalescing, in bases\n"
        "                     (default: 20.0; 0 = disabled)\n"
        "  -max-instances #   Per-family recruited-instance cap (default: 10000).\n"
        "                     Raise for genomes with extreme-copy families (e.g. Alu);\n"
        "                     low caps truncate reported copies= and MDL savings.\n"
        "  -max-seed-hits #   Raw seed-hit cap per family (default: 50000; auto-raised\n"
        "                     to >= -max-instances)\n\n"
        "Optional (output):\n"
        "  -instances <file>  Output instance BED\n"
        "  -stats <file>      Output family statistics TSV\n"
        "  -trace-dir <dir>   Dump FASTA+BED+TSV after each refine stage (diagnostics)\n"
        "  -split-audit <file> Write split-stage decision audit TSV\n"
        "  -rescue-audit <file> Write recall-rescue target/candidate audit TSV\n"
        "  -external-tools <mode> External tool policy: off|auto|require (default: off)\n"
        "  -external-qc <file> Write optional seqkit stats TSV for final FASTA\n"
        "  -seqkit <path>      Path to seqkit executable for -external-qc\n"
        "  -v / -vv           Verbosity level\n"
    );
    exit(1);
}

int main(int argc, char *argv[])
{
    char *sequence_file = NULL;
    char *output_file = NULL;
    char *instances_file = NULL;
    char *stats_file = NULL;
    char *freq_file = NULL;
    char *freq_output = NULL;
    char *trace_dir = NULL;
    char *split_audit_file = NULL;
    char *rescue_audit_file = NULL;
    char *external_tools_mode_str = NULL;
    char *external_qc_file = NULL;
    char *seqkit_path = NULL;
    int num_threads = 1;
    int verbose = 0;
    float coalesce_factor = 20.0f;
    ExternalToolsMode external_tools_mode = EXTERNAL_TOOLS_OFF;
    int recall_rescue = 0;
    int recruit_short = 0;
    int rescue_full_genome = 0;
    int rescue_l_delta = 1;
    int rescue_maxrepeats = 2000;
    int rescue_min_gap = 200;
    int x;

    /* Parse required arguments */
    if (!co_get_string(argc, argv, "-sequence", &sequence_file) ||
        !co_get_string(argc, argv, "-output", &output_file)) {
        usage();
    }

    /* Parse optional arguments */
    co_get_string(argc, argv, "-instances", &instances_file);
    co_get_string(argc, argv, "-stats", &stats_file);
    co_get_string(argc, argv, "-freq-output", &freq_output);
    co_get_string(argc, argv, "-freq", &freq_file);
    co_get_string(argc, argv, "-trace-dir", &trace_dir);
    co_get_string(argc, argv, "-split-audit", &split_audit_file);
    co_get_string(argc, argv, "-rescue-audit", &rescue_audit_file);
    co_get_string(argc, argv, "-external-tools", &external_tools_mode_str);
    co_get_string(argc, argv, "-external-qc", &external_qc_file);
    co_get_string(argc, argv, "-seqkit", &seqkit_path);
    if (co_has_option(argc, argv, "-split-audit") && !split_audit_file) {
        fprintf(stderr, "ERROR: -split-audit requires a file path\n");
        return 1;
    }
    if (co_has_option(argc, argv, "-rescue-audit") && !rescue_audit_file) {
        fprintf(stderr, "ERROR: -rescue-audit requires a file path\n");
        return 1;
    }
    if (co_has_option(argc, argv, "-external-tools") &&
        !external_tools_mode_str) {
        fprintf(stderr, "ERROR: -external-tools requires off|auto|require\n");
        return 1;
    }
    if (external_tools_mode_str &&
        !external_tools_mode_parse(external_tools_mode_str,
                                   &external_tools_mode)) {
        fprintf(stderr, "ERROR: -external-tools must be off, auto, or require\n");
        return 1;
    }
    if (co_has_option(argc, argv, "-external-qc") && !external_qc_file) {
        fprintf(stderr, "ERROR: -external-qc requires a file path\n");
        return 1;
    }
    if (co_has_option(argc, argv, "-seqkit") && !seqkit_path) {
        fprintf(stderr, "ERROR: -seqkit requires a file path\n");
        return 1;
    }
    if (external_qc_file && external_tools_mode == EXTERNAL_TOOLS_OFF) {
        external_tools_mode = EXTERNAL_TOOLS_AUTO;
    }
    g_refine_split_audit_path = split_audit_file;
    co_get_bool(argc, argv, "-recall-rescue", &recall_rescue);
    co_get_bool(argc, argv, "-recruit-short", &recruit_short);
    co_get_bool(argc, argv, "-rescue-full-genome", &rescue_full_genome);
    if (co_has_option(argc, argv, "-rescue-l-delta") &&
        !co_get_int(argc, argv, "-rescue-l-delta", &rescue_l_delta)) {
        fprintf(stderr, "ERROR: -rescue-l-delta requires an integer\n");
        return 1;
    }
    if (co_get_int(argc, argv, "-rescue-l-delta", &rescue_l_delta)) {
        if (rescue_l_delta < 0 || rescue_l_delta > 8) {
            fprintf(stderr, "ERROR: -rescue-l-delta must be in [0, 8]\n");
            return 1;
        }
    }
    if (co_has_option(argc, argv, "-rescue-maxrepeats") &&
        !co_get_int(argc, argv, "-rescue-maxrepeats", &rescue_maxrepeats)) {
        fprintf(stderr, "ERROR: -rescue-maxrepeats requires an integer\n");
        return 1;
    }
    if (co_get_int(argc, argv, "-rescue-maxrepeats", &rescue_maxrepeats)) {
        if (rescue_maxrepeats < 1) {
            fprintf(stderr, "ERROR: -rescue-maxrepeats must be positive\n");
            return 1;
        }
    }
    if (co_has_option(argc, argv, "-rescue-min-gap") &&
        !co_get_int(argc, argv, "-rescue-min-gap", &rescue_min_gap)) {
        fprintf(stderr, "ERROR: -rescue-min-gap requires an integer\n");
        return 1;
    }
    if (co_get_int(argc, argv, "-rescue-min-gap", &rescue_min_gap)) {
        if (rescue_min_gap < 1) {
            fprintf(stderr, "ERROR: -rescue-min-gap must be positive\n");
            return 1;
        }
    }
    if (trace_dir) {
        /* Create trace directory if it doesn't exist */
        if (mkdir(trace_dir, 0755) != 0 && errno != EEXIST) {
            fprintf(stderr, "WARNING: could not create trace-dir '%s'\n", trace_dir);
            trace_dir = NULL;
        } else {
            fprintf(stderr, "Trace dumps enabled → %s\n", trace_dir);
        }
    }
    co_get_int(argc, argv, "-threads", &num_threads);
    if (num_threads < 1) num_threads = 1;

    glen_t chunk_size = DISCOVER_SPLIT_THRESHOLD_DEFAULT;
    {
        int chunk_mb = 0;
        if (co_get_int(argc, argv, "-chunk-size", &chunk_mb)) {
            if (chunk_mb < 10) {
                fprintf(stderr, "ERROR: -chunk-size must be >= 10 (Mb)\n");
                return 1;
            }
            chunk_size = (glen_t)chunk_mb * 1000000LL;
        }
    }

    /* Sampling parameters for very large genomes */
    glen_t sample_size = 1000LL * 1000000LL;    /* default 1Gb */
    glen_t window_size = 1000LL * 1000LL;       /* default 1Mb */
    unsigned int sample_seed = 42;
    char *sample_output = NULL;
    co_get_string(argc, argv, "-sample-output", &sample_output);
    {
        int sample_mb = 0;
        if (co_get_int(argc, argv, "-sample-size", &sample_mb)) {
            if (sample_mb < 100) {
                fprintf(stderr, "ERROR: -sample-size must be >= 100 (Mb)\n");
                return 1;
            }
            sample_size = (glen_t)sample_mb * 1000000LL;
        }
        int window_kb = 0;
        if (co_get_int(argc, argv, "-window-size", &window_kb)) {
            if (window_kb < 100 || window_kb > 10000) {
                fprintf(stderr, "ERROR: -window-size must be in [100, 10000] (kb)\n");
                return 1;
            }
            window_size = (glen_t)window_kb * 1000LL;
        }
        int seed_val = 0;
        if (co_get_int(argc, argv, "-seed", &seed_val))
            sample_seed = (unsigned int)seed_val;
    }

    /* MDL mode selection */
    char *mdl_mode_str = NULL;
    co_get_string(argc, argv, "-mdl-mode", &mdl_mode_str);
    if (mdl_mode_str) {
        if (strcmp(mdl_mode_str, "none") == 0)
            g_mdl_mode = MDL_MODE_NONE;
        else if (strcmp(mdl_mode_str, "exact") == 0)
            g_mdl_mode = MDL_MODE_EXACT;
        else if (strcmp(mdl_mode_str, "upper") == 0)
            g_mdl_mode = MDL_MODE_UPPER;
        else {
            fprintf(stderr, "ERROR: unknown -mdl-mode '%s' (use none|exact|upper)\n",
                    mdl_mode_str);
            return 1;
        }
    }

    verbose = co_get_bool(argc, argv, "-vv", &x)  ? 2 :
              co_get_bool(argc, argv, "-v", &x)    ? 1 : 0;

    /* ================================================================
     * Step 1: Load genome
     * ================================================================ */
    if (verbose) fprintf(stderr, "Loading genome from %s...\n", sequence_file);
    Genome *genome = genome_load(sequence_file);
    if (genome == NULL) {
        fprintf(stderr, "ERROR: Failed to load genome\n");
        return 1;
    }
    if (verbose) genome_print_stats(genome);

    /* Track full genome metadata for BED remapping and display */
    Genome *full_genome = genome;
    SeqSegment *sample_segments = NULL;
    int num_sample_segments = 0;

    /* ================================================================
     * Step 1b: Sample genome if too large
     * ================================================================ */
    if (genome->raw_length > sample_size) {
        fprintf(stderr, "Genome (%.1fG bases) exceeds sample-size (%.1fG), activating sampling...\n",
                (double)genome->raw_length / 1e9, (double)sample_size / 1e9);

        Genome *sampled = genome_sample_windows(genome, sample_size, window_size,
                                                 &sample_segments, &num_sample_segments,
                                                 sample_seed, verbose);
        if (!sampled) {
            fprintf(stderr, "ERROR: Genome sampling failed\n");
            genome_free(genome);
            return 1;
        }

        /* Write sampled genome to FASTA if requested */
        if (sample_output)
            genome_write_fasta(sampled, sample_output, sample_segments,
                               num_sample_segments);

        /* Free original sequence buffer to reclaim memory (keep metadata) */
        free(genome->sequence);
        genome->sequence = NULL;

        genome = sampled;  /* pipeline uses sampled genome from here */
    }

    int k = default_k(genome->length);
    fprintf(stderr, "Using k=%d, mdl-mode=%s\n", k,
            g_mdl_mode == MDL_MODE_NONE  ? "none" :
            g_mdl_mode == MDL_MODE_EXACT ? "exact" : "upper");

    /* ================================================================
     * Step 2: Discover consensus families (seed-and-extend)
     * ================================================================ */
    DiscoverParams dparams;
    discover_default_params(&dparams);
    dparams.freq_file = freq_file;
    dparams.freq_output = freq_output;
    dparams.VERBOSE = verbose;

    /* Override discovery parameters from CLI */
    int tmp_int;
    float tmp_float;
    if (co_get_int(argc, argv, "-l", &tmp_int))          dparams.l = tmp_int;
    if (co_get_int(argc, argv, "-L", &tmp_int))          dparams.L = tmp_int;
    if (co_get_int(argc, argv, "-minthresh", &tmp_int))   dparams.MINTHRESH = tmp_int;
    if (co_get_int(argc, argv, "-goodlength", &tmp_int))  dparams.GOODLENGTH = tmp_int;
    if (co_get_int(argc, argv, "-maxgap", &tmp_int))      dparams.MAXOFFSET = tmp_int;
    if (co_get_int(argc, argv, "-match", &tmp_int))        dparams.MATCH = tmp_int;
    if (co_get_int(argc, argv, "-mismatch", &tmp_int))     dparams.MISMATCH = tmp_int;
    if (co_get_int(argc, argv, "-gap", &tmp_int))          dparams.GAP = tmp_int;
    if (co_get_int(argc, argv, "-cappenalty", &tmp_int))    dparams.CAPPENALTY = tmp_int;
    if (co_get_int(argc, argv, "-minimprovement", &tmp_int)) dparams.MINIMPROVEMENT = tmp_int;
    if (co_get_int(argc, argv, "-stopafter", &tmp_int))     dparams.WHEN_TO_STOP = tmp_int;
    if (co_get_float(argc, argv, "-maxentropy", &tmp_float)) dparams.MAXENTROPY = tmp_float;
    if (co_get_int(argc, argv, "-tandemdist", &tmp_int))    dparams.TANDEMDIST = tmp_int;
    if (co_get_int(argc, argv, "-maxoccurrences", &tmp_int)) dparams.MAXN = tmp_int;
    if (co_get_int(argc, argv, "-maxrepeats", &tmp_int))     dparams.MAXR = tmp_int;

    /* Override refinement-stage parameters from CLI */
    if (co_get_float(argc, argv, "-max-divergence", &tmp_float)) {
        if (tmp_float < 0.0f || tmp_float > 1.0f) {
            fprintf(stderr, "ERROR: -max-divergence must be in [0, 1]\n");
            return 1;
        }
        g_align_max_divergence = tmp_float;
    }
    if (co_get_int(argc, argv, "-refine-gap", &tmp_int))
        g_align_gap = tmp_int;
    if (co_get_int(argc, argv, "-max-instances", &tmp_int)) {
        if (tmp_int < 2) {
            fprintf(stderr, "ERROR: -max-instances must be >= 2\n");
            return 1;
        }
        g_align_max_instances = tmp_int;
    }
    if (co_get_int(argc, argv, "-max-seed-hits", &tmp_int)) {
        if (tmp_int < 2) {
            fprintf(stderr, "ERROR: -max-seed-hits must be >= 2\n");
            return 1;
        }
        g_align_max_seed_hits = tmp_int;
    }
    /* Keep the seed-hit ceiling at least as large as the instance ceiling,
     * otherwise raising -max-instances alone is throttled upstream. */
    if (g_align_max_seed_hits < g_align_max_instances)
        g_align_max_seed_hits = g_align_max_instances;
    if (co_get_int(argc, argv, "-refine-maxoffset", &tmp_int)) {
        if (tmp_int < 1 || tmp_int > ALIGN_MAXOFFSET_LIMIT) {
            fprintf(stderr, "ERROR: -refine-maxoffset must be in [1, %d]\n",
                    ALIGN_MAXOFFSET_LIMIT);
            return 1;
        }
        g_align_maxoffset = tmp_int;
    }
    {
        int dp_cells_m = 0;
        if (co_get_int(argc, argv, "-max-dp-cells", &dp_cells_m)) {
            if (dp_cells_m < 1) {
                fprintf(stderr, "ERROR: -max-dp-cells must be positive\n");
                return 1;
            }
            g_refine_max_dp_cells = (int64_t)dp_cells_m;
        }
    }
    if (co_get_float(argc, argv, "-coalesce-factor", &tmp_float)) {
        if (tmp_float < 0.0f) {
            fprintf(stderr, "ERROR: -coalesce-factor must be >= 0 (use 0 to disable)\n");
            return 1;
        }
        coalesce_factor = tmp_float;
    }

    fprintf(stderr, "Discovering repeat families...\n");
    CandidateList *candidates;
    if (genome->raw_length > chunk_size) {
        candidates = discover_chunked(genome, &dparams, verbose, num_threads, chunk_size);
    } else {
        candidates = discover_families(genome, &dparams, num_threads);
    }
    if (candidates == NULL) {
        fprintf(stderr, "ERROR: Discovery failed\n");
        if (sample_segments) {
            free(sample_segments);
            genome_free_chunk(genome);
            genome_free(full_genome);
        } else {
            genome_free(genome);
        }
        return 1;
    }

    if (recall_rescue) {
        RecallRescueOptions rescue_options;
        rescue_options.full_genome = rescue_full_genome;
        rescue_options.l_delta = rescue_l_delta;
        rescue_options.max_repeats = rescue_maxrepeats;
        rescue_options.min_gap = rescue_min_gap;
        rescue_options.verbose = verbose;
        rescue_options.num_threads = num_threads;
        rescue_options.chunk_size = chunk_size;
        rescue_options.audit_file = rescue_audit_file;

        RecallRescueCallbacks rescue_callbacks;
        rescue_callbacks.create_chunk = genome_create_chunk;
        rescue_callbacks.free_chunk = genome_free_chunk;
        rescue_callbacks.discover_chunked = discover_chunked;
        rescue_callbacks.remap_instance_coordinates = remap_instance_coordinates;

        if (recall_rescue_run(genome, candidates, &dparams,
                              &rescue_options, &rescue_callbacks) < 0) {
            candidates_free(candidates);
            if (sample_segments) {
                free(sample_segments);
                genome_free_chunk(genome);
                genome_free(full_genome);
            } else {
                genome_free(genome);
            }
            return 1;
        }
    }

    if (verbose) candidates_print_stats(candidates);
    /* Trace stage 01: discover */
    if (trace_dir) trace_dump(trace_dir, 1, "discover", candidates, genome);

    /* ================================================================
     * Step 3: Build k-mer table (needed for refine pipeline)
     * ================================================================ */
    if (verbose) fprintf(stderr, "Building k-mer table for refinement (k=%d)...\n", k);
    KmerTable *kt = kmer_count(genome, k, DEFAULT_TANDEMDIST, num_threads);
    if (kt == NULL) {
        fprintf(stderr, "ERROR: K-mer counting failed\n");
        candidates_free(candidates);
        genome_free(genome);
        return 1;
    }

    /* Trim at MINTHRESH=2 to match discovery threshold */
    int kmer_minthresh = dparams.MINTHRESH;
    if (kmer_minthresh < 2) kmer_minthresh = 2;
    int64_t kmer_removed = kmer_trim(kt, kmer_minthresh);
    fprintf(stderr, "Trimmed %" PRId64 " k-mers with freq < %d, %" PRId64 " remain\n",
            kmer_removed, kmer_minthresh, kt->num_entries);

    /* ================================================================
     * Step 3b: Build k-mer position index (needed by refine/merge)
     * ================================================================ */
    if (verbose) fprintf(stderr, "Building k-mer position index...\n");
    kmer_build_positions(kt, genome, num_threads);

    /* ================================================================
     * Step 4: Compact — remove families with < 2 instances or short consensus
     * ================================================================ */
    {
        int write_idx = 0;
        int compact_removed = 0;
        int min_final_consensus = (k * 2 > dparams.GOODLENGTH) ? k * 2 : dparams.GOODLENGTH;

        for (int i = 0; i < candidates->num_families; i++) {
            CandidateFamily *f = &candidates->families[i];
            if (f->num_instances >= 2 && f->consensus_length >= min_final_consensus) {
                if (write_idx != i)
                    candidates->families[write_idx] = candidates->families[i];
                candidates->families[write_idx].id = (mdl_uid_t)write_idx;
                write_idx++;
            } else {
                free(f->consensus);
                free(f->instances);
                f->consensus = NULL;
                f->instances = NULL;
                compact_removed++;
            }
        }
        candidates->num_families = write_idx;
        fprintf(stderr, "Compacted: removed %d dead/short families, %d remain\n",
                compact_removed, write_idx);
    }
    /* Trace stage 02: compact */
    if (trace_dir) trace_dump(trace_dir, 2, "compact", candidates, genome);

    if (candidates->num_families == 0) {
        fprintf(stderr, "No families survived compaction. Writing empty output.\n");
        output_fasta(output_file, candidates);
        candidates_free(candidates);
        kmer_free(kt);
        if (sample_segments) {
            free(sample_segments);
            genome_free_chunk(genome);
            genome_free(full_genome);
        } else {
            genome_free(genome);
        }
        return 0;
    }

    if (verbose) candidates_print_stats(candidates);

    /* Phase 5 (F') REVERTED 2026-05-01: Step 4b (align_refine_all pre-pass)
     * caused 12x instance explosion that made refine_assemble_fragments
     * O(n²) hang for 27+ hours on TAIR10 nuclear (single-thread, 11 GB RSS,
     * no completion). Step 4c (RMBlast short-family recruit) was bundled with
     * 4b and is reverted together. chr4 family-level recall also regressed
     * (-7.5pp 80×80, -60pp on 10-99 copy bin).
     * Functions align_blast_recruit_short_families and the early
     * align_refine_all call are kept in source for reference but unwired. */

    /* ================================================================
     * Step 5: Merge redundant families
     * ================================================================ */
    if (verbose) fprintf(stderr, "Merging redundant families...\n");
    int n_merges = refine_merge_families(candidates, genome, kt, k, verbose, num_threads);
    fprintf(stderr, "Merged %d families, %d remain\n",
            n_merges, candidates->num_families);
    if (verbose && n_merges > 0) candidates_print_stats(candidates);
    /* Trace stage 03: merge */
    if (trace_dir) trace_dump(trace_dir, 3, "merge", candidates, genome);

    /* ================================================================
     * Step 6: Split families with bimodal divergence
     * ================================================================ */
    if (verbose) fprintf(stderr, "Splitting bimodal families...\n");
    /* Single split pass. An iterated 4-pass variant was benchmarked and
     * rejected: it produced no family-level recall gain while inflating the
     * library ~22% and pushing the compression ratio above 1 (dl_library >
     * total_savings). Cause: each extra Otsu pass mints overlapping sub-families
     * that clear the relaxed split MDL gate, then enter the library via the
     * standalone-fallback admit branch (mdl.c) carrying model_cost but ~zero
     * exclusive savings. Keep this single-pass. */
    int n_splits = refine_split_families(candidates, genome, kt, k,
                                         genome->length, verbose,
                                         candidates->num_families,
                                         num_threads);
    fprintf(stderr, "Split %d families, %d total\n",
            n_splits, candidates->num_families);
    if (verbose && n_splits > 0) candidates_print_stats(candidates);
    /* Trace stage 04: split */
    if (trace_dir) trace_dump(trace_dir, 4, "split", candidates, genome);

    /* ================================================================
     * Step 6b: Fragment assembly (assemble adjacent TE fragments)
     * ================================================================ */
    if (verbose) fprintf(stderr, "Assembling TE fragments...\n");
    int n_assembled = refine_assemble_fragments(candidates, genome, kt, k,
                                                 genome->length, verbose,
                                                 num_threads);
    fprintf(stderr, "Assembled %d fragment pairs, %d families remain\n",
            n_assembled, candidates->num_families);
    if (verbose && n_assembled > 0) candidates_print_stats(candidates);
    /* Trace stage 05: assemble */
    if (trace_dir) trace_dump(trace_dir, 5, "assemble", candidates, genome);

    /* ================================================================
     * Step 6c: Short-family RMBlast recruitment (opt-in, needs rmblastn)
     * Banded-DP recruitment under-recruits divergent copies of SHORT elements
     * (SINEs/MITEs/short fragments, <500 bp) — the dominant recall gap on
     * TAIR10 (<500 bp families: 36% missed at 80x80 vs 4% for >=6 kb). rmblastn
     * is far more sensitive for short divergent matches; this batch-recruits
     * extra non-overlapping instances (div<=max-divergence gated, capped per
     * family) so weakly-supported real short families clear the MDL/standalone
     * admission gates instead of being pruned. Placed AFTER assembly so the
     * boosted instance counts never feed the assembly sweep (the 12x-explosion
     * that caused the 2026-05-01 revert lived in the bundled align_refine_all
     * pre-pass, not here; the assembly SWEEP_BUDGET guard now backstops it too).
     * ================================================================ */
    if (recruit_short) {
        if (verbose) fprintf(stderr, "Short-family RMBlast recruitment...\n");
        int n_recruited = align_blast_recruit_short_families(candidates, genome,
                                                             k, num_threads, verbose);
        if (n_recruited < 0)
            fprintf(stderr, "WARNING: -recruit-short requested but rmblastn "
                    "not available (set $RMBLASTN_BIN or add to PATH); skipped.\n");
        else
            fprintf(stderr, "Short-family recruitment: added %d instances\n",
                    n_recruited);
    }

    /* ================================================================
     * Step 7: MDL scoring and library selection
     * ================================================================ */
    if (verbose) fprintf(stderr, "MDL scoring and library selection...\n");
    MDLResult mdl_result = mdl_select_library(candidates, genome->length);

    fprintf(stderr, "MDL Results:\n");
    fprintf(stderr, "  Accepted families:    %d / %d\n",
            mdl_result.num_accepted, candidates->num_families);
    fprintf(stderr, "  Bases covered:        %" PRId64 " / %" PRId64 " (%.1f%%)\n",
            mdl_result.bases_covered, (int64_t)genome->length,
            100.0 * mdl_result.bases_covered / genome->length);
    fprintf(stderr, "  Library cost:         %.0f bits\n", mdl_result.dl_library);
    fprintf(stderr, "  Total savings:        %.0f bits\n", mdl_result.total_savings);
    fprintf(stderr, "  Compression ratio:    %.4f\n", mdl_result.compression_ratio);
    /* Trace stage 06: mdl_select — only accepted families */
    if (trace_dir) trace_dump_filtered(trace_dir, 6, "mdl_select", candidates, genome, 1);

    /* ================================================================
     * Step 8: Prune marginal families
     * ================================================================ */
    int n_pruned = refine_prune_families(candidates, genome->length, verbose,
                                          mdl_result.num_accepted);
    if (n_pruned > 0)
        fprintf(stderr, "Pruned %d marginal families\n", n_pruned);
    /* Trace stage 07: prune — only accepted families */
    if (trace_dir) trace_dump_filtered(trace_dir, 7, "prune", candidates, genome, 1);

    /* ================================================================
     * Step 8b: Tandem-instance coalescing
     * Pure reporting transform — merges adjacent same-family same-strand
     * instances within 1.5 × consensus_length.  Aligns mdl-repeat output
     * with RM-style truth annotation convention (which merges adjacent
     * tandem copies into one interval).  No effect on MDL scoring (done).
     * ================================================================ */
    int n_coalesced = 0;
    if (coalesce_factor > 0.0f) {
        n_coalesced = refine_coalesce_tandem_instances(candidates, coalesce_factor, verbose);
        if (n_coalesced > 0)
            fprintf(stderr, "Coalesced %d tandem instance pairs\n", n_coalesced);
    }
    /* Trace stage 08: coalesce — only accepted families */
    if (trace_dir) trace_dump_filtered(trace_dir, 8, "coalesce", candidates, genome, 1);

    /* ================================================================
     * Step 8c: Drop chimeric / over-extended long families
     * Long + high-divergence consensi (assembly chains, discovery
     * over-extensions to the L-cap, coalesced tandems) whose copies match only
     * a fragment.  Runs last so it catches monsters from every upstream stage.
     * Length is scope; divergence is the decision.  Recall-neutral.
     * ================================================================ */
    int n_chimeric = refine_drop_chimeric_long(candidates, verbose);
    if (n_chimeric > 0)
        fprintf(stderr, "Dropped %d chimeric/over-extended long families\n",
                n_chimeric);

    /* #1: warn if any accepted family is pinned at the instance ceiling — its
     * reported copies= and MDL savings are truncated (the consensus itself is
     * fine; it was built during discovery). This counts the precise condition
     * (num_instances == cap), not mere seed-hit saturation. */
    {
        int n_capped = 0;
        for (int fi = 0; fi < candidates->num_families; fi++)
            if (candidates->families[fi].num_instances >= g_align_max_instances)
                n_capped++;
        if (n_capped > 0)
            fprintf(stderr,
                    "WARNING: %d famil%s at the instance cap "
                    "(-max-instances=%d); their copies= and MDL savings are "
                    "truncated. Raise -max-instances for extreme-copy families.\n",
                    n_capped, n_capped == 1 ? "y is" : "ies are",
                    g_align_max_instances);
    }

    /* ================================================================
     * Step 9: Write outputs
     * ================================================================ */
    output_fasta(output_file, candidates);

    if (instances_file) {
        if (sample_segments) {
            /* Remap instance coordinates from sample → original genome */
            remap_instance_coordinates(candidates, sample_segments, num_sample_segments);
            output_bed(instances_file, candidates, full_genome);
        } else {
            output_bed(instances_file, candidates, genome);
        }
    }

    if (stats_file)
        output_stats(stats_file, candidates);

    if (external_qc_file) {
        ExternalQcConfig external_qc;
        memset(&external_qc, 0, sizeof(external_qc));
        external_qc.mode = external_tools_mode;
        external_qc.seqkit_path = seqkit_path;
        external_qc.qc_output_path = external_qc_file;
        external_qc.timeout_sec = 300;
        if (external_qc_run_seqkit_stats(&external_qc, output_file) != 0) {
            candidates_free(candidates);
            kmer_free(kt);
            if (sample_segments) {
                free(sample_segments);
                genome_free_chunk(genome);
                genome_free(full_genome);
            } else {
                genome_free(genome);
            }
            return 1;
        }
    }

    /* Cleanup */
    candidates_free(candidates);
    kmer_free(kt);
    if (sample_segments) {
        free(sample_segments);
        genome_free_chunk(genome);   /* free sampled genome (shares strings) */
        genome_free(full_genome);    /* free original metadata + strings */
    } else {
        genome_free(genome);
    }
    fprintf(stderr, "Done.\n");
    return 0;
}
