#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rescue.h"

typedef struct {
    int matched_family_id;
    double consensus_identity;
    double contained_instance_fraction;
    double length_ratio;
} RescueDuplicateEvidence;

typedef struct {
    int candidates;
    int near_duplicates;
    double max_consensus_identity;
    double max_contained_instance_fraction;
    double max_length_ratio;
} RescueAppendStats;

typedef struct {
    int seq_index;
    gpos_t start;
    gpos_t end;
} RawInterval;

static double candidate_instance_contained_fraction(const CandidateFamily *query,
                                                    const CandidateFamily *target)
{
    if (!query || !target || query->num_instances <= 0 ||
        target->num_instances <= 0)
        return 0.0;

    int contained = 0;
    for (int i = 0; i < query->num_instances; i++) {
        const Instance *qi = &query->instances[i];
        if (qi->aligned_length <= 0) continue;
        gpos_t qs = qi->position;
        gpos_t qe = qs + (gpos_t)qi->aligned_length;
        int matched = 0;

        for (int j = 0; j < target->num_instances; j++) {
            const Instance *tj = &target->instances[j];
            if (qi->seq_index != tj->seq_index) continue;
            if (tj->aligned_length <= 0) continue;

            gpos_t ts = tj->position;
            gpos_t te = ts + (gpos_t)tj->aligned_length;
            gpos_t os = qs > ts ? qs : ts;
            gpos_t oe = qe < te ? qe : te;
            if (oe <= os) continue;
            if ((oe - os) * 100 >= (gpos_t)qi->aligned_length * 80) {
                matched = 1;
                break;
            }
        }

        if (matched) contained++;
    }

    return (double)contained / (double)query->num_instances;
}

static double consensus_window_identity_score(const char *shorter, int short_len,
                                              const char *longer, int long_len,
                                              int reverse_complement)
{
    if (!shorter || !longer || short_len <= 0 || long_len <= 0 ||
        short_len > long_len)
        return 0.0;

    int best_matches = 0;
    for (int off = 0; off <= long_len - short_len; off++) {
        int matches = 0;
        int comparable = 0;
        for (int i = 0; i < short_len; i++) {
            char a = shorter[i];
            char b = reverse_complement
                   ? (char)dna_complement(longer[off + short_len - 1 - i])
                   : longer[off + i];
            if (a == DNA_N || b == DNA_N) continue;
            comparable++;
            if (a == b) matches++;
        }
        if (comparable > 0 && matches > best_matches)
            best_matches = matches;
    }

    return (double)best_matches / (double)short_len;
}

static double candidate_consensus_identity(const CandidateFamily *a,
                                           const CandidateFamily *b,
                                           double *length_ratio_out)
{
    if (!a || !b || !a->consensus || !b->consensus ||
        a->consensus_length <= 0 || b->consensus_length <= 0)
        return 0.0;

    int min_len = a->consensus_length < b->consensus_length
                ? a->consensus_length : b->consensus_length;
    int max_len = a->consensus_length > b->consensus_length
                ? a->consensus_length : b->consensus_length;
    double length_ratio = (double)min_len / (double)max_len;
    if (length_ratio_out) *length_ratio_out = length_ratio;
    if (length_ratio < 0.80)
        return 0.0;

    const char *shorter = a->consensus;
    const char *longer = b->consensus;
    int short_len = a->consensus_length;
    int long_len = b->consensus_length;
    if (short_len > long_len) {
        shorter = b->consensus;
        longer = a->consensus;
        short_len = b->consensus_length;
        long_len = a->consensus_length;
    }

    double fwd = consensus_window_identity_score(shorter, short_len,
                                                longer, long_len, 0);
    double rev = consensus_window_identity_score(shorter, short_len,
                                                longer, long_len, 1);
    return fwd > rev ? fwd : rev;
}

static int rescue_candidate_is_duplicate(const CandidateList *dst,
                                         const CandidateFamily *query,
                                         RescueDuplicateEvidence *evidence)
{
    if (!dst || !query) return 0;
    RescueDuplicateEvidence best;
    best.matched_family_id = -1;
    best.consensus_identity = 0.0;
    best.contained_instance_fraction = 0.0;
    best.length_ratio = 0.0;
    int duplicate = 0;

    for (int i = 0; i < dst->num_families; i++) {
        double length_ratio = 0.0;
        double identity = candidate_consensus_identity(query, &dst->families[i],
                                                       &length_ratio);
        double contained = candidate_instance_contained_fraction(query,
                                                                 &dst->families[i]);
        double score = identity * contained * length_ratio;
        double best_score = best.consensus_identity *
                            best.contained_instance_fraction *
                            best.length_ratio;
        if (score > best_score) {
            best.matched_family_id = dst->families[i].id;
            best.consensus_identity = identity;
            best.contained_instance_fraction = contained;
            best.length_ratio = length_ratio;
        }
        if (length_ratio >= 0.80 && identity >= 0.80 && contained >= 0.80) {
            best.matched_family_id = dst->families[i].id;
            best.consensus_identity = identity;
            best.contained_instance_fraction = contained;
            best.length_ratio = length_ratio;
            duplicate = 1;
            break;
        }
    }

    if (evidence) *evidence = best;
    return duplicate;
}

static FILE *rescue_audit_open(const char *path)
{
    if (!path) return NULL;
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "WARNING: could not open rescue audit '%s'\n", path);
        return NULL;
    }
    fprintf(fp, "record_type\tdecision\treason\tsource_family_id\t"
                "output_family_id\tmatched_family_id\tseq_index\traw_start\t"
                "raw_end\tlength\tconsensus_length\tnum_instances\t"
                "consensus_identity\tcontained_instance_fraction\tlength_ratio\n");
    return fp;
}

static void rescue_audit_write_segment(FILE *fp, int source_id,
                                       const SeqSegment *seg)
{
    if (!fp || !seg) return;
    fprintf(fp, "target_segment\ttarget\tuncovered_gap\t%d\t-1\t-1\t"
                "%d\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t0\t0\t"
                "0.000000\t0.000000\t0.000000\n",
            source_id, seg->seq_index, (int64_t)seg->raw_start,
            (int64_t)seg->raw_end, (int64_t)seg->seg_length);
}

static void rescue_audit_write_candidate(FILE *fp, const char *decision,
                                         const char *reason,
                                         int source_family_id,
                                         int output_family_id,
                                         int matched_family_id,
                                         const CandidateFamily *fam,
                                         const RescueDuplicateEvidence *evidence)
{
    if (!fp || !fam) return;
    double identity = evidence ? evidence->consensus_identity : 0.0;
    double contained = evidence ? evidence->contained_instance_fraction : 0.0;
    double length_ratio = evidence ? evidence->length_ratio : 0.0;
    fprintf(fp, "candidate\t%s\t%s\t%d\t%d\t%d\t-1\t0\t0\t0\t%d\t%d\t"
                "%.6f\t%.6f\t%.6f\n",
            decision ? decision : "unknown",
            reason ? reason : "unknown",
            source_family_id, output_family_id, matched_family_id,
            fam->consensus_length, fam->num_instances,
            identity, contained, length_ratio);
}

static int append_candidate_list(CandidateList *dst, CandidateList *src,
                                 uint32_t discovery_flags,
                                 int *out_duplicates,
                                 FILE *rescue_audit_fp,
                                 RescueAppendStats *rescue_stats)
{
    if (!dst || !src) return 0;
    if (out_duplicates) *out_duplicates = 0;
    if (rescue_stats) memset(rescue_stats, 0, sizeof(*rescue_stats));
    if (src->num_families <= 0) {
        candidates_free(src);
        return 0;
    }

    int needed = dst->num_families + src->num_families;
    if (needed > dst->cap_families) {
        int new_cap = dst->cap_families > 0 ? dst->cap_families : 1;
        while (new_cap < needed) new_cap *= 2;
        CandidateFamily *tmp = realloc(dst->families,
                                       (size_t)new_cap * sizeof(CandidateFamily));
        if (!tmp) return -1;
        dst->families = tmp;
        dst->cap_families = new_cap;
    }

    int appended = 0;
    for (int i = 0; i < src->num_families; i++) {
        RescueDuplicateEvidence evidence;
        evidence.matched_family_id = -1;
        evidence.consensus_identity = 0.0;
        evidence.contained_instance_fraction = 0.0;
        evidence.length_ratio = 0.0;
        int is_rescue = (discovery_flags & CAND_QF_RESCUE_DISCOVERY) != 0;
        int is_duplicate = is_rescue &&
            rescue_candidate_is_duplicate(dst, &src->families[i], &evidence);
        if (is_rescue && rescue_stats) {
            rescue_stats->candidates++;
            if (evidence.consensus_identity >
                rescue_stats->max_consensus_identity)
                rescue_stats->max_consensus_identity =
                    evidence.consensus_identity;
            if (evidence.contained_instance_fraction >
                rescue_stats->max_contained_instance_fraction)
                rescue_stats->max_contained_instance_fraction =
                    evidence.contained_instance_fraction;
            if (evidence.length_ratio > rescue_stats->max_length_ratio)
                rescue_stats->max_length_ratio = evidence.length_ratio;
            if (!is_duplicate &&
                evidence.length_ratio >= 0.80 &&
                evidence.consensus_identity >= 0.80 &&
                evidence.contained_instance_fraction >= 0.50)
                rescue_stats->near_duplicates++;
        }
        if (is_duplicate) {
            rescue_audit_write_candidate(rescue_audit_fp, "filter",
                                         "duplicate_existing",
                                         src->families[i].id, -1,
                                         evidence.matched_family_id,
                                         &src->families[i], &evidence);
            free(src->families[i].consensus);
            free(src->families[i].instances);
            memset(&src->families[i], 0, sizeof(src->families[i]));
            if (out_duplicates) (*out_duplicates)++;
            continue;
        }

        CandidateFamily *out = &dst->families[dst->num_families];
        *out = src->families[i];
        out->id = (mdl_uid_t)dst->num_families;
        out->discovery_flags |= discovery_flags;
        rescue_audit_write_candidate(rescue_audit_fp, "append", "kept",
                                     src->families[i].id, out->id,
                                     evidence.matched_family_id, out,
                                     &evidence);
        memset(&src->families[i], 0, sizeof(src->families[i]));
        dst->num_families++;
        appended++;
    }

    free(src->families);
    src->families = NULL;
    src->num_families = 0;
    candidates_free(src);
    return appended;
}

static int cmp_raw_interval(const void *a, const void *b)
{
    const RawInterval *ia = (const RawInterval *)a;
    const RawInterval *ib = (const RawInterval *)b;
    if (ia->seq_index != ib->seq_index)
        return ia->seq_index - ib->seq_index;
    if (ia->start < ib->start) return -1;
    if (ia->start > ib->start) return  1;
    if (ia->end < ib->end) return -1;
    if (ia->end > ib->end) return  1;
    return 0;
}

static int cmp_seg_coord_asc(const void *a, const void *b)
{
    const SeqSegment *sa = (const SeqSegment *)a;
    const SeqSegment *sb = (const SeqSegment *)b;
    if (sa->seq_index != sb->seq_index)
        return sa->seq_index - sb->seq_index;
    if (sa->raw_start < sb->raw_start) return -1;
    if (sa->raw_start > sb->raw_start) return  1;
    if (sa->raw_end < sb->raw_end) return -1;
    if (sa->raw_end > sb->raw_end) return  1;
    return 0;
}

static int build_rescue_segments_from_gaps(const Genome *genome,
                                           const CandidateList *cl,
                                           glen_t min_segment_len,
                                           glen_t flank,
                                           SeqSegment **out_segments,
                                           int *out_num_segments)
{
    *out_segments = NULL;
    *out_num_segments = 0;

    int64_t n_inst = 0;
    for (int f = 0; f < cl->num_families; f++)
        n_inst += cl->families[f].num_instances;

    RawInterval *iv = NULL;
    if (n_inst > 0) {
        iv = malloc((size_t)n_inst * sizeof(RawInterval));
        if (!iv) return -1;
    }

    int n_iv = 0;
    for (int f = 0; f < cl->num_families; f++) {
        const CandidateFamily *fam = &cl->families[f];
        for (int i = 0; i < fam->num_instances; i++) {
            const Instance *inst = &fam->instances[i];
            int seq = inst->seq_index;
            if (seq < 0 || seq >= genome->num_sequences)
                seq = genome_get_seq_index(genome, inst->position);
            if (seq < 0 || seq >= genome->num_sequences) continue;

            gpos_t seq_start = (seq == 0) ? 0 : genome->boundaries[seq - 1];
            gpos_t seq_end = genome->boundaries[seq];
            if (seq == genome->num_sequences - 1) seq_end -= 1;

            gpos_t s = inst->position - PADLENGTH;
            gpos_t e = s + (gpos_t)inst->aligned_length;
            if (s < seq_start) s = seq_start;
            if (e > seq_end) e = seq_end;
            if (e <= s) continue;

            iv[n_iv].seq_index = seq;
            iv[n_iv].start = s;
            iv[n_iv].end = e;
            n_iv++;
        }
    }

    if (n_iv > 0)
        qsort(iv, (size_t)n_iv, sizeof(RawInterval), cmp_raw_interval);

    int seg_cap = 128;
    int n_seg = 0;
    SeqSegment *segments = malloc((size_t)seg_cap * sizeof(SeqSegment));
    if (!segments) { free(iv); return -1; }

    int iv_idx = 0;
    for (int seq = 0; seq < genome->num_sequences; seq++) {
        gpos_t seq_start = (seq == 0) ? 0 : genome->boundaries[seq - 1];
        gpos_t seq_end = genome->boundaries[seq];
        if (seq == genome->num_sequences - 1) seq_end -= 1;
        gpos_t cursor = seq_start;

        while (iv_idx < n_iv && iv[iv_idx].seq_index < seq)
            iv_idx++;

        while (iv_idx < n_iv && iv[iv_idx].seq_index == seq) {
            gpos_t cov_s = iv[iv_idx].start;
            gpos_t cov_e = iv[iv_idx].end;
            iv_idx++;
            while (iv_idx < n_iv && iv[iv_idx].seq_index == seq &&
                   iv[iv_idx].start <= cov_e) {
                if (iv[iv_idx].end > cov_e) cov_e = iv[iv_idx].end;
                iv_idx++;
            }

            if (cov_s > cursor && cov_s - cursor >= min_segment_len) {
                gpos_t s = cursor - flank;
                gpos_t e = cov_s + flank;
                if (s < seq_start) s = seq_start;
                if (e > seq_end) e = seq_end;
                if (e > s) {
                    if (n_seg >= seg_cap) {
                        seg_cap *= 2;
                        SeqSegment *tmp = realloc(segments,
                            (size_t)seg_cap * sizeof(SeqSegment));
                        if (!tmp) { free(segments); free(iv); return -1; }
                        segments = tmp;
                    }
                    segments[n_seg].seq_index = seq;
                    segments[n_seg].raw_start = s;
                    segments[n_seg].raw_end = e;
                    segments[n_seg].seg_length = e - s;
                    n_seg++;
                }
            }
            if (cov_e > cursor) cursor = cov_e;
        }

        if (seq_end > cursor && seq_end - cursor >= min_segment_len) {
            gpos_t s = cursor - flank;
            gpos_t e = seq_end;
            if (s < seq_start) s = seq_start;
            if (e > s) {
                if (n_seg >= seg_cap) {
                    seg_cap *= 2;
                    SeqSegment *tmp = realloc(segments,
                        (size_t)seg_cap * sizeof(SeqSegment));
                    if (!tmp) { free(segments); free(iv); return -1; }
                    segments = tmp;
                }
                segments[n_seg].seq_index = seq;
                segments[n_seg].raw_start = s;
                segments[n_seg].raw_end = e;
                segments[n_seg].seg_length = e - s;
                n_seg++;
            }
        }
    }

    free(iv);
    if (n_seg == 0) {
        free(segments);
        return 0;
    }

    qsort(segments, (size_t)n_seg, sizeof(SeqSegment), cmp_seg_coord_asc);
    int write_idx = 0;
    for (int i = 0; i < n_seg; i++) {
        if (write_idx > 0 &&
            segments[i].seq_index == segments[write_idx - 1].seq_index &&
            segments[i].raw_start <= segments[write_idx - 1].raw_end) {
            if (segments[i].raw_end > segments[write_idx - 1].raw_end) {
                segments[write_idx - 1].raw_end = segments[i].raw_end;
                segments[write_idx - 1].seg_length =
                    segments[write_idx - 1].raw_end -
                    segments[write_idx - 1].raw_start;
            }
        } else {
            if (write_idx != i)
                segments[write_idx] = segments[i];
            write_idx++;
        }
    }
    n_seg = write_idx;

    *out_segments = segments;
    *out_num_segments = n_seg;
    return n_seg;
}

int recall_rescue_run(const Genome *genome, CandidateList *candidates,
                      const DiscoverParams *base_params,
                      const RecallRescueOptions *options,
                      const RecallRescueCallbacks *callbacks)
{
    if (!genome || !candidates || !base_params || !options || !callbacks)
        return -1;

    DiscoverParams rescue_params = *base_params;
    int base_l = base_params->l;
    if (base_l <= 0)
        base_l = (int)ceil(1.0 + log((double)genome->length) / log(4.0));
    rescue_params.l = base_l - options->l_delta;
    if (rescue_params.l < 8) rescue_params.l = 8;
    rescue_params.MAXR = options->max_repeats;
    rescue_params.freq_file = NULL;
    rescue_params.freq_output = NULL;

    fprintf(stderr, "Recall rescue discovery: l=%d (base=%d, delta=%d), "
            "MAXR=%d, mode=%s\n",
            rescue_params.l, base_l, options->l_delta, rescue_params.MAXR,
            options->full_genome ? "full-genome" : "targeted");
    FILE *rescue_audit_fp = rescue_audit_open(options->audit_file);

    CandidateList *rescued = NULL;
    SeqSegment *rescue_segments = NULL;
    int n_rescue_segments = 0;
    Genome *rescue_genome = NULL;

    if (!options->full_genome) {
        if (!callbacks->create_chunk || !callbacks->free_chunk ||
            !callbacks->remap_instance_coordinates) {
            fprintf(stderr, "ERROR: Recall rescue targeted mode missing callbacks\n");
            if (rescue_audit_fp) fclose(rescue_audit_fp);
            return -1;
        }

        int n_targets = build_rescue_segments_from_gaps(
            genome, candidates, (glen_t)options->min_gap,
            (glen_t)rescue_params.L, &rescue_segments,
            &n_rescue_segments);
        if (n_targets < 0) {
            fprintf(stderr, "ERROR: could not build targeted rescue segments\n");
            if (rescue_audit_fp) fclose(rescue_audit_fp);
            return -1;
        }
        if (n_targets > 0) {
            rescue_genome = callbacks->create_chunk(genome, rescue_segments,
                                                    n_rescue_segments);
            if (!rescue_genome) {
                fprintf(stderr, "ERROR: could not create targeted rescue genome\n");
                free(rescue_segments);
                if (rescue_audit_fp) fclose(rescue_audit_fp);
                return -1;
            }
            fprintf(stderr, "Recall rescue targeted %d uncovered segments "
                    "(%.1fM bases, min_gap=%d)\n",
                    n_rescue_segments,
                    (double)rescue_genome->raw_length / 1e6,
                    options->min_gap);
            for (int s = 0; s < n_rescue_segments; s++)
                rescue_audit_write_segment(rescue_audit_fp, s,
                                           &rescue_segments[s]);
            rescued = discover_families(rescue_genome, &rescue_params,
                                        options->num_threads);
            if (rescued)
                callbacks->remap_instance_coordinates(rescued, rescue_segments,
                                                      n_rescue_segments);
        } else {
            fprintf(stderr, "Recall rescue targeted 0 uncovered segments; "
                    "skipping secondary discovery\n");
            rescued = calloc(1, sizeof(CandidateList));
            if (rescued) {
                rescued->cap_families = 1;
                rescued->families = calloc(1, sizeof(CandidateFamily));
            }
        }
    } else if (genome->raw_length > options->chunk_size) {
        if (!callbacks->discover_chunked) {
            fprintf(stderr, "ERROR: Recall rescue full-genome mode missing chunk callback\n");
            if (rescue_audit_fp) fclose(rescue_audit_fp);
            return -1;
        }
        rescued = callbacks->discover_chunked(genome, &rescue_params,
                                              options->verbose,
                                              options->num_threads,
                                              options->chunk_size);
    } else {
        rescued = discover_families(genome, &rescue_params,
                                    options->num_threads);
    }

    if (!rescued) {
        fprintf(stderr, "ERROR: Recall rescue discovery failed\n");
        if (rescue_audit_fp) fclose(rescue_audit_fp);
        if (rescue_genome) callbacks->free_chunk(rescue_genome);
        free(rescue_segments);
        return -1;
    }

    int rescue_duplicates = 0;
    RescueAppendStats rescue_stats;
    int appended = append_candidate_list(candidates, rescued,
                                         CAND_QF_RESCUE_DISCOVERY,
                                         &rescue_duplicates,
                                         rescue_audit_fp,
                                         &rescue_stats);
    if (rescue_audit_fp) fclose(rescue_audit_fp);
    if (rescue_genome) callbacks->free_chunk(rescue_genome);
    free(rescue_segments);
    if (appended < 0) {
        fprintf(stderr, "ERROR: could not append recall rescue candidates\n");
        candidates_free(rescued);
        return -1;
    }
    fprintf(stderr, "Recall rescue appended %d candidate families "
            "(filtered %d duplicates), %d total\n",
            appended, rescue_duplicates, candidates->num_families);
    fprintf(stderr, "Recall rescue evidence: candidates=%d, "
            "max_identity=%.3f, max_containment=%.3f, "
            "max_length_ratio=%.3f, near_duplicates=%d\n",
            rescue_stats.candidates,
            rescue_stats.max_consensus_identity,
            rescue_stats.max_contained_instance_fraction,
            rescue_stats.max_length_ratio,
            rescue_stats.near_duplicates);

    return appended;
}
