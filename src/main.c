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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

static int default_k(glen_t len)
{
    return (int)ceil(1.0 + log((double)len) / log(4.0));
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
        "  -maxrepeats #      Max families to discover (default: 100000)\n\n"
        "Optional (refinement):\n"
        "  -threads #         Number of threads for refinement (default: 1)\n"
        "  -mdl-mode <mode>   MDL position encoding: none|exact|upper (default: exact)\n\n"
        "Optional (output):\n"
        "  -instances <file>  Output instance BED\n"
        "  -stats <file>      Output family statistics TSV\n"
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
    int num_threads = 1;
    int verbose = 0;
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
    /* If -freq matched -freq-output, clear it */
    if (freq_file && freq_output && freq_file == freq_output)
        freq_file = NULL;
    co_get_int(argc, argv, "-threads", &num_threads);
    if (num_threads < 1) num_threads = 1;

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

    fprintf(stderr, "Discovering repeat families...\n");
    CandidateList *candidates = discover_families(genome, &dparams);
    if (candidates == NULL) {
        fprintf(stderr, "ERROR: Discovery failed\n");
        genome_free(genome);
        return 1;
    }
    if (verbose) candidates_print_stats(candidates);

    /* ================================================================
     * Step 3: Build k-mer table (needed for refine pipeline)
     * ================================================================ */
    if (verbose) fprintf(stderr, "Building k-mer table for refinement (k=%d)...\n", k);
    KmerTable *kt = kmer_count(genome, k, DEFAULT_TANDEMDIST);
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
                candidates->families[write_idx].id = (uid_t)write_idx;
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

    if (candidates->num_families == 0) {
        fprintf(stderr, "No families survived compaction. Writing empty output.\n");
        output_fasta(output_file, candidates);
        candidates_free(candidates);
        kmer_free(kt);
        genome_free(genome);
        return 0;
    }

    if (verbose) candidates_print_stats(candidates);

    /* ================================================================
     * Step 5: Merge redundant families
     * ================================================================ */
    if (verbose) fprintf(stderr, "Merging redundant families...\n");
    int n_merges = refine_merge_families(candidates, genome, kt, k, verbose, num_threads);
    fprintf(stderr, "Merged %d families, %d remain\n",
            n_merges, candidates->num_families);
    if (verbose && n_merges > 0) candidates_print_stats(candidates);

    /* ================================================================
     * Step 6: Split families with bimodal divergence
     * ================================================================ */
    if (verbose) fprintf(stderr, "Splitting bimodal families...\n");
    int n_splits = refine_split_families(candidates, genome, kt, k,
                                         genome->length, verbose,
                                         candidates->num_families);
    fprintf(stderr, "Split %d families, %d total\n",
            n_splits, candidates->num_families);
    if (verbose && n_splits > 0) candidates_print_stats(candidates);

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

    /* ================================================================
     * Step 8: Prune marginal families
     * ================================================================ */
    int n_pruned = refine_prune_families(candidates, genome->length, verbose,
                                          mdl_result.num_accepted);
    if (n_pruned > 0)
        fprintf(stderr, "Pruned %d marginal families\n", n_pruned);

    /* ================================================================
     * Step 8b: Post-prune recovery pass
     * Pruning reduces R, which lowers per-instance cost ceil(log2(R)).
     * Re-score rejected families with final R — some may become viable.
     * ================================================================ */
    {
        /* Count current accepted */
        int final_R = 0;
        for (int i = 0; i < candidates->num_families; i++)
            if (candidates->families[i].mdl_score > 0 &&
                candidates->families[i].num_instances >= 2)
                final_R++;

        if (final_R < mdl_result.num_accepted) {
            /* R decreased from pruning — re-score rejected families */
            if (final_R < 2) final_R = 2;
            int n_recovered = 0;

            for (int i = 0; i < candidates->num_families; i++) {
                CandidateFamily *fam = &candidates->families[i];
                if (fam->mdl_score > 0) continue;  /* already accepted */
                if (fam->num_instances < 2) continue;
                if (fam->consensus == NULL) continue;

                /* Re-score with reduced R */
                mdl_score_family(fam, genome->length, final_R);

                if (fam->mdl_score > 0) {
                    n_recovered++;
                    if (verbose)
                        fprintf(stderr, "  Recovered F%d: mdl_score=%.1f "
                                "(R=%d)\n", fam->id, fam->mdl_score,
                                final_R);
                }
            }

            if (n_recovered > 0)
                fprintf(stderr, "Recovery pass: recovered %d families "
                        "(R: %d -> %d)\n",
                        n_recovered, mdl_result.num_accepted, final_R + n_recovered);
        }
    }

    /* ================================================================
     * Step 9: Write outputs
     * ================================================================ */
    output_fasta(output_file, candidates);

    if (instances_file)
        output_bed(instances_file, candidates, genome);

    if (stats_file)
        output_stats(stats_file, candidates);

    /* Cleanup */
    candidates_free(candidates);
    kmer_free(kt);
    genome_free(genome);
    fprintf(stderr, "Done.\n");
    return 0;
}
