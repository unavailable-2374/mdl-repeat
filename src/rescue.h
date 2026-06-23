#ifndef MDL_RESCUE_H
#define MDL_RESCUE_H

#include "types.h"
#include "genome.h"
#include "candidates.h"
#include "discover.h"

typedef struct {
    int    seq_index;
    gpos_t raw_start;
    gpos_t raw_end;
    glen_t seg_length;
} SeqSegment;

typedef struct {
    int full_genome;
    int l_delta;
    int max_repeats;
    int min_gap;
    int verbose;
    int num_threads;
    glen_t chunk_size;
    const char *audit_file;
} RecallRescueOptions;

typedef Genome *(*RescueCreateChunkFn)(const Genome *genome,
                                       const SeqSegment *segments,
                                       int num_segments);
typedef void (*RescueFreeChunkFn)(Genome *chunk);
typedef CandidateList *(*RescueDiscoverChunkedFn)(const Genome *genome,
                                                  const DiscoverParams *params,
                                                  int verbose,
                                                  int num_threads,
                                                  glen_t chunk_size);
typedef void (*RescueRemapInstanceCoordinatesFn)(CandidateList *candidates,
                                                 const SeqSegment *segments,
                                                 int num_segments);

typedef struct {
    RescueCreateChunkFn create_chunk;
    RescueFreeChunkFn free_chunk;
    RescueDiscoverChunkedFn discover_chunked;
    RescueRemapInstanceCoordinatesFn remap_instance_coordinates;
} RecallRescueCallbacks;

int recall_rescue_run(const Genome *genome, CandidateList *candidates,
                      const DiscoverParams *base_params,
                      const RecallRescueOptions *options,
                      const RecallRescueCallbacks *callbacks);

#endif /* MDL_RESCUE_H */
