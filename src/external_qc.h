#ifndef MDL_EXTERNAL_QC_H
#define MDL_EXTERNAL_QC_H

#include "tool_runner.h"

typedef enum {
    EXTERNAL_TOOLS_OFF = 0,
    EXTERNAL_TOOLS_AUTO,
    EXTERNAL_TOOLS_REQUIRE
} ExternalToolsMode;

typedef struct {
    ExternalToolsMode mode;
    const char *seqkit_path;
    const char *qc_output_path;
    int timeout_sec;
} ExternalQcConfig;

int external_tools_mode_parse(const char *text, ExternalToolsMode *out);
const char *external_tools_mode_name(ExternalToolsMode mode);

/*
 * Run optional non-mutating FASTA QC on the final library.
 * Returns 0 when QC succeeds or is intentionally skipped.
 * Returns 1 only when mode requires the tool/QC result and that requirement fails.
 */
int external_qc_run_seqkit_stats(const ExternalQcConfig *cfg,
                                 const char *library_fasta);

#endif
