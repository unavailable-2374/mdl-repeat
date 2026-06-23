#include "external_qc.h"

#include <stdio.h>
#include <string.h>

int external_tools_mode_parse(const char *text, ExternalToolsMode *out)
{
    if (!text || !out) return 0;
    if (strcmp(text, "off") == 0) {
        *out = EXTERNAL_TOOLS_OFF;
        return 1;
    }
    if (strcmp(text, "auto") == 0) {
        *out = EXTERNAL_TOOLS_AUTO;
        return 1;
    }
    if (strcmp(text, "require") == 0) {
        *out = EXTERNAL_TOOLS_REQUIRE;
        return 1;
    }
    return 0;
}

const char *external_tools_mode_name(ExternalToolsMode mode)
{
    switch (mode) {
        case EXTERNAL_TOOLS_OFF: return "off";
        case EXTERNAL_TOOLS_AUTO: return "auto";
        case EXTERNAL_TOOLS_REQUIRE: return "require";
        default: return "unknown";
    }
}

int external_qc_run_seqkit_stats(const ExternalQcConfig *cfg,
                                 const char *library_fasta)
{
    if (!cfg || !cfg->qc_output_path || !library_fasta)
        return 0;
    if (cfg->mode == EXTERNAL_TOOLS_OFF)
        return 0;

    char *argv[] = {
        "seqkit",
        "stats",
        "-T",
        (char *)library_fasta,
        NULL
    };

    ToolCommand cmd;
    memset(&cmd, 0, sizeof(cmd));
    cmd.name = "seqkit";
    cmd.path_override = cfg->seqkit_path;
    cmd.argv = argv;
    cmd.stdout_path = cfg->qc_output_path;
    cmd.timeout_sec = cfg->timeout_sec > 0 ? cfg->timeout_sec : 300;
    cmd.required = (cfg->mode == EXTERNAL_TOOLS_REQUIRE);

    ToolResult result = tool_run(&cmd);
    if (result.status == TOOL_OK) {
        fprintf(stderr, "External QC: seqkit stats wrote %s\n",
                cfg->qc_output_path);
        return 0;
    }

    if (cfg->mode == EXTERNAL_TOOLS_REQUIRE) {
        fprintf(stderr, "ERROR: required external QC seqkit failed "
                        "(status=%s, exit=%d): %s\n",
                tool_status_name(result.status), result.exit_code,
                result.error_message);
        return 1;
    }

    fprintf(stderr, "WARNING: optional external QC seqkit skipped/failed "
                    "(status=%s, exit=%d): %s\n",
            tool_status_name(result.status), result.exit_code,
            result.error_message);
    return 0;
}
