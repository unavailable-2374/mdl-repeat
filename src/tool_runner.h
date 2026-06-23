#ifndef MDL_TOOL_RUNNER_H
#define MDL_TOOL_RUNNER_H

#include <stddef.h>

typedef enum {
    TOOL_OK = 0,
    TOOL_NOT_FOUND,
    TOOL_EXIT_NONZERO,
    TOOL_TIMEOUT,
    TOOL_IO_ERROR
} ToolStatus;

typedef struct {
    const char *name;
    const char *path_override;
    char *const *argv;
    const char *stdin_path;
    const char *stdout_path;
    const char *stderr_path;
    int timeout_sec;
    int required;
} ToolCommand;

typedef struct {
    ToolStatus status;
    int exit_code;
    double wall_sec;
    char resolved_path[4096];
    char error_message[512];
} ToolResult;

const char *tool_status_name(ToolStatus status);
int tool_find(const char *name, const char *path_override,
              char *resolved, size_t resolved_size);
ToolResult tool_run(const ToolCommand *cmd);

#endif
