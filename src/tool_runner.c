#define _POSIX_C_SOURCE 200809L

#include "tool_runner.h"

#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

static int has_slash(const char *s)
{
    return s && strchr(s, '/') != NULL;
}

static int is_executable(const char *path)
{
    return path && path[0] && access(path, X_OK) == 0;
}

static void set_error(ToolResult *r, const char *prefix, const char *detail)
{
    if (!r) return;
    snprintf(r->error_message, sizeof(r->error_message), "%s%s%s",
             prefix ? prefix : "", detail ? ": " : "", detail ? detail : "");
}

static double elapsed_sec(struct timespec start, struct timespec end)
{
    return (double)(end.tv_sec - start.tv_sec) +
           (double)(end.tv_nsec - start.tv_nsec) / 1e9;
}

static void sleep_poll_interval(void)
{
    struct timespec ts;
    ts.tv_sec = 0;
    ts.tv_nsec = 10 * 1000 * 1000;
    nanosleep(&ts, NULL);
}

const char *tool_status_name(ToolStatus status)
{
    switch (status) {
        case TOOL_OK: return "ok";
        case TOOL_NOT_FOUND: return "not_found";
        case TOOL_EXIT_NONZERO: return "exit_nonzero";
        case TOOL_TIMEOUT: return "timeout";
        case TOOL_IO_ERROR: return "io_error";
        default: return "unknown";
    }
}

int tool_find(const char *name, const char *path_override,
              char *resolved, size_t resolved_size)
{
    if (!resolved || resolved_size == 0) return 0;
    resolved[0] = '\0';

    const char *candidate = path_override && path_override[0]
                          ? path_override : name;
    if (!candidate || !candidate[0]) return 0;

    if (has_slash(candidate)) {
        if (!is_executable(candidate)) return 0;
        snprintf(resolved, resolved_size, "%s", candidate);
        return 1;
    }

    const char *path = getenv("PATH");
    if (!path || !path[0]) path = ".";

    char *path_copy = malloc(strlen(path) + 1);
    if (!path_copy) return 0;
    strcpy(path_copy, path);

    int found = 0;
    for (char *dir = strtok(path_copy, ":");
         dir;
         dir = strtok(NULL, ":")) {
        if (!dir[0]) dir = ".";
        char buf[4096];
        int n = snprintf(buf, sizeof(buf), "%s/%s", dir, candidate);
        if (n < 0 || (size_t)n >= sizeof(buf)) continue;
        if (is_executable(buf)) {
            snprintf(resolved, resolved_size, "%s", buf);
            found = 1;
            break;
        }
    }

    free(path_copy);
    return found;
}

static int redirect_fd(const char *path, int target_fd, int flags, mode_t mode)
{
    if (!path) return 0;
    int fd = open(path, flags, mode);
    if (fd < 0) return -1;
    if (dup2(fd, target_fd) < 0) {
        close(fd);
        return -1;
    }
    close(fd);
    return 0;
}

ToolResult tool_run(const ToolCommand *cmd)
{
    ToolResult r;
    memset(&r, 0, sizeof(r));
    r.status = TOOL_IO_ERROR;
    r.exit_code = -1;

    if (!cmd || !cmd->name || !cmd->argv || !cmd->argv[0]) {
        set_error(&r, "invalid tool command", NULL);
        return r;
    }

    if (!tool_find(cmd->name, cmd->path_override,
                   r.resolved_path, sizeof(r.resolved_path))) {
        r.status = TOOL_NOT_FOUND;
        set_error(&r, "tool not found", cmd->path_override ?
                  cmd->path_override : cmd->name);
        return r;
    }

    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC, &start);

    pid_t pid = fork();
    if (pid < 0) {
        set_error(&r, "fork failed", strerror(errno));
        return r;
    }

    if (pid == 0) {
        if (redirect_fd(cmd->stdin_path, STDIN_FILENO, O_RDONLY, 0) != 0)
            _exit(126);
        if (redirect_fd(cmd->stdout_path, STDOUT_FILENO,
                        O_WRONLY | O_CREAT | O_TRUNC, 0644) != 0)
            _exit(126);
        if (redirect_fd(cmd->stderr_path, STDERR_FILENO,
                        O_WRONLY | O_CREAT | O_TRUNC, 0644) != 0)
            _exit(126);

        execv(r.resolved_path, cmd->argv);
        _exit(errno == ENOENT ? 127 : 126);
    }

    int status = 0;
    while (1) {
        pid_t w = waitpid(pid, &status, WNOHANG);
        if (w == pid) break;
        if (w < 0) {
            r.status = TOOL_IO_ERROR;
            set_error(&r, "waitpid failed", strerror(errno));
            return r;
        }

        clock_gettime(CLOCK_MONOTONIC, &now);
        if (cmd->timeout_sec > 0 &&
            elapsed_sec(start, now) > (double)cmd->timeout_sec) {
            kill(pid, SIGKILL);
            waitpid(pid, &status, 0);
            clock_gettime(CLOCK_MONOTONIC, &now);
            r.wall_sec = elapsed_sec(start, now);
            r.status = TOOL_TIMEOUT;
            set_error(&r, "tool timed out", cmd->name);
            return r;
        }
        sleep_poll_interval();
    }

    clock_gettime(CLOCK_MONOTONIC, &now);
    r.wall_sec = elapsed_sec(start, now);

    if (WIFEXITED(status)) {
        r.exit_code = WEXITSTATUS(status);
        if (r.exit_code == 0) {
            r.status = TOOL_OK;
        } else if (r.exit_code == 126 || r.exit_code == 127) {
            r.status = TOOL_IO_ERROR;
            set_error(&r, "tool execution failed", cmd->name);
        } else {
            r.status = TOOL_EXIT_NONZERO;
            snprintf(r.error_message, sizeof(r.error_message),
                     "tool exited with code %d", r.exit_code);
        }
    } else if (WIFSIGNALED(status)) {
        r.exit_code = 128 + WTERMSIG(status);
        r.status = TOOL_EXIT_NONZERO;
        snprintf(r.error_message, sizeof(r.error_message),
                 "tool terminated by signal %d", WTERMSIG(status));
    } else {
        r.status = TOOL_IO_ERROR;
        set_error(&r, "unknown tool status", cmd->name);
    }

    (void)cmd->required;
    return r;
}
