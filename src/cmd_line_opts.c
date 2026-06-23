#include "cmd_line_opts.h"

#include <stdlib.h>
#include <string.h>

static int co_find_option(int argc, char** argv, const char* text) {
	int ii;
	for(ii = 1; ii < argc; ++ii) {
		if(strcmp(argv[ii], text) == 0) {
			return ii;
		}
	}

	return -1;
}

int co_get_int(int argc, char** argv, const char* text, int* res) {
	int ii = co_find_option(argc, argv, text);
	if(ii < 0 || ii + 1 >= argc) {
		return 0;
	}

	*res = atoi(argv[ii+1]);
	return 1;
}

int co_get_bool(int argc, char** argv, const char* text, int* res) {
	int ii = co_find_option(argc, argv, text);
	if(ii < 0) {
		return 0;
	}

	if(res) {
		*res = 1;
	}
	return 1;
}

int co_get_float(int argc, char** argv, const char* text, float* res) {
	int ii = co_find_option(argc, argv, text);
	if(ii < 0 || ii + 1 >= argc) {
		return 0;
	}

	*res = atof(argv[ii+1]);
	return 1;
}

int co_get_string(int argc, char** argv, const char* text, char* * res) {
	int ii = co_find_option(argc, argv, text);
	if(ii < 0 || ii + 1 >= argc || argv[ii+1][0] == '-') {
		return 0;
	}

	*res = argv[ii+1];
	return 1;
}

int co_has_option(int argc, char** argv, const char* text) {
	return co_find_option(argc, argv, text) >= 0;
}
