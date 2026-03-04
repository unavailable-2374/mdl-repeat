#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/mdl.h"

int main(void)
{
    /* Test vectors from Rissanen's universal integer code.
     * DESIGN_DOC shows 12.344 for n=100, but that omits the 4th iterated log
     * term log2(1.450) = 0.536 > 0. Correct value: 12.880. */
    struct { int64_t n; double expected; } tests[] = {
        {1,   1.518},
        {2,   2.518},
        {3,   3.768},
        {10,  7.364},
        {100, 12.880},
    };

    int n_tests = sizeof(tests) / sizeof(tests[0]);
    int passed = 0;

    printf("L_int test vectors:\n");
    for (int i = 0; i < n_tests; i++) {
        double result = L_int(tests[i].n);
        double diff = fabs(result - tests[i].expected);
        int ok = diff < 0.01; /* within 0.01 bits */

        printf("  L_int(%" PRId64 ") = %.3f (expected %.3f) %s\n",
               tests[i].n, result, tests[i].expected,
               ok ? "PASS" : "FAIL");

        if (ok) passed++;
    }

    printf("\n%d/%d tests passed\n", passed, n_tests);
    return (passed == n_tests) ? 0 : 1;
}
