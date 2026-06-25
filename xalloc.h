/* xalloc.h — allocate-or-die for leaf/scratch paths in a CLI tool.
 * Entry points (load_uint8, reslice_trilinear) still return error codes;
 * these cover deep scratch allocations where OOM is unrecoverable anyway. */
#ifndef SPMCOREG_XALLOC_H
#define SPMCOREG_XALLOC_H
#include <stdlib.h>
#include <stdio.h>

static inline void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) { fprintf(stderr, "spmcoreg: out of memory (%zu bytes)\n", n); exit(1); }
    return p;
}
static inline void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n, s);
    if (!p) { fprintf(stderr, "spmcoreg: out of memory (%zu x %zu bytes)\n", n, s); exit(1); }
    return p;
}
#endif
