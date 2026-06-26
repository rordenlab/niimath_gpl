/* xalloc.h — allocate-or-bail for leaf/scratch paths.
 * Entry points (load_uint8, reslice_trilinear) still return error codes; these
 * cover deep scratch allocations where threading a return code through every
 * caller (cost/smooth/powell) would be churn. On OOM we longjmp back to the
 * guard armed in coreg_run() so the caller gets an error return instead of a
 * process exit — required for the long-lived WASM Web Worker, where exit(1)
 * kills the whole niimath instance. If no guard is armed (shouldn't happen via
 * the library entry points) we fall back to exit(1). */
#ifndef SPMCOREG_XALLOC_H
#define SPMCOREG_XALLOC_H
#include <stdlib.h>
#include <stdio.h>
#include <setjmp.h>

/* Thread-local when OpenMP is in play (the project runs one subject per thread),
 * matching powell_newuoa.c's PNL_TLOCAL / allineate.c's AL_TLOCAL. Plain globals
 * otherwise (also avoids __thread on Windows/non-OMP builds, as those files do). */
#ifdef _OPENMP
#define SC_TLOCAL __thread
#else
#define SC_TLOCAL
#endif

extern SC_TLOCAL jmp_buf spmcoreg_oom_jmp;   /* armed per-thread by coreg_run() */
extern SC_TLOCAL int spmcoreg_oom_armed;

static inline void spmcoreg_oom(void) {
    if (spmcoreg_oom_armed) { spmcoreg_oom_armed = 0; longjmp(spmcoreg_oom_jmp, 1); }
    exit(1);   /* no guard installed: preserve the original allocate-or-die CLI behavior */
}
static inline void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) { fprintf(stderr, "spmcoreg: out of memory (%zu bytes)\n", n); spmcoreg_oom(); }
    return p;
}
static inline void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n, s);
    if (!p) { fprintf(stderr, "spmcoreg: out of memory (%zu x %zu bytes)\n", n, s); spmcoreg_oom(); }
    return p;
}
#endif
