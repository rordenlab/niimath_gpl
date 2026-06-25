/* cost.h — SPM optfun port: smooth joint histogram + cost (mi/nmi/ecc/ncc). */
#ifndef SPMCOREG_COST_H
#define SPMCOREG_COST_H
#include "loaduint8.h"
#include "hist2.h"

typedef enum { COST_MI, COST_NMI, COST_ECC, COST_NCC } cost_fun_t;

/* Objective o(P) minimised by Powell, == SPM optfun. fwhm is the 2-vector
 * histogram smoothing (default [7 7]). */
double coreg_cost(const Vol *VG, const Vol *VF, const double *P, int np,
                  double samp, cost_fun_t cf, const double fwhm[2]);

/* Cost from an already-built 256x256 histogram (allocates; used by tests). */
double coreg_cost_from_H(const double H256[65536], cost_fun_t cf, const double fwhm[2]);

/* Reusable per-pass scratch: histogram-smoothing kernels + buffers, computed/sized
 * once from fwhm and reused every eval (the hot path allocates and recomputes nothing). */
typedef struct {
    double *tmp, *out, *s1, *s2;   /* smoothing + marginal buffers */
    int Nr, Nc;                    /* smoothed-histogram dims */
    double krn1[64], krn2[64];     /* normalized separable smoothing kernels */
    int L1, L2;                    /* kernel lengths */
} CostScratch;
int  cost_scratch_init(CostScratch *cs, const double fwhm[2]);
void cost_scratch_free(CostScratch *cs);

/* Same objective using precomputed base samples + scratch (the optimizer hot path). */
double coreg_cost_cached(const BaseSamples *base, const Vol *VF, const double *P, int np,
                         cost_fun_t cf, const double fwhm[2], CostScratch *cs);

int cost_parse(const char *s, cost_fun_t *out); /* "mi"/"nmi"/"ecc"/"ncc" -> enum */

#endif
