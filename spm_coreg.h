/* spm_coreg.h — top-level estimate driver (spm_coreg.m). */
#ifndef SPMCOREG_SPM_COREG_H
#define SPMCOREG_SPM_COREG_H
#include "cost.h"

typedef struct {
    cost_fun_t cf;        /* default NMI */
    double sep[4]; int nsep;   /* default [4 2] */
    double fwhm[2];       /* histogram smoothing, default [7 7] */
    double tol[6];        /* per-param accuracy, default [.02 .02 .02 .001 .001 .001] */
    double params0[6];    /* starting estimate, default 0 */
    int dither;           /* loaduint8 dither (0 = validation/reproducible) */
    int coarse_downsample;/* 1 = block-average coarse passes (fast); 0 = SPM sparse sampling (faithful) */
} coreg_flags;

void coreg_default_flags(coreg_flags *f);

/* Validate flags (sep/fwhm/tol/dither). Returns 0 if usable, 1 otherwise. */
int coreg_validate_flags(const coreg_flags *f);

/* Estimate rigid-body params registering src(moving) -> ref(stationary).
 * Returns 0 on success; xk holds the 6 params (VF.mat\spm_matrix(xk)*VG.mat). */
int coreg_estimate(const char *ref, const char *src, const coreg_flags *f, double xk[6]);

/* In-memory variant (niimath chain): ref/src are already-loaded images, not paths.
 * Does not take ownership of either nim. */
int coreg_estimate_nim(const nifti_image *ref, const nifti_image *src, const coreg_flags *f, double xk[6]);

#endif
