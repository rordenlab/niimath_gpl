/* hist2.h — SPM joint-histogram (spm_hist2 / hist2.c) port. */
#ifndef SPMCOREG_HIST2_H
#define SPMCOREG_HIST2_H
#include "loaduint8.h"

/* Build the 256x256 joint histogram H (row index = f value, col = g value;
 * H[ivf + ivg*256]) for parameter vector P (length np) at sampling step samp (mm).
 * M = inv(VF.mat) * spm_matrix(P) * VG.mat ;  s = samp ./ VG.vox.  H must hold 65536. */
void coreg_hist2(const Vol *VG, const Vol *VF, const double *P, int np,
                 double samp, double H[65536]);

/* full partial-volume variant (both intensity axes) — dither experiment only */
void coreg_hist2_fullpv(const Vol *VG, const Vol *VF, const double *P, int np,
                        double samp, double H[65536]);

/* --- Optimization: precomputed base (stationary) samples ---
 * The base sample positions (rx,ry,rz) and their binned values (ivg) depend only
 * on VG + samp + the dither table, NOT on the transform. Build once per pass and
 * reuse across all optimizer evals: the per-eval loop then skips the VG gather.
 * Produces bit-identical histograms to coreg_hist2. */
typedef struct {
    float *rx,*ry,*rz;      /* base sample coords (1-based, jittered) */
    unsigned char *ivg;     /* floor(samp_vol(VG)+0.5), 0..255 */
    long n;
    double mat[4][4];       /* VG SPM mat (for M = inv(VF.mat)*spm_matrix(P)*mat) */
} BaseSamples;

int  base_samples_build(const Vol *VG, double samp, BaseSamples *bs); /* 0 ok */
void base_samples_free(BaseSamples *bs);

/* Per-pass reusable thread buffer for the parallel histogram build. Allocated
 * once per pass (rides on CostScratch), reused across all evals so the hot path
 * never allocates. nt<=1 / buf==NULL (always so on non-OpenMP/WASM) => the
 * cached builder takes its serial fallback. */
typedef struct { double *buf; int nt; } Hist2Scratch;
int  hist2_scratch_init(Hist2Scratch *hs);   /* 0 ok (also ok = serial: buf=NULL,nt=1) */
void hist2_scratch_free(Hist2Scratch *hs);

/* hs may be NULL (serial). When non-NULL with hs->nt>1 and a large sample count,
 * the build is parallelized across hs->nt threads using hs->buf. */
void coreg_hist2_cached(const BaseSamples *bs, const Vol *VF, const double *P, int np,
                        double H[65536], Hist2Scratch *hs);

#endif
