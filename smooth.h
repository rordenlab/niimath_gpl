/* smooth.h — SPM spm_smoothkern + smooth_uint8 (separable conv) port. */
#ifndef SPMCOREG_SMOOTH_H
#define SPMCOREG_SMOOTH_H
#include "loaduint8.h"

/* Fill krn (caller-allocated, >= 2*ceil(2*fwhm)+1) for x = -lim..lim.
 * t=1: Gaussian (x) 1st-degree B-spline; t=0: 0th-degree. Returns length. */
int spm_smoothkern(double fwhm, double *krn, int t);

/* In-place separable smoothing of V->u8 with per-axis fwhm (in voxels),
 * matching SPM smooth_uint8 / spm_conv_vol (x,y truncated; z renormalized). */
void smooth_uint8(Vol *V, const double fwhm[3]);

#endif
