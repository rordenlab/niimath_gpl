/* reslice.h — trilinear reslice (SPM spm_reslice interp=1) port. */
#ifndef SPMCOREG_RESLICE_H
#define SPMCOREG_RESLICE_H
#include "loaduint8.h"

/* Resample source file into the ref grid using params P (length np), trilinear.
 * Voxel->voxel M = inv(VF.mat)*spm_matrix(P)*VG.mat (same as coreg).
 * Returns malloc'd float[ref nvox] (x fastest). Out-of-FOV voxels = 0 if fill_zero,
 * else NaN (matching spm_bsplins degree 1). VG provides the output grid/mat.
 * nn!=0 => nearest-neighbour. */
float *reslice_trilinear(const Vol *VG, const char *src, const double *P, int np,
                         int nn, int fill_zero);

/* Same reslice from an already-loaded source image (niimath in-memory chain).
 * Does not take ownership of nim. */
float *reslice_trilinear_nim(const Vol *VG, const nifti_image *nim, const double *P,
                             int np, int nn, int fill_zero);

#endif
