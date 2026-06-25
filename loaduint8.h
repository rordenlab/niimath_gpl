/* loaduint8.h — SPM loaduint8 port: read NIfTI -> robust-max -> uint8 volume. */
#ifndef SPMCOREG_LOADUINT8_H
#define SPMCOREG_LOADUINT8_H
#include "nifti_io.h"

typedef struct {
    int nx, ny, nz;
    long nvox;
    double mat[4][4];   /* SPM V.mat (1-based voxel -> world) */
    double vox[3];      /* voxel sizes = sqrt(sum(mat(1:3,1:3).^2)) per column */
    unsigned char *u8;  /* nvox bytes, x fastest */
} Vol;

/* When dither != 0, add reproducible sub-quantization noise (production).
 * dither == 0 matches the dither-off golden used for validation. */
int load_uint8(const char *fname, Vol *V, int dither);

/* Same quantization on an already-loaded image (niimath in-memory chain adapter).
 * Does not take ownership of nim; `tag` labels diagnostics. Returns 0 on success. */
int load_uint8_nim(const nifti_image *nim, const char *tag, Vol *V, int dither);

void vol_free(Vol *V);

/* Shared geometry/datatype helpers (also used by reslice.c, main.c). */
double vol_raw_at(const void *data, int datatype, long k); /* raw voxel -> double */
int    vol_datatype_ok(int datatype);                      /* 1 if a supported scalar type */
/* SPM V.mat (1-based) from sform, falling back to qform. Returns 0 ok, 1 = no usable affine. */
int    spm_vox_mat(const nifti_image *nim, double M[4][4]);

/* Block-average downsample by per-axis integer factor F (>=1). Builds a smaller
 * Vol with the correct SPM affine/voxel sizes (block centres). Returns 0 on
 * success (out owns new u8; vol_free it), 1 if a factor leaves <1 voxel. */
int    vol_downsample(const Vol *in, const int F[3], Vol *out);

#endif
