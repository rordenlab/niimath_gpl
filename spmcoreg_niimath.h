/* spmcoreg_niimath.h — niimath chain-op glue for the GPL spm_coreg module.
 * Part of the GPL-2 payload (niimath_gpl); only built with -DHAVE_GPL. */
#ifndef NII_SPMCOREG_H
#define NII_SPMCOREG_H
#include "nifti_io.h"

#ifdef __cplusplus
extern "C" {
#endif

/* `-spmcoreg ref [opts]` chain op. argv[*ac] is the "-spmcoreg" token on entry;
 * parses `ref` + sub-options and advances *ac past everything consumed. Registers
 * the in-memory `nim` (moving) to `ref` (stationary): reslices onto the ref grid by
 * default, or rewrites nim's sform/qform with `-estimate`. `fin` is the source path
 * (fallback only). Returns 0 on success, nonzero on error. */
int nii_spmcoreg(nifti_image *nim, char *fin, int *ac, int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
