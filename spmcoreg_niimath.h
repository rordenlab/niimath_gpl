/* spmcoreg_niimath.h — niimath chain-op glue for the GPL spm_coreg module.
 * Part of the GPL-2 payload (niimath_gpl); only built with -DHAVE_GPL. */
#ifndef NII_SPMCOREG_H
#define NII_SPMCOREG_H
#include "nifti_io.h"

#ifdef __cplusplus
extern "C" {
#endif

/* `-spm_coreg ref [opts]` chain op. argv[*ac] is the "-spm_coreg" token on entry;
 * parses `ref` + sub-options and advances *ac past everything consumed. Registers
 * the in-memory `nim` (moving) to `ref` (stationary): reslices onto the ref grid by
 * default, or rewrites nim's sform/qform with `-estimate`. `fin` is the source path
 * (fallback only). `dt_double` nonzero means the chain runs in -dt double (main64):
 * reslice then requires -dt float and is rejected; -estimate works either way.
 * Returns 0 on success, nonzero on error. */
int nii_spmcoreg(nifti_image *nim, char *fin, int *ac, int argc, char *argv[], int dt_double);

/* `-spm_deface tmpl mask [opts]` chain op. Like -deface but registers with SPM
 * spm_coreg; `nim` is the image being defaced (fixed reference, stays in its grid).
 * argv[*ac] is the "-spm_deface" token on entry; advances *ac past what it consumes.
 * `dt_double` nonzero (-dt double / main64) rejects deface (it writes float32 data).
 * Returns 0 on success, nonzero on error. */
int nii_spm_deface(nifti_image *nim, int *ac, int argc, char *argv[], int dt_double);

#ifdef __cplusplus
}
#endif

#endif
