/* spmcoreg_niimath.c — niimath chain-op glue for the GPL spm_coreg module.
 * Part of the GPL-2 payload (niimath_gpl); only built with -DHAVE_GPL.
 *
 * The GPL module only computes the rigid transform (estimate). niimath's BSD
 * code applies it: reslice via nii_reslice_affine() (shared with -allineate);
 * -estimate rewrites the header only (Phase 3b). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spmcoreg_niimath.h"
#include "spm_coreg.h"
#include "matrix.h"
#include "loaduint8.h"   /* spm_vox_mat */
#ifdef HAVE_ALLINEATE
#include "allineate.h"   /* nii_reslice_affine, mat44, AL_INTERP_* (BSD reslice) */
#endif

/* Recover the 0-based NIfTI affine (index->world) from spm_vox_mat's 1-based M:
 * SPM index = NIfTI index + 1, so the 3x3 is identical and the translation gains
 * back the column sums that spm_vox_mat subtracted. */
static void nii0_from_spm(const double Mspm[4][4], double Mnii[4][4]) {
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 3; c++) Mnii[r][c] = Mspm[r][c];
        Mnii[r][3] = Mspm[r][3] + (Mspm[r][0] + Mspm[r][1] + Mspm[r][2]);
    }
}

/* C = A (double 4x4) * B (nifti_dmat44). */
static nifti_dmat44 premul(const double A[4][4], nifti_dmat44 B) {
    nifti_dmat44 C;
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++) {
            double s = 0; for (int k = 0; k < 4; k++) s += A[r][k] * B.m[k][c];
            C.m[r][c] = s;
        }
    return C;
}

/* parse a flag value as a finite double; pos>0, or >=0 if allow_zero. 0 on success. */
static int parse_num(const char *s, double *out, int allow_zero) {
    char *end; double v = strtod(s, &end);
    if (end == s || *end != '\0' || !isfinite(v) || v < 0.0 || (!allow_zero && v == 0.0)) return 1;
    *out = v; return 0;
}
/* strict integer parse restricted to {0,1}. 0 on success. */
static int parse_bit(const char *s, int *out) {
    char *end; long v = strtol(s, &end, 10);
    if (end == s || *end != '\0' || (v != 0 && v != 1)) return 1;
    *out = (int)v; return 0;
}

int nii_spmcoreg(nifti_image *nim, char *fin, int *ac, int argc, char *argv[]) {
    (void)fin;
    /* argv[*ac] == "-spmcoreg"; next token is the (required) reference image. */
    if (*ac + 1 >= argc) { fprintf(stderr, "-spmcoreg requires a reference image argument\n"); return 1; }
    (*ac)++;
    const char *ref = argv[*ac];

    coreg_flags f; coreg_default_flags(&f);
    int verbose = 0, estimate = 0;
    int nn = 0;          /* 0 = trilinear (interp=1), 1 = nearest  (Phase 3) */
    int fill_zero = 1;   /* out-of-FOV: 1 = zero (FSL-style), 0 = NaN (Phase 3) */

    /* sub-options: consume while the next top-level token is a recognized flag.
     * value tokens (sep/fwhm numbers, names) are grabbed explicitly inside each case. */
    while (*ac + 1 < argc && argv[*ac + 1][0] == '-') {
        (*ac)++;
        const char *o = argv[*ac];
        if (!strcmp(o, "-cost")) {
            if (*ac + 1 >= argc || cost_parse(argv[++(*ac)], &f.cf)) { fprintf(stderr, "-cost must be nmi|mi|ecc|ncc|ls\n"); return 1; }
        } else if (!strcmp(o, "-sep")) {
            f.nsep = 0;
            while (*ac + 1 < argc && f.nsep < 4) { double v; if (parse_num(argv[*ac + 1], &v, 0)) break; f.sep[f.nsep++] = v; (*ac)++; }
            if (f.nsep < 1) { fprintf(stderr, "-sep needs >=1 positive value\n"); return 1; }
        } else if (!strcmp(o, "-fwhm")) {
            if (*ac + 2 >= argc || parse_num(argv[++(*ac)], &f.fwhm[0], 1) || parse_num(argv[++(*ac)], &f.fwhm[1], 1)) {
                fprintf(stderr, "-fwhm needs two values >=0\n"); return 1; }
        } else if (!strcmp(o, "-dither")) {
            if (*ac + 1 >= argc || parse_bit(argv[++(*ac)], &f.dither)) { fprintf(stderr, "-dither must be 0 or 1\n"); return 1; }
        } else if (!strcmp(o, "-coarse")) {
            const char *v = (*ac + 1 < argc) ? argv[++(*ac)] : "";
            if (!strcmp(v, "downsample")) f.coarse_downsample = 1;
            else if (!strcmp(v, "sparse")) f.coarse_downsample = 0;
            else { fprintf(stderr, "-coarse must be sparse|downsample\n"); return 1; }
        } else if (!strcmp(o, "-verbose")) {
            if (*ac + 1 >= argc || parse_bit(argv[++(*ac)], &verbose)) { fprintf(stderr, "-verbose must be 0 or 1\n"); return 1; }
        } else if (!strcmp(o, "-estimate")) {
            estimate = 1;
        } else if (!strcmp(o, "-interp")) {
            const char *v = (*ac + 1 < argc) ? argv[++(*ac)] : "";
            if (!strcmp(v, "trilinear")) nn = 0;
            else if (!strcmp(v, "nearest")) nn = 1;
            else { fprintf(stderr, "-interp must be trilinear|nearest\n"); return 1; }
        } else if (!strcmp(o, "-fill")) {
            const char *v = (*ac + 1 < argc) ? argv[++(*ac)] : "";
            if (!strcmp(v, "zero")) fill_zero = 1;
            else if (!strcmp(v, "nan")) fill_zero = 0;
            else { fprintf(stderr, "-fill must be zero|nan\n"); return 1; }
        } else {
            (*ac)--;   /* not ours: hand the token back to the niimath chain */
            break;
        }
    }
    nifti_image *refnim = nifti_image_read(ref, 1);
    if (!refnim) { fprintf(stderr, "-spmcoreg: cannot read reference %s\n", ref); return 1; }
    double xk[6];
    int rc = coreg_estimate_nim(refnim, nim, &f, xk);   /* ref=stationary, nim=moving */
    if (rc) { fprintf(stderr, "-spmcoreg: estimate failed\n"); nifti_image_free(refnim); return 1; }

    if (verbose) {
        printf("params  %.6f %.6f %.6f  %.6f %.6f %.6f\n", xk[0], xk[1], xk[2], xk[3], xk[4], xk[5]);
        double A[4][4]; spm_matrix(xk, 6, A);
        printf("spm_matrix(params) [world transform applied to source]:\n");
        for (int r = 0; r < 4; r++) printf("  %9.5f %9.5f %9.5f %9.5f\n", A[r][0], A[r][1], A[r][2], A[r][3]);
    }

    /* -estimate: header-only. Keep source data/dims; left-multiply its world
     * affines by inv(spm_matrix(xk)) (SPM's spm_get_space(VF, spm_matrix(xk)\VF.mat)). */
    if (estimate) {
        double A[4][4], Ainv[4][4];
        spm_matrix(xk, 6, A);
        if (mat44_inv(A, Ainv)) { fprintf(stderr, "-spmcoreg: singular transform\n"); nifti_image_free(refnim); return 1; }
        nifti_image_free(refnim);
        if (nim->sform_code <= 0 && nim->qform_code <= 0) {
            fprintf(stderr, "-spmcoreg -estimate: source has no sform/qform to update\n"); return 1; }
        if (nim->sform_code > 0) {
            nim->sto_xyz = premul(Ainv, nim->sto_xyz);
            nim->sto_ijk = nifti_dmat44_inverse(nim->sto_xyz);
        }
        if (nim->qform_code > 0) {
            nim->qto_xyz = premul(Ainv, nim->qto_xyz);
            nim->qto_ijk = nifti_dmat44_inverse(nim->qto_xyz);
            double dx, dy, dz, qfac;
            nifti_dmat44_to_quatern(nim->qto_xyz, &nim->quatern_b, &nim->quatern_c, &nim->quatern_d,
                                    &nim->qoffset_x, &nim->qoffset_y, &nim->qoffset_z, &dx, &dy, &dz, &qfac);
            nim->pixdim[0] = qfac;   /* rigid -> dx,dy,dz unchanged; only orientation/offset move */
        }
        return 0;   /* data + dims untouched */
    }

    /* Apply the transform. Default: reslice the source onto the ref grid using
     * BSD interpolation; voxel->voxel (0-based) gam = inv(VF)*spm_matrix(xk)*VG. */
#ifdef HAVE_ALLINEATE
    double VFs[4][4], VGs[4][4];
    if (spm_vox_mat(nim, VFs) || spm_vox_mat(refnim, VGs)) {
        fprintf(stderr, "-spmcoreg: source/ref lacks a usable affine\n"); nifti_image_free(refnim); return 1; }
    double VFn[4][4], VGn[4][4]; nii0_from_spm(VFs, VFn); nii0_from_spm(VGs, VGn);
    double A[4][4], VFi[4][4], t[4][4], M[4][4];
    spm_matrix(xk, 6, A);
    if (mat44_inv(VFn, VFi)) { fprintf(stderr, "-spmcoreg: singular source affine\n"); nifti_image_free(refnim); return 1; }
    mat44_mul(VFi, A, t); mat44_mul(t, VGn, M);
    mat44 gam; memset(&gam, 0, sizeof gam); gam.m[3][3] = 1.0f;
    for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) gam.m[r][c] = (float)M[r][c];
    int interp = nn ? AL_INTERP_NN : AL_INTERP_LINEAR;
    float fillv = fill_zero ? 0.0f : (float)NAN;
    rc = nii_reslice_affine(nim, refnim, gam, interp, fillv);
    nifti_image_free(refnim);
    if (rc) { fprintf(stderr, "-spmcoreg: reslice failed\n"); return 1; }
#else
    (void)nn; (void)fill_zero;
    nifti_image_free(refnim);
    fprintf(stderr, "-spmcoreg: reslice requires the allineate module (build without AL=0); use -estimate instead\n");
    return 1;
#endif
    return 0;
}
