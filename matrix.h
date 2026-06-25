/* matrix.h — SPM spm_matrix port + 4x4 affine math (double, row-major). */
#ifndef SPMCOREG_MATRIX_H
#define SPMCOREG_MATRIX_H

/* Build SPM affine A = T*R*Z*S from parameter vector P of length n (n in 1..12).
 * P is padded with [0 0 0 0 0 0 1 1 1 0 0 0] to 12, matching spm_matrix.m. */
void spm_matrix(const double *P, int n, double A[4][4]);

/* C = A * B (4x4). */
void mat44_mul(const double A[4][4], const double B[4][4], double C[4][4]);

/* Inv = inverse(A) via Gauss-Jordan with partial pivoting. Returns 0 on success. */
int  mat44_inv(const double A[4][4], double Inv[4][4]);

#endif
