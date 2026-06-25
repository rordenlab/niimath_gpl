/* powell.h — SPM spm_powell.m port (Powell + Brent line search). */
#ifndef SPMCOREG_POWELL_H
#define SPMCOREG_POWELL_H

/* Objective: returns f(p), p has length n. ctx is user data. */
typedef double (*objfun)(const double *p, int n, void *ctx);

/* Minimise f from starting point p (length n, updated in place).
 * xi: n*n initial direction set, column-major (direction j = xi[j*n + 0..n-1]).
 * tolsc: per-parameter stopping scale. Returns f at the minimum. */
double spm_powell(double *p, int n, double *xi, const double *tolsc,
                  objfun f, void *ctx);

#endif
