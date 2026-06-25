/* precision.h — smooth-stage buffer precision (`sflt`).
 *
 * Default = double (faithful to SPM, and FASTER on real 1mm data). float-smooth
 * was tested: its output is uint8, but the +/-1 rounding shifts the smoothed
 * volumes enough to move the coarse minimum and disturb the coarse-to-fine
 * handoff, costing ~25% more optimizer evals (480 vs 378 on the large pair) —
 * verified, not assumed. -DSMOOTH32 makes `sflt` = float: it cuts the smoothing
 * buffers (the peak-RAM hog on large volumes, ~37%) at that eval-count/speed cost.
 * Use -DSMOOTH32 only when memory-constrained on large data.
 *
 * NOTE: histogram (H), cost, transform params and 4x4 matrices are ALWAYS double. */
#ifndef SPMCOREG_PRECISION_H
#define SPMCOREG_PRECISION_H

#ifdef SMOOTH32
typedef float sflt;
#else
typedef double sflt;
#endif

#endif
