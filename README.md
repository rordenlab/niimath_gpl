# niimath_gpl

GPL-2 payload for niimath's optional `-spmcoreg` and `-spm_deface` commands — a C port of SPM's `spm_coreg` rigid-body coregistration (John Ashburner, Wellcome Centre for Human Neuroimaging). Used as the `src/GPL` git submodule of [niimath](https://github.com/rordenlab/niimath).

niimath itself is BSD-2-Clause. Building it with this module (`make GPL=1` / `cmake -DENABLE_GPL=ON`) produces a **GPL-2** combined-work binary; the default niimath build does not include this code and stays BSD-2-Clause.

The module computes only the rigid transform (estimate). niimath's own BSD code applies it (reslice / mask warp / header rewrite). Files: `cost hist2 loaduint8 matrix powell reslice smooth spm_coreg` (+ headers) and the niimath glue `spmcoreg_niimath`.

See `LICENSE` (GNU GPL-2) and `NOTICE` (SPM attribution).
