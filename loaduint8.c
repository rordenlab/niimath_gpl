/* loaduint8.c — faithful port of SPM spm_coreg.m:loaduint8 (general path).
 * Robust max via 2048-bin histogram (cumsum/sum > 0.9999), scale to uint8.
 * Dither off (validation) is byte-exact vs SPM with the random term zeroed. */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "loaduint8.h"
#include "nifti_io.h"

static int is_int_type(int dt) {
    switch (dt) {
        case 2: case 4: case 8: case 256: case 512: case 768:
        case 1024: case 1280: return 1;  /* u8,i16,i32,i8,u16,u32,i64,u64 */
        default: return 0;               /* float32/64 -> 0 */
    }
}

/* Supported scalar real datatypes (reject complex/RGB/float128/etc.). */
int vol_datatype_ok(int dt) {
    switch (dt) {
        case 2: case 4: case 8: case 16: case 64:
        case 256: case 512: case 768: case 1024: case 1280: return 1;
        default: return 0;
    }
}

/* read raw voxel k as double (no scaling). Caller must vol_datatype_ok(dt) first. */
double vol_raw_at(const void *d, int dt, long k) {
    switch (dt) {
        case 2:    return ((const unsigned char *)d)[k];
        case 4:    return ((const short *)d)[k];
        case 8:    return ((const int *)d)[k];
        case 16:   return ((const float *)d)[k];
        case 64:   return ((const double *)d)[k];
        case 256:  return ((const signed char *)d)[k];
        case 512:  return ((const unsigned short *)d)[k];
        case 768:  return ((const unsigned int *)d)[k];
        case 1024: return (double)((const int64_t *)d)[k];
        case 1280: return (double)((const uint64_t *)d)[k];
        default:   return 0.0;  /* unreachable: validated upstream */
    }
}

/* An affine is usable if every entry is finite, each voxel-axis column has
 * nonzero length, and the 3x3 is non-singular (invertible). */
static int affine_usable(const double s[4][4]) {
    for (int r=0;r<4;r++) for (int c=0;c<4;c++) if (!isfinite(s[r][c])) return 0;
    for (int c=0;c<3;c++) {
        double n = sqrt(s[0][c]*s[0][c]+s[1][c]*s[1][c]+s[2][c]*s[2][c]);
        if (!(n > 1e-6)) return 0;
    }
    double det = s[0][0]*(s[1][1]*s[2][2]-s[1][2]*s[2][1])
               - s[0][1]*(s[1][0]*s[2][2]-s[1][2]*s[2][0])
               + s[0][2]*(s[1][0]*s[2][1]-s[1][1]*s[2][0]);
    return isfinite(det) && fabs(det) > 1e-12;
}

/* SPM V.mat (1-based voxel -> world). Prefer a USABLE sform; fall back to a
 * usable qform. Returns 0 on success, 1 if neither is present and usable. */
int spm_vox_mat(const nifti_image *nim, double M[4][4]) {
    const double (*s)[4] = NULL;
    if (nim->sform_code > 0 && affine_usable(nim->sto_xyz.m))      s = nim->sto_xyz.m;
    else if (nim->qform_code > 0 && affine_usable(nim->qto_xyz.m)) s = nim->qto_xyz.m;
    else return 1;
    for (int r=0;r<4;r++) {
        for (int c=0;c<3;c++) M[r][c]=s[r][c];
        M[r][3]=s[r][3]-(s[r][0]+s[r][1]+s[r][2]);
    }
    return 0;
}

/* Quantize an already-loaded nifti_image to a uint8 Vol (no file I/O, no
 * ownership of nim — caller frees). Shared by load_uint8 (file path) and the
 * niimath in-memory chain adapter. `tag` labels diagnostics. */
int load_uint8_nim(const nifti_image *nim, const char *tag, Vol *V, int dither) {
    int dt = nim->datatype;
    /* reject what we cannot register: unsupported datatype, oversized/empty dims, 4D+, no usable affine */
    if (!vol_datatype_ok(dt)) {
        fprintf(stderr,"%s: unsupported datatype %d (scalar real types only)\n",tag,dt);
        return 1; }
    if (nim->nx<=0 || nim->ny<=0 || nim->nz<=0 ||
        nim->nx>INT_MAX || nim->ny>INT_MAX || nim->nz>INT_MAX) {
        fprintf(stderr,"%s: dimensions out of range\n",tag);
        return 1; }
    long nvox = (long)nim->nx*nim->ny*nim->nz;
    if ((long long)nvox != (long long)nim->nvox) {
        fprintf(stderr,"%s: 4D/multi-volume input not supported (use a single 3D volume)\n",tag);
        return 1; }
    if (!isfinite(nim->scl_slope) || !isfinite(nim->scl_inter)) {
        fprintf(stderr,"%s: non-finite scl_slope/scl_inter\n",tag);
        return 1; }
    double slope = (nim->scl_slope==0.0) ? 1.0 : nim->scl_slope;
    double inter = nim->scl_inter;
    double acc   = is_int_type(dt) ? fabs(slope) : 0.0;
    const void *d = nim->data;

    /* pass 1: global max/min of scaled finite voxels; mx gets +acc once */
    double mx=-INFINITY, mn=INFINITY; long nfin=0;
    for (long k=0;k<nvox;k++) {
        double f = vol_raw_at(d,dt,k)*slope+inter;
        if (!isfinite(f)) continue;
        nfin++;
        if (f>mx) mx=f; if (f<mn) mn=f;
    }
    mx += acc;
    if (nfin==0 || !isfinite(mx) || !isfinite(mn) || mx<=mn) {
        fprintf(stderr,"%s: image has no finite dynamic range (empty/constant)\n",tag);
        return 1; }

    /* pass 2: 2048-bin histogram */
    const int nh=2048;
    long *h = calloc(nh,sizeof(long));
    if (!h) return 1;
    double a = (mx-mn)/(nh-1);
    double b = (nh-1)/(mx-mn);
    long ntot=0;
    for (long k=0;k<nvox;k++) {
        double f = vol_raw_at(d,dt,k)*slope+inter;
        if (!isfinite(f)) continue;
        int idx = (int)round((f + (a - mn))*b);   /* MATLAB round = half away from 0 */
        if (idx<1) idx=1; if (idx>nh) idx=nh;      /* guard (theory: in [1,nh]) */
        h[idx-1]++; ntot++;
    }
    /* first bin where cumsum/sum > 0.9999 (1-based); fallback nh */
    long cum=0; int tmp1=nh;
    double thr=0.9999*(double)ntot;
    for (int i=0;i<nh;i++){ cum+=h[i]; if ((double)cum>thr){ tmp1=i+1; break; } }
    free(h);
    /* mx = (mn*nh - mx + tmp1*(mx-mn))/(nh-1)  (mx on RHS = old mx) */
    mx = (mn*nh - mx + tmp1*(mx-mn))/(nh-1);

    /* final: scale to uint8 */
    unsigned char *u8 = malloc(nvox);
    if (!u8) return 1;
    double sc = 255.0/(mx-mn);
    for (long k=0;k<nvox;k++) {
        double f = vol_raw_at(d,dt,k)*slope+inter;
        double r = 0.0;
        if (dither && acc!=0.0) {
            /* deterministic sub-quantization dither in [0,acc): mask 21 bits -> [0,1) */
            unsigned long s = (unsigned long)k*2654435761UL;
            r = (double)((s>>11) & 0x1FFFFFUL)/(double)(1UL<<21) * acc;
        }
        double v = round((f+r-mn)*sc);
        if (v<0) v=0; if (v>255) v=255;
        u8[k]=(unsigned char)v;
    }

    if (spm_vox_mat(nim, V->mat)) {
        fprintf(stderr,"%s: no usable sform/qform affine\n",tag);
        free(u8); return 1; }
    for (int c=0;c<3;c++)
        V->vox[c]=sqrt(V->mat[0][c]*V->mat[0][c]+V->mat[1][c]*V->mat[1][c]+V->mat[2][c]*V->mat[2][c]);
    V->nx=nim->nx; V->ny=nim->ny; V->nz=nim->nz; V->nvox=nvox; V->u8=u8;
    return 0;
}

int load_uint8(const char *fname, Vol *V, int dither) {
    nifti_image *nim = nifti_image_read(fname, 1);
    if (!nim) return 1;
    int r = load_uint8_nim(nim, fname, V, dither);
    nifti_image_free(nim);
    return r;
}

void vol_free(Vol *V){ if(V&&V->u8){ free(V->u8); V->u8=NULL; } }

int vol_downsample(const Vol *in, const int F[3], Vol *out) {
    int Fx=F[0]<1?1:F[0], Fy=F[1]<1?1:F[1], Fz=F[2]<1?1:F[2];
    /* trailing voxels not divisible by the factor are cropped (block decimation,
     * as SPM does); only used for the opt-in -coarse downsample fast path. */
    int ncx=in->nx/Fx, ncy=in->ny/Fy, ncz=in->nz/Fz;
    if (ncx<1 || ncy<1 || ncz<1) return 1;
    long nx=in->nx, ny=in->ny, ncv=(long)ncx*ncy*ncz, blk=(long)Fx*Fy*Fz;
    unsigned char *u = malloc(ncv);
    if (!u) return 1;
    for (int K=0;K<ncz;K++) for (int J=0;J<ncy;J++) for (int I=0;I<ncx;I++) {
        long sum=0;
        for (int kk=0;kk<Fz;kk++) for (int jj=0;jj<Fy;jj++) for (int ii=0;ii<Fx;ii++)
            sum += in->u8[((long)I*Fx+ii) + nx*(((long)J*Fy+jj) + ny*(long)(K*Fz+kk))];
        u[I + (long)ncx*(J + (long)ncy*K)] = (unsigned char)((double)sum/blk + 0.5);
    }
    /* new SPM affine = M * S, where S maps coarse 1-based index i' to old index
     * F*i' - (F-1)/2 (block centre). Columns scale by F; offset shifts by -(F-1)/2. */
    for (int r=0;r<4;r++) {
        out->mat[r][0]=in->mat[r][0]*Fx;
        out->mat[r][1]=in->mat[r][1]*Fy;
        out->mat[r][2]=in->mat[r][2]*Fz;
        out->mat[r][3]=in->mat[r][3]
            - (in->mat[r][0]*(Fx-1) + in->mat[r][1]*(Fy-1) + in->mat[r][2]*(Fz-1))/2.0;
    }
    for (int c=0;c<3;c++)
        out->vox[c]=sqrt(out->mat[0][c]*out->mat[0][c]+out->mat[1][c]*out->mat[1][c]+out->mat[2][c]*out->mat[2][c]);
    out->nx=ncx; out->ny=ncy; out->nz=ncz; out->nvox=ncv; out->u8=u;
    return 0;
}

#ifdef LOADUINT8_TEST
#include <stdio.h>
#include <assert.h>
static int check(const char *nii, const char *golden, const char *tag){
    Vol V; if(load_uint8(nii,&V,0)){fprintf(stderr,"load fail %s\n",nii);return 2;}
    FILE *f=fopen(golden,"rb"); if(!f){fprintf(stderr,"no golden %s\n",golden);return 2;}
    unsigned char *g=malloc(V.nvox); long n=fread(g,1,V.nvox,f); fclose(f);
    if(n!=V.nvox){fprintf(stderr,"golden size %ld != %ld\n",n,V.nvox);return 2;}
    long diff=0,maxd=0; for(long k=0;k<V.nvox;k++){int dd=abs((int)V.u8[k]-(int)g[k]); if(dd){diff++; if(dd>maxd)maxd=dd;}}
    printf("%s: nvox=%ld mismatches=%ld maxabs=%ld\n",tag,V.nvox,diff,maxd);
    free(g); vol_free(&V);
    return diff==0?0:1;
}
/* Phase 1 adapter check: load_uint8_nim(read(f)) must be bit-identical to load_uint8(f). */
static int check_nim(const char *nii){
    Vol Vf, Vn;
    if(load_uint8(nii,&Vf,0)){fprintf(stderr,"file load fail %s\n",nii);return 2;}
    nifti_image *nim=nifti_image_read(nii,1);
    if(!nim||load_uint8_nim(nim,nii,&Vn,0)){fprintf(stderr,"nim load fail %s\n",nii);return 2;}
    nifti_image_free(nim);
    int ok = (Vf.nvox==Vn.nvox) && (memcmp(Vf.u8,Vn.u8,Vf.nvox)==0)
             && (memcmp(Vf.mat,Vn.mat,sizeof(Vf.mat))==0);
    printf("nim-adapter %s: %s\n",nii,ok?"identical":"MISMATCH");
    vol_free(&Vf); vol_free(&Vn);
    return ok?0:1;
}
int main(void){
    int r1=check("golden/6_anat-T1w_acq-tfl.nii","golden/vg_uint8_raw.bin","VG");
    int r2=check("golden/avg152T1.nii","golden/vf_uint8_raw.bin","VF");
    int r3=check_nim("golden/6_anat-T1w_acq-tfl.nii");
    int r4=check_nim("golden/avg152T1.nii");
    if(r1==0&&r2==0&&r3==0&&r4==0){printf("M2 GATE PASS\n");return 0;}
    printf("M2 GATE FAIL\n");return 1;
}
#endif
