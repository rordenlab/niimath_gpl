/* reslice.c — trilinear/NN reslice into ref grid (SPM spm_reslice interp=1). */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "reslice.h"
#include "nifti_io.h"
#include "matrix.h"
#include "loaduint8.h"   /* vol_raw_at, vol_datatype_ok, spm_vox_mat */

#define TINY 5e-2

/* SPM bsplines.c mirror(): reflect 0-based index into [0,m-1]. */
static int mirror(int i,int m){ i=abs(i); if(i<m)return i; if(m==1)return 0;
    int m2=(m-1)*2; i%=m2; return (i<m)?i:m2-i; }

/* Reslice from an already-loaded source image (no file I/O, no ownership of nim). */
float *reslice_trilinear_nim(const Vol *VG, const nifti_image *nim, const double *P,
                             int np, int nn, int fill_zero) {
    float fillv = fill_zero ? 0.0f : NAN;
    int dt=nim->datatype;
    if (!vol_datatype_ok(dt) || nim->nx<=0 || nim->ny<=0 || nim->nz<=0 ||
        nim->nx>INT_MAX || nim->ny>INT_MAX || nim->nz>INT_MAX) {
        fprintf(stderr,"reslice: source unsupported datatype/dimensions\n"); return NULL; }
    int sx=nim->nx, sy=nim->ny, sz=nim->nz;
    long snv=(long)sx*sy*sz;
    if ((long long)snv!=(long long)nim->nvox) {
        fprintf(stderr,"reslice: source is 4D/multi-volume\n"); return NULL; }
    if (!isfinite(nim->scl_slope) || !isfinite(nim->scl_inter)) {
        fprintf(stderr,"reslice: non-finite source scl_slope/scl_inter\n"); return NULL; }
    double slope=(nim->scl_slope==0.0)?1.0:nim->scl_slope, inter=nim->scl_inter;
    double *S=malloc(sizeof(double)*snv);                 /* scaled source */
    if(!S) return NULL;
    for(long k=0;k<snv;k++) S[k]=vol_raw_at(nim->data,dt,k)*slope+inter;
    double VFmat[4][4];
    if (spm_vox_mat(nim,VFmat)) { free(S); return NULL; }

    /* M = inv(VFmat) * spm_matrix(P) * VG.mat  (ref vox -> src vox) */
    double A[4][4],VFi[4][4],tmp[4][4],M[4][4];
    spm_matrix(P,np,A);
    if (mat44_inv(VFmat,VFi)) { fprintf(stderr,"reslice: singular source affine\n"); free(S); return NULL; }
    mat44_mul(VFi,A,tmp); mat44_mul(tmp,VG->mat,M);

    int nx=VG->nx, ny=VG->ny, nz=VG->nz;
    float *out=malloc(sizeof(float)*VG->nvox);
    if(!out){ free(S); return NULL; }
    for(int iz=0;iz<nz;iz++) for(int iy=0;iy<ny;iy++) for(int ix=0;ix<nx;ix++){
        double rx=ix+1, ry=iy+1, rz=iz+1;
        double px=M[0][0]*rx+M[0][1]*ry+M[0][2]*rz+M[0][3];
        double py=M[1][0]*rx+M[1][1]*ry+M[1][2]*rz+M[1][3];
        double pz=M[2][0]*rx+M[2][1]*ry+M[2][2]*rz+M[2][3];
        long o=(long)ix+(long)nx*(iy+(long)ny*iz);
        if (px<1-TINY||px>sx+TINY||py<1-TINY||py>sy+TINY||pz<1-TINY||pz>sz+TINY){ out[o]=fillv; continue; }
        if (nn) {
            int jx=mirror((int)floor(px+0.5)-1,sx), jy=mirror((int)floor(py+0.5)-1,sy), jz=mirror((int)floor(pz+0.5)-1,sz);
            out[o]=(float)S[jx+(long)sx*(jy+(long)sy*jz)];
            continue;
        }
        int x0=(int)floor(px), y0=(int)floor(py), z0=(int)floor(pz);
        double fx=px-x0, fy=py-y0, fz=pz-z0;
        int x0i=mirror(x0-1,sx), x1i=mirror(x0,sx);   /* 0-based neighbours, mirror boundary */
        int y0i=mirror(y0-1,sy), y1i=mirror(y0,sy);
        int z0i=mirror(z0-1,sz), z1i=mirror(z0,sz);
        #define V(xi,yi,zi) S[(xi)+(long)sx*((yi)+(long)sy*(zi))]
        double c00=V(x0i,y0i,z0i)*(1-fx)+V(x1i,y0i,z0i)*fx;
        double c10=V(x0i,y1i,z0i)*(1-fx)+V(x1i,y1i,z0i)*fx;
        double c01=V(x0i,y0i,z1i)*(1-fx)+V(x1i,y0i,z1i)*fx;
        double c11=V(x0i,y1i,z1i)*(1-fx)+V(x1i,y1i,z1i)*fx;
        #undef V
        double c0=c00*(1-fy)+c10*fy, c1=c01*(1-fy)+c11*fy;
        out[o]=(float)(c0*(1-fz)+c1*fz);
    }
    free(S);
    return out;
}

float *reslice_trilinear(const Vol *VG, const char *src, const double *P, int np,
                         int nn, int fill_zero) {
    nifti_image *nim=nifti_image_read(src,1); if(!nim) return NULL;
    float *out=reslice_trilinear_nim(VG,nim,P,np,nn,fill_zero);
    nifti_image_free(nim);
    return out;
}

#ifdef RESLICE_TEST
#include <stdio.h>
#include <string.h>
#include <assert.h>
int main(void){
    Vol VG; if(load_uint8("golden/6_anat-T1w_acq-tfl.nii",&VG,0)){fprintf(stderr,"load fail\n");return 2;}
    FILE *fp=fopen("golden/x_final.txt","r"); char tag[16]; double xk[6];
    fscanf(fp,"%15s",tag); for(int i=0;i<6;i++) fscanf(fp,"%lf",&xk[i]); fclose(fp);
    float *r=reslice_trilinear(&VG,"golden/avg152T1.nii",xk,6,0,0); /* NaN fill: matches SPM golden */
    fp=fopen("golden/ravg_trilin.bin","rb"); float *g=malloc(sizeof(float)*VG.nvox);
    long n=fread(g,sizeof(float),VG.nvox,fp); fclose(fp); assert(n==VG.nvox);
    long both_nan=0, nan_mismatch=0, finite=0; double worst=0;
    for(long k=0;k<VG.nvox;k++){
        int cn=isnan(r[k]), gn=isnan(g[k]);
        if(cn&&gn){both_nan++; continue;}
        if(cn!=gn){nan_mismatch++; continue;}
        finite++; double d=fabs((double)r[k]-(double)g[k]); if(d>worst)worst=d;
    }
    free(r);free(g);vol_free(&VG);
    printf("reslice: finite=%ld both_NaN=%ld NaN_mismatch=%ld  max|C-golden|=%.3e\n",
           finite,both_nan,nan_mismatch,worst);
    assert(nan_mismatch==0 && worst<1e-5);
    printf("M8 GATE PASS\n");
    return 0;
}
#endif
