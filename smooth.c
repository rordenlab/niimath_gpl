/* smooth.c — faithful port of spm_smoothkern.m + spm_conv_vol.c (separable). */
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "xalloc.h"
#include "smooth.h"
#include "precision.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int spm_smoothkern(double fwhm, double *krn, int t) {
    int lim = (int)ceil(2*fwhm);
    double s = (fwhm/sqrt(8*log(2.0))); s = s*s + DBL_EPSILON;
    int n = 2*lim+1;
    if (t==1) {
        double w1 = 0.5*sqrt(2/s), w2 = -0.5/s, w3 = sqrt(s/2/M_PI);
        for (int idx=0; idx<n; idx++) {
            double x = idx - lim;
            double k = 0.5*( erf(w1*(x+1))*(x+1) + erf(w1*(x-1))*(x-1) - 2*erf(w1*x)*x )
                     + w3*( exp(w2*(x+1)*(x+1)) + exp(w2*(x-1)*(x-1)) - 2*exp(w2*x*x) );
            krn[idx] = k<0 ? 0 : k;
        }
    } else { /* t==0 */
        double w1 = 1/sqrt(2*s);
        for (int idx=0; idx<n; idx++) {
            double x = idx - lim;
            double k = 0.5*(erf(w1*(x+0.5)) - erf(w1*(x-0.5)));
            krn[idx] = k<0 ? 0 : k;
        }
    }
    return n;
}

#ifdef _OPENMP
static void conv_axis_parallel(const sflt *in, sflt *out, int nx,int ny,int nz,
                               int axis, const double *krn, int klen, int koff, int renorm,
                               int L, int stride, int n_a, long nlines, int nt) {
    sflt *buffs = xmalloc(sizeof(sflt)*(size_t)L*(size_t)nt);
    #pragma omp parallel for schedule(static)
    for (long line=0; line<nlines; line++) {
        int b = (int)(line / n_a);
        int a = (int)(line - (long)b*n_a);
        sflt *buff = buffs + (size_t)omp_get_thread_num()*(size_t)L;
        long base;
        if (axis==0)      base = (long)a*nx + (long)b*nx*ny;       /* a=y,b=z */
        else if (axis==1) base = a + (long)b*nx*ny;                /* a=x,b=z */
        else              base = a + (long)b*nx;                   /* a=x,b=y */
        for (int o=0;o<L;o++) buff[o]=in[base+(long)o*stride];
        for (int o=0;o<L;o++) {
            int fstart = ((o-koff >= L) ? o-L-koff+1 : 0);
            int fend   = ((o-(koff+klen) < 0) ? o-koff+1 : klen);
            double sum1=0.0, sum2=0.0;
            for (int k=fstart;k<fend;k++){ sum1 += buff[o-koff-k]*krn[k]; if(renorm) sum2+=krn[k]; }
            out[base+(long)o*stride] = renorm ? (sum2!=0.0 ? sum1/sum2 : 0.0) : sum1;
        }
    }
    free(buffs);
}
#endif

/* 1D conv along one axis, replicating spm_conv_vol's convxy(x,y; renorm=0)
 * and convxyz z-pass (renorm=1). koff = -lim. */
static void conv_axis(const sflt *in, sflt *out, int nx,int ny,int nz,
                      int axis, const double *krn, int klen, int koff, int renorm) {
    int L = (axis==0)?nx : (axis==1)?ny : nz;
    int stride = (axis==0)?1 : (axis==1)?nx : nx*ny;
    int n_a = (axis==0)?ny:nx;        /* outer dims to iterate */
    int n_b = (axis==2)?ny:nz;
#ifdef _OPENMP
    int nt = omp_get_max_threads();
    long nlines = (long)n_a*n_b;
    if (nt > 1 && nlines > 64) {
        conv_axis_parallel(in,out,nx,ny,nz,axis,krn,klen,koff,renorm,L,stride,n_a,nlines,nt);
        return;
    }
#endif
    sflt *buff = xmalloc(sizeof(sflt)*L);
    for (int b=0;b<n_b;b++) {
        for (int a=0;a<n_a;a++) {
            /* base offset of this line */
            long base;
            if (axis==0)      base = (long)a*nx + (long)b*nx*ny;       /* a=y,b=z */
            else if (axis==1) base = a + (long)b*nx*ny;                /* a=x,b=z */
            else              base = a + (long)b*nx;                   /* a=x,b=y */
            for (int o=0;o<L;o++) buff[o]=in[base+(long)o*stride];
            for (int o=0;o<L;o++) {
                int fstart = ((o-koff >= L) ? o-L-koff+1 : 0);
                int fend   = ((o-(koff+klen) < 0) ? o-koff+1 : klen);
                double sum1=0.0, sum2=0.0;
                for (int k=fstart;k<fend;k++){ sum1 += buff[o-koff-k]*krn[k]; if(renorm) sum2+=krn[k]; }
                out[base+(long)o*stride] = renorm ? (sum2!=0.0 ? sum1/sum2 : 0.0) : sum1;
            }
        }
    }
    free(buff);
}

void smooth_uint8(Vol *V, const double fwhm[3]) {
    int nx=V->nx, ny=V->ny, nz=V->nz; long nvox=V->nvox;
    /* kernels (normalized to sum 1, as smooth_uint8 does) */
    double *kr[3]; int kl[3], koff[3];
    for (int d=0; d<3; d++) {
        int lim=(int)ceil(2*fwhm[d]);
        kr[d]=xmalloc(sizeof(double)*(2*lim+1));
        kl[d]=spm_smoothkern(fwhm[d], kr[d], 1);
        double s=0; for(int i=0;i<kl[d];i++) s+=kr[d][i];
        for(int i=0;i<kl[d];i++) kr[d][i]/=s;
        koff[d]=-lim;
    }
    sflt *a=xmalloc(sizeof(sflt)*nvox), *b=xmalloc(sizeof(sflt)*nvox);
    for (long k=0;k<nvox;k++) a[k]=V->u8[k];
    conv_axis(a,b,nx,ny,nz,0,kr[0],kl[0],koff[0],0); /* x, no renorm */
    conv_axis(b,a,nx,ny,nz,1,kr[1],kl[1],koff[1],0); /* y, no renorm */
    conv_axis(a,b,nx,ny,nz,2,kr[2],kl[2],koff[2],1); /* z, renorm */
    for (long k=0;k<nvox;k++) {
        double v=floor(b[k]+0.5);          /* RINT = floor(x+0.5) */
        if (v<0) v=0; if (v>255) v=255;
        V->u8[k]=(unsigned char)v;
    }
    free(a); free(b); for(int d=0;d<3;d++) free(kr[d]);
}

#ifdef SMOOTH_TEST
#include <stdio.h>
#include <assert.h>
int main(void){
    /* (a) kernel values vs golden */
    FILE *f=fopen("golden/smoothkern_golden.txt","r"); if(!f){fprintf(stderr,"no kern golden\n");return 2;}
    int t; double fwhm; double worst=0; int rows=0;
    double krn[64];
    while (fscanf(f,"%d %lf",&t,&fwhm)==2) {
        int lim=(int)ceil(2*fwhm), n=2*lim+1;
        double gv[17]; for(int i=0;i<17;i++) if(fscanf(f,"%lf",&gv[i])!=1) return 3;
        spm_smoothkern(fwhm,krn,t); (void)n;
        /* golden x range is -8..8 (17 vals); compare where |x|<=lim (kernel support) */
        for (int gi=0;gi<17;gi++){ int x=gi-8; if(abs(x)<=lim){ double d=fabs(krn[x+lim]-gv[gi]); if(d>worst)worst=d; } }
        rows++;
    }
    fclose(f);
    printf("smoothkern rows=%d  max|krn-golden|=%.3e\n",rows,worst);
    assert(worst<1e-12);

    /* (b) forced fwhm=[3 3 3] smoothing of template vs golden volume */
    Vol V; if(load_uint8("golden/avg152T1.nii",&V,0)){fprintf(stderr,"load fail\n");return 2;}
    double fw[3]={3,3,3}; smooth_uint8(&V,fw);
    f=fopen("golden/vf_smooth3.bin","rb"); unsigned char *g=malloc(V.nvox);
    long nn=fread(g,1,V.nvox,f); fclose(f); assert(nn==V.nvox);
    long diff=0,maxd=0; for(long k=0;k<V.nvox;k++){int dd=abs((int)V.u8[k]-(int)g[k]); if(dd){diff++; if(dd>maxd)maxd=dd;}}
    printf("smooth3 vol: nvox=%ld mismatches=%ld maxabs=%ld\n",V.nvox,diff,maxd);
    free(g); vol_free(&V);
    assert(diff==0);
    printf("M3 GATE PASS\n");
    return 0;
}
#endif
