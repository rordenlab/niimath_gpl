/* cost.c — faithful port of spm_coreg.m:optfun (histogram smoothing + cost). */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "xalloc.h"
#include "cost.h"
#include "hist2.h"
#include "smooth.h"

int cost_parse(const char *s, cost_fun_t *out) {
    if (!strcmp(s,"mi"))  {*out=COST_MI;  return 0;}
    if (!strcmp(s,"nmi")) {*out=COST_NMI; return 0;}
    if (!strcmp(s,"ecc")) {*out=COST_ECC; return 0;}
    if (!strcmp(s,"ncc")) {*out=COST_NCC; return 0;}
    /* ls/pearson: allineate's names for normalized cross-correlation (Pearson) —
     * the within-modality cost. Identical to ncc; pair with -fwhm 0 0 for the
     * pure direct-Pearson (no histogram smoothing). */
    if (!strcmp(s,"ls") || !strcmp(s,"pearson")) {*out=COST_NCC; return 0;}
    return 1;
}

static void hist_sizes(const double fwhm[2], int *Nr, int *Nc) {
    *Nc = 256 + 2*(int)ceil(2*fwhm[0]);   /* growth along g axis (dim2) */
    *Nr = 256 + 2*(int)ceil(2*fwhm[1]);   /* growth along f axis (dim1) */
}

int cost_scratch_init(CostScratch *cs, const double fwhm[2]) {
    cs->tmp=cs->out=cs->s1=cs->s2=NULL;   /* safe for cost_scratch_free on early return */
    /* guard the fixed krn1[64]/krn2[64] for direct (non-CLI) callers */
    for (int i=0;i<2;i++) if (!(isfinite(fwhm[i]) && fwhm[i]>=0) || 2*(int)ceil(2*fwhm[i])+1 > 64) return 1;
    hist_sizes(fwhm, &cs->Nr, &cs->Nc);
    /* normalized separable smoothing kernels — fwhm is constant for the run */
    cs->L1 = spm_smoothkern(fwhm[0], cs->krn1, 1);   /* row kernel (dim2 = g) */
    cs->L2 = spm_smoothkern(fwhm[1], cs->krn2, 1);   /* col kernel (dim1 = f) */
    double s;
    s=0; for(int i=0;i<cs->L1;i++) s+=cs->krn1[i]; for(int i=0;i<cs->L1;i++) cs->krn1[i]/=s;
    s=0; for(int i=0;i<cs->L2;i++) s+=cs->krn2[i]; for(int i=0;i<cs->L2;i++) cs->krn2[i]/=s;
    cs->tmp = xmalloc((size_t)256*cs->Nc*sizeof(double));   /* every element written each eval */
    cs->out = xmalloc((size_t)cs->Nr*cs->Nc*sizeof(double));
    cs->s1  = xmalloc((size_t)cs->Nc*sizeof(double));
    cs->s2  = xmalloc((size_t)cs->Nr*sizeof(double));
    return 0;
}
void cost_scratch_free(CostScratch *cs) {
    if(!cs) return;
    free(cs->tmp); free(cs->out); free(cs->s1); free(cs->s2);
    cs->tmp=cs->out=cs->s1=cs->s2=NULL;
}

/* Smooth H (256x256) by conv2 with the precomputed separable kernels in cs,
 * into cs->out (Nr x Nc). All elements are fully written, so cs buffers are
 * reused across evals without zeroing. */
static void smooth_hist_into(const double H[65536], CostScratch *cs) {
    int Nc=cs->Nc, Nr=cs->Nr, L1=cs->L1, L2=cs->L2;
    const double *krn1=cs->krn1, *krn2=cs->krn2;
    double *tmp=cs->tmp, *out=cs->out;
    /* conv2 full along dim2 (g): tmp is 256 x Nc (column-major: tmp[f + c*256]) */
    for (int f=0; f<256; f++)
        for (int n=0; n<Nc; n++) {
            double acc=0.0;
            int jlo = n-(L1-1); if(jlo<0) jlo=0;
            int jhi = n;        if(jhi>255) jhi=255;
            for (int j=jlo;j<=jhi;j++) acc += H[f + j*256]*krn1[n-j];
            tmp[f + n*256] = acc;
        }
    /* conv2 full along dim1 (f): out is Nr x Nc (column-major: out[m + c*Nr]) */
    for (int n=0; n<Nc; n++)
        for (int m=0; m<Nr; m++) {
            double acc=0.0;
            int ilo = m-(L2-1); if(ilo<0) ilo=0;
            int ihi = m;        if(ihi>255) ihi=255;
            for (int i=ilo;i<=ihi;i++) acc += tmp[i + n*256]*krn2[m-i];
            out[m + n*Nr] = acc;
        }
}

/* cost on cs->out (Nr x Nc, column-major H[m + n*Nr]); marginals in cs->s1/s2 */
static double cost_on(CostScratch *cs, cost_fun_t cf) {
    int Nr=cs->Nr, Nc=cs->Nc; double *H=cs->out, *s1=cs->s1, *s2=cs->s2;
    long NN=(long)Nr*Nc;
    double eps=DBL_EPSILON, sh=0.0;
    for (long k=0;k<NN;k++){ H[k]+=eps; sh+=H[k]; }
    for (long k=0;k<NN;k++) H[k]/=sh;
    for (int n=0;n<Nc;n++){ double a=0; for(int m=0;m<Nr;m++) a+=H[m+n*Nr]; s1[n]=a; }
    for (int m=0;m<Nr;m++){ double a=0; for(int n=0;n<Nc;n++) a+=H[m+n*Nr]; s2[m]=a; }

    double o=0.0;
    if (cf==COST_MI || cf==COST_ECC) {
        double mi=0.0;
        for (int n=0;n<Nc;n++) for(int m=0;m<Nr;m++){
            double h=H[m+n*Nr]; mi += h*log2(h/(s2[m]*s1[n]));
        }
        if (cf==COST_MI) o=-mi;
        else { double d=0; for(int n=0;n<Nc;n++) d+=s1[n]*log2(s1[n]);
                            for(int m=0;m<Nr;m++) d+=s2[m]*log2(s2[m]);
               o = -(-2*mi/d); }
    } else if (cf==COST_NMI) {
        double num=0, den=0;
        for(int n=0;n<Nc;n++) num+=s1[n]*log2(s1[n]);
        for(int m=0;m<Nr;m++) num+=s2[m]*log2(s2[m]);
        for(int n=0;n<Nc;n++) for(int m=0;m<Nr;m++){ double h=H[m+n*Nr]; den+=h*log2(h); }
        o = -(num/den);
    } else { /* NCC */
        double m1=0,m2=0;
        for(int m=0;m<Nr;m++) m1+=s2[m]*(m+1);
        for(int n=0;n<Nc;n++) m2+=s1[n]*(n+1);
        double sig1=0,sig2=0;
        for(int m=0;m<Nr;m++) sig1+=s2[m]*((m+1)-m1)*((m+1)-m1);
        for(int n=0;n<Nc;n++) sig2+=s1[n]*((n+1)-m2)*((n+1)-m2);
        sig1=sqrt(sig1); sig2=sqrt(sig2);
        double ncc=0;
        for(int n=0;n<Nc;n++) for(int m=0;m<Nr;m++)
            ncc += H[m+n*Nr]*((m+1)-m1)*((n+1)-m2);
        o = -(ncc/(sig1*sig2));
    }
    return o;
}

/* smooth + cost using a caller-owned scratch (no allocation) */
static double cost_from_H_scratch(const double H256[65536], cost_fun_t cf,
                                  const double fwhm[2], CostScratch *cs) {
    if (fwhm[0]==0.0 && fwhm[1]==0.0) {
        /* no histogram smoothing: cost runs directly on the 256x256 histogram
         * (cs sized Nr=Nc=256). Skips the conv + the length-1-kernel path;
         * bit-identical to the general case, whose 1-tap kernel is exactly 1.0. */
        memcpy(cs->out, H256, 65536*sizeof(double));
        return cost_on(cs, cf);
    }
    smooth_hist_into(H256, cs);
    return cost_on(cs, cf);
}

double coreg_cost_from_H(const double H256[65536], cost_fun_t cf, const double fwhm[2]) {
    CostScratch cs;
    if (cost_scratch_init(&cs, fwhm)) return INFINITY;   /* invalid fwhm: worst cost */
    double o = cost_from_H_scratch(H256, cf, fwhm, &cs);
    cost_scratch_free(&cs);
    return o;
}

long g_cost_evals = 0;   /* profiling: number of objective evaluations */

double coreg_cost(const Vol *VG, const Vol *VF, const double *P, int np,
                  double samp, cost_fun_t cf, const double fwhm[2]) {
    g_cost_evals++;
    double H[65536];
    coreg_hist2(VG, VF, P, np, samp, H);
    return coreg_cost_from_H(H, cf, fwhm);
}

double coreg_cost_cached(const BaseSamples *base, const Vol *VF, const double *P, int np,
                         cost_fun_t cf, const double fwhm[2], CostScratch *cs) {
    g_cost_evals++;
    double H[65536];
    coreg_hist2_cached(base, VF, P, np, H);
    return cost_from_H_scratch(H, cf, fwhm, cs);   /* hot path: no per-eval allocation */
}

#ifdef COST_TEST
#include <stdio.h>
#include <assert.h>
int main(void){
    Vol VG,VF; double sep=2, fwhm[2]={7,7};
    if(load_uint8("golden/6_anat-T1w_acq-tfl.nii",&VG,0)||load_uint8("golden/avg152T1.nii",&VF,0)){
        fprintf(stderr,"load fail\n");return 2;}
    double fg[3],ff[3];
    for(int i=0;i<3;i++){double v=VG.vox[i]; fg[i]=sqrt(fmax(sep*sep-v*v,0))/v; v=VF.vox[i]; ff[i]=sqrt(fmax(sep*sep-v*v,0))/v;}
    smooth_uint8(&VG,fg); smooth_uint8(&VF,ff);
    double P[6]={0,0,0,0,0,0}; double worst=0;
    const char *names[4]={"mi","ecc","nmi","ncc"};        /* golden file order */
    cost_fun_t order[4]={COST_MI,COST_ECC,COST_NMI,COST_NCC};
    for (int si=0; si<2; si++) {
        double samp = si==0?4:2;
        char fn[64]; snprintf(fn,sizeof fn,"golden/cost_x0_s%d.txt",(int)samp);
        FILE *f=fopen(fn,"r"); assert(f);
        char tag[16]; double gval;
        while (fscanf(f,"%15s %lf",tag,&gval)==2) {
            cost_fun_t cf=0; for(int i=0;i<4;i++) if(!strcmp(tag,names[i])) cf=order[i];
            double o=coreg_cost(&VG,&VF,P,6,samp,cf,fwhm);
            double d=fabs(o-gval); if(d>worst)worst=d;
            printf("samp=%d %-4s C=%.15g golden=%.15g d=%.3e\n",(int)samp,tag,o,gval,d);
        }
        fclose(f);
    }
    vol_free(&VG); vol_free(&VF);
    printf("M5 max|cost-golden|=%.3e\n",worst);
    assert(worst<1e-9);
    printf("M5 GATE PASS\n");
    return 0;
}
#endif
