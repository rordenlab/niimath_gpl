/* hist2.c — faithful port of SPM src/hist2.c (the spm_hist2 MEX core).
 * Float arithmetic for sampling/coords, double for M and H, matching SPM. */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "hist2.h"
#include "matrix.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* trilinear sample of uint8 volume of dims d at (x,y,z) [1-based, float] */
static float samp_vol(const int d[3], const unsigned char *f, float x, float y, float z) {
    int ix=(int)floorf(x), iy=(int)floorf(y), iz=(int)floorf(z);
    float dx1=x-ix, dx2=1.0f-dx1;
    float dy1=y-iy, dy2=1.0f-dy1;
    float dz1=z-iz, dz2=1.0f-dz1;
    const unsigned char *ff = f + ix-1 + d[0]*(iy-1 + (long)d[1]*(iz-1));
    float k222=ff[0],    k122=ff[1];
    float k212=ff[d[0]], k112=ff[d[0]+1];
    float vf;
    if (iz < d[2]-1) {
        ff += (long)d[0]*d[1];
        float k221=ff[0],    k121=ff[1];
        float k211=ff[d[0]], k111=ff[d[0]+1];
        vf = (((k222*dx2+k122*dx1)*dy2 + (k212*dx2+k112*dx1)*dy1))*dz2
           + (((k221*dx2+k121*dx1)*dy2 + (k211*dx2+k111*dx1)*dy1))*dz1;
    } else {
        vf = (((k222*dx2+k122*dx1)*dy2 + (k212*dx2+k112*dx1)*dy1));
    }
    return vf;
}

static const float ran[97] = {
0.656619f,0.891183f,0.488144f,0.992646f,0.373326f,0.531378f,0.181316f,0.501944f,0.422195f,
0.660427f,0.673653f,0.95733f,0.191866f,0.111216f,0.565054f,0.969166f,0.0237439f,0.870216f,
0.0268766f,0.519529f,0.192291f,0.715689f,0.250673f,0.933865f,0.137189f,0.521622f,0.895202f,
0.942387f,0.335083f,0.437364f,0.471156f,0.14931f,0.135864f,0.532498f,0.725789f,0.398703f,
0.358419f,0.285279f,0.868635f,0.626413f,0.241172f,0.978082f,0.640501f,0.229849f,0.681335f,
0.665823f,0.134718f,0.0224933f,0.262199f,0.116515f,0.0693182f,0.85293f,0.180331f,0.0324186f,
0.733926f,0.536517f,0.27603f,0.368458f,0.0128863f,0.889206f,0.866021f,0.254247f,0.569481f,
0.159265f,0.594364f,0.3311f,0.658613f,0.863634f,0.567623f,0.980481f,0.791832f,0.152594f,
0.833027f,0.191863f,0.638987f,0.669f,0.772088f,0.379818f,0.441585f,0.48306f,0.608106f,
0.175996f,0.00202556f,0.790224f,0.513609f,0.213229f,0.10345f,0.157337f,0.407515f,0.407757f,
0.0526927f,0.941815f,0.149972f,0.384374f,0.311059f,0.168534f,0.896648f };

void coreg_hist2(const Vol *VG, const Vol *VF, const double *P, int np,
                 double samp, double H[65536]) {
    /* M = inv(VF.mat) * spm_matrix(P) * VG.mat */
    double A[4][4], VFi[4][4], tmp[4][4], M[4][4];
    spm_matrix(P, np, A);
    if (mat44_inv(VF->mat, VFi)) { memset(H,0,65536*sizeof(double)); return; }
    mat44_mul(VFi, A, tmp);
    mat44_mul(tmp, VG->mat, M);

    int dg[3]={VG->nx,VG->ny,VG->nz}, df[3]={VF->nx,VF->ny,VF->nz};
    const unsigned char *g=VG->u8, *f=VF->u8;
    float s0[3], s[3];
    for (int i=0;i<3;i++){ s0[i]=(float)(samp/VG->vox[i]); s[i]=(dg[i]>1)?s0[i]:0.0f; }

    memset(H,0,65536*sizeof(double));
    int iran=0;
    for (float z=1.0f; z<=dg[2]-s[2]; z+=s0[2])
    for (float y=1.0f; y< dg[1]-s[1]; y+=s0[1])
    for (float x=1.0f; x< dg[0]-s[0]; x+=s0[0]) {
        float rx = x + ran[iran=(iran+1)%97]*s[0];
        float ry = y + ran[iran=(iran+1)%97]*s[1];
        float rz = z + ran[iran=(iran+1)%97]*s[2];
        float xp = (float)(M[0][0]*rx + M[0][1]*ry + M[0][2]*rz + M[0][3]);
        float yp = (float)(M[1][0]*rx + M[1][1]*ry + M[1][2]*rz + M[1][3]);
        float zp = (float)(M[2][0]*rx + M[2][1]*ry + M[2][2]*rz + M[2][3]);
        if (zp>=1.0f && zp<=df[2] && yp>=1.0f && yp<df[1] && xp>=1.0f && xp<df[0]) {
            float vf = samp_vol(df, f, xp,yp,zp);
            int ivf = (int)floorf(vf);
            int ivg = (int)floorf(samp_vol(dg, g, rx,ry,rz)+0.5f);
            H[ivf+ivg*256] += (1-(vf-ivf));
            if (ivf<255) H[ivf+1+ivg*256] += (vf-ivf);
        }
    }
}

/* Full partial-volume variant: distributes each sample over BOTH intensity axes
 * (the commented-out path in SPM's hist2.c). For the dither experiment only. */
void coreg_hist2_fullpv(const Vol *VG, const Vol *VF, const double *P, int np,
                        double samp, double H[65536]) {
    double A[4][4], VFi[4][4], tmp[4][4], M[4][4];
    spm_matrix(P, np, A);
    if (mat44_inv(VF->mat, VFi)) { memset(H,0,65536*sizeof(double)); return; }
    mat44_mul(VFi, A, tmp); mat44_mul(tmp, VG->mat, M);
    int dg[3]={VG->nx,VG->ny,VG->nz}, df[3]={VF->nx,VF->ny,VF->nz};
    const unsigned char *g=VG->u8, *f=VF->u8;
    float s0[3], s[3];
    for (int i=0;i<3;i++){ s0[i]=(float)(samp/VG->vox[i]); s[i]=(dg[i]>1)?s0[i]:0.0f; }
    memset(H,0,65536*sizeof(double));
    int iran=0;
    for (float z=1.0f; z<=dg[2]-s[2]; z+=s0[2])
    for (float y=1.0f; y< dg[1]-s[1]; y+=s0[1])
    for (float x=1.0f; x< dg[0]-s[0]; x+=s0[0]) {
        float rx = x + ran[iran=(iran+1)%97]*s[0];
        float ry = y + ran[iran=(iran+1)%97]*s[1];
        float rz = z + ran[iran=(iran+1)%97]*s[2];
        float xp = (float)(M[0][0]*rx + M[0][1]*ry + M[0][2]*rz + M[0][3]);
        float yp = (float)(M[1][0]*rx + M[1][1]*ry + M[1][2]*rz + M[1][3]);
        float zp = (float)(M[2][0]*rx + M[2][1]*ry + M[2][2]*rz + M[2][3]);
        if (zp>=1.0f && zp<=df[2] && yp>=1.0f && yp<df[1] && xp>=1.0f && xp<df[0]) {
            float vf = samp_vol(df, f, xp,yp,zp);
            float vg = samp_vol(dg, g, rx,ry,rz);
            int ivf=(int)floorf(vf), ivg=(int)floorf(vg);
            float df_=vf-ivf, dg_=vg-ivg;
            H[ivf+ivg*256] += (1-df_)*(1-dg_);
            if (ivf<255) H[ivf+1+ivg*256] += df_*(1-dg_);
            if (ivg<255){ H[ivf+(ivg+1)*256] += (1-df_)*dg_;
                          if (ivf<255) H[ivf+1+(ivg+1)*256] += df_*dg_; }
        }
    }
}

/* upper bound on samples; the large test pair is ~2M, so 100M is generous.
 * Also catches a degenerate -sep so tiny that samp/vox underflows the loop step. */
#define MAX_BASE_SAMPLES 100000000L

int base_samples_build(const Vol *VG, double samp, BaseSamples *bs) {
    bs->rx=bs->ry=bs->rz=NULL; bs->ivg=NULL; bs->n=0;   /* safe for base_samples_free on error */
    int dg[3]={VG->nx,VG->ny,VG->nz};
    const unsigned char *g=VG->u8;
    float s0[3], s[3];
    for (int i=0;i<3;i++){ s0[i]=(float)(samp/VG->vox[i]); s[i]=(dg[i]>1)?s0[i]:0.0f; }
    /* Pass 1: count exactly (no cap estimate -> no overflow). The bound also
     * terminates a degenerate step (s0 underflowed to 0 -> non-advancing loop). */
    long n=0;
    for (float z=1.0f; z<=dg[2]-s[2]; z+=s0[2])
    for (float y=1.0f; y< dg[1]-s[1]; y+=s0[1])
    for (float x=1.0f; x< dg[0]-s[0]; x+=s0[0])
        if (++n > MAX_BASE_SAMPLES) { fprintf(stderr,"coreg: degenerate sampling step (too many samples)\n"); return 1; }
    if (n==0) { fprintf(stderr,"coreg: no samples (check dims/sep)\n"); return 1; }
    bs->rx=malloc(sizeof(float)*n); bs->ry=malloc(sizeof(float)*n);
    bs->rz=malloc(sizeof(float)*n); bs->ivg=malloc(n);
    if(!bs->rx||!bs->ry||!bs->rz||!bs->ivg){ base_samples_free(bs); return 1; }
    /* Pass 2: fill. Identical loop + iran reset -> identical samples to the
     * original single-pass coreg_hist2 (validated bit-identical against it). */
    long k=0; int iran=0;
    for (float z=1.0f; z<=dg[2]-s[2]; z+=s0[2])
    for (float y=1.0f; y< dg[1]-s[1]; y+=s0[1])
    for (float x=1.0f; x< dg[0]-s[0]; x+=s0[0]) {
        float rx = x + ran[iran=(iran+1)%97]*s[0];
        float ry = y + ran[iran=(iran+1)%97]*s[1];
        float rz = z + ran[iran=(iran+1)%97]*s[2];
        bs->rx[k]=rx; bs->ry[k]=ry; bs->rz[k]=rz;
        bs->ivg[k]=(unsigned char)(int)floorf(samp_vol(dg,g,rx,ry,rz)+0.5f);
        k++;
    }
    bs->n=k;
    for(int r=0;r<4;r++) for(int c=0;c<4;c++) bs->mat[r][c]=VG->mat[r][c];
    return 0;
}

void base_samples_free(BaseSamples *bs){
    if(!bs) return;
    free(bs->rx); free(bs->ry); free(bs->rz); free(bs->ivg);
    bs->rx=bs->ry=bs->rz=NULL; bs->ivg=NULL; bs->n=0;
}

/* Accumulate samples [lo,hi) into Hdst. Shared by the serial path and each
 * OpenMP thread (which owns a private Hdst chunk-buffer). */
static void hist2_accum(const BaseSamples *bs, const unsigned char *f, const int df[3],
                        const double M[4][4], long lo, long hi, double *Hdst) {
    for (long i=lo;i<hi;i++) {
        float rx=bs->rx[i], ry=bs->ry[i], rz=bs->rz[i];
        float xp = (float)(M[0][0]*rx + M[0][1]*ry + M[0][2]*rz + M[0][3]);
        float yp = (float)(M[1][0]*rx + M[1][1]*ry + M[1][2]*rz + M[1][3]);
        float zp = (float)(M[2][0]*rx + M[2][1]*ry + M[2][2]*rz + M[2][3]);
        if (zp>=1.0f && zp<=df[2] && yp>=1.0f && yp<df[1] && xp>=1.0f && xp<df[0]) {
            float vf = samp_vol(df, f, xp,yp,zp);
            int ivf = (int)floorf(vf);
            int ivg = bs->ivg[i];
            Hdst[ivf+ivg*256] += (1-(vf-ivf));
            if (ivf<255) Hdst[ivf+1+ivg*256] += (vf-ivf);
        }
    }
}

/* Per-pass thread-buffer for the parallel build. On non-OpenMP builds nt is
 * forced to 1 and nothing is allocated (serial fallback, bit-identical to the
 * pre-OpenMP code). A malloc failure also degrades to serial rather than
 * failing the pass — the registration still runs single-threaded. */
int hist2_scratch_init(Hist2Scratch *hs) {
    hs->buf = NULL; hs->nt = 1;
#ifdef _OPENMP
    int nt = omp_get_max_threads();
    if (nt > 1) {
        hs->buf = malloc((size_t)nt*65536*sizeof(double));   /* zeroed per-slice each eval */
        if (hs->buf) hs->nt = nt;   /* else stay serial (nt=1, buf=NULL) */
    }
#endif
    return 0;
}
void hist2_scratch_free(Hist2Scratch *hs) {
    if (!hs) return;
    free(hs->buf); hs->buf = NULL; hs->nt = 1;
}

void coreg_hist2_cached(const BaseSamples *bs, const Vol *VF, const double *P, int np,
                        double H[65536], Hist2Scratch *hs) {
    double A[4][4], VFi[4][4], tmp[4][4], M[4][4];
    spm_matrix(P, np, A);
    if (mat44_inv(VF->mat, VFi)) { memset(H,0,65536*sizeof(double)); return; }
    mat44_mul(VFi, A, tmp); mat44_mul(tmp, bs->mat, M);
    int df[3]={VF->nx,VF->ny,VF->nz};
    const unsigned char *f=VF->u8;
    memset(H,0,65536*sizeof(double));
#ifdef _OPENMP
    /* Parallelize only when the sample count pays for thread setup (mirrors
     * allineate's >100000 guard). Each thread owns a contiguous static chunk
     * and a private 512 KB histogram slice (in the per-pass hs->buf — no
     * per-eval allocation), zeroing only its own slice, then we reduce in fixed
     * thread order.
     *
     * CRITICAL: chunk by the ACTUAL team size (omp_get_num_threads()), NOT the
     * requested hs->nt. num_threads(nt) is only an upper bound — the runtime may
     * launch fewer (OMP_THREAD_LIMIT, OMP_DYNAMIC, resource limits). Thread IDs
     * are then 0..actual-1; chunking by the requested nt would leave the chunks
     * for the unlaunched IDs unprocessed -> silently dropped samples. Chunking by
     * the actual team size makes IDs 0..nteam-1 tile [0,bs->n) exactly, and the
     * reduction runs over `got` (the actual size), so no stale/unwritten slice is
     * summed. Buffer safety: slice index t < nteam <= num_threads(nt) cap == the
     * hs->nt the buffer was sized for, so t is always in bounds.
     *
     * Deterministic for a given ACTUAL team size (fixed chunks + fixed reduction
     * order); float add order differs from serial and across team sizes, but stays
     * within the SPM golden tolerance (the histogram feeds a smoothed cost). */
    if (hs && hs->buf && hs->nt > 1 && bs->n > 100000) {
        int got = 1;
        #pragma omp parallel num_threads(hs->nt)
        {
            int nteam = omp_get_num_threads();      /* ACTUAL team size, may be < hs->nt */
            int t = omp_get_thread_num();
            if (t == 0) got = nteam;                /* publish once; region-end barrier makes it visible */
            double *Hl = hs->buf + (size_t)t*65536;
            memset(Hl, 0, 65536*sizeof(double));    /* each thread zeroes its OWN slice */
            long lo = bs->n * (long long)t     / nteam;   /* wide intermediate: bs->n*t can exceed 32-bit long (Windows LLP64) */
            long hi = bs->n * (long long)(t+1) / nteam;
            hist2_accum(bs, f, df, M, lo, hi, Hl);
        }
        for (int t=0;t<got;t++){ const double *Hl=hs->buf+(size_t)t*65536;
            for (int k=0;k<65536;k++) H[k]+=Hl[k]; }
        return;
    }
#endif
    (void)hs;
    hist2_accum(bs, f, df, M, 0, bs->n, H);
}

#ifdef HIST2_TEST
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "smooth.h"
static double check(const Vol *VG, const Vol *VF, double samp, const char *golden){
    double H[65536]; double P[6]={0,0,0,0,0,0};
    coreg_hist2(VG,VF,P,6,samp,H);
    FILE *fp=fopen(golden,"rb"); if(!fp){fprintf(stderr,"no %s\n",golden);exit(2);}
    double *G=malloc(65536*sizeof(double)); long n=fread(G,sizeof(double),65536,fp); fclose(fp);
    assert(n==65536);
    double worst=0, sumH=0, sumG=0;
    for(int i=0;i<65536;i++){ double d=fabs(H[i]-G[i]); if(d>worst)worst=d; sumH+=H[i]; sumG+=G[i]; }
    printf("%s: sumH=%.6f sumG=%.6f max|H-G|=%.3e\n",golden,sumH,sumG,worst);
    free(G); return worst;
}
int main(void){
    Vol VG,VF; double sep=2;
    if(load_uint8("golden/6_anat-T1w_acq-tfl.nii",&VG,0)||load_uint8("golden/avg152T1.nii",&VF,0)){
        fprintf(stderr,"load fail\n");return 2;}
    /* smooth with fwhmg = sqrt(max(sep^2-vox^2,0))/vox (== 0 here, identity) */
    double fg[3],ff[3];
    for(int i=0;i<3;i++){ double v=VG.vox[i]; fg[i]=sqrt(fmax(sep*sep-v*v,0))/v; v=VF.vox[i]; ff[i]=sqrt(fmax(sep*sep-v*v,0))/v; }
    smooth_uint8(&VG,fg); smooth_uint8(&VF,ff);
    double w4=check(&VG,&VF,4,"golden/H_x0_s4.bin");
    double w2=check(&VG,&VF,2,"golden/H_x0_s2.bin");
    /* cached path must be BIT-IDENTICAL to coreg_hist2 across several transforms */
    double Ptest[5][6]={{0,0,0,0,0,0},{3,-20,-5,0.02,-0.01,-0.04},{-8,4,2,-0.05,0.1,0.03},
                        {1,1,1,0,0,0},{0,0,0,0.2,-0.15,0.1}};
    double wcache=0;
    for(double samp=4; samp>=2; samp-=2){
        BaseSamples bs; assert(base_samples_build(&VG,samp,&bs)==0);
        for(int t=0;t<5;t++){
            double Ho[65536],Hc[65536];
            coreg_hist2(&VG,&VF,Ptest[t],6,samp,Ho);
            coreg_hist2_cached(&bs,&VF,Ptest[t],6,Hc,NULL);   /* NULL hs: serial, bit-identical */
            for(int i=0;i<65536;i++){ double d=fabs(Ho[i]-Hc[i]); if(d>wcache)wcache=d; }
        }
        base_samples_free(&bs);
    }
    printf("cached-vs-original max|H-H|=%.3e (must be 0)\n",wcache);
    vol_free(&VG); vol_free(&VF);
    assert(w4<1e-6 && w2<1e-6);
    assert(wcache==0.0);
    printf("M4 GATE PASS\n");
    return 0;
}
#endif
