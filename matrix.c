/* matrix.c — faithful port of SPM spm_matrix.m (forward) + 4x4 affine math. */
#include <math.h>
#include <string.h>
#include "matrix.h"

void mat44_mul(const double A[4][4], const double B[4][4], double C[4][4]) {
    double T[4][4];
    for (int i=0;i<4;i++)
        for (int j=0;j<4;j++) {
            double s=0.0;
            for (int k=0;k<4;k++) s += A[i][k]*B[k][j];
            T[i][j]=s;
        }
    memcpy(C,T,sizeof T);
}

void spm_matrix(const double *P, int n, double A[4][4]) {
    /* pad: q = [0 0 0 0 0 0 1 1 1 0 0 0] */
    double p[12] = {0,0,0,0,0,0,1,1,1,0,0,0};
    for (int i=0;i<n && i<12;i++) p[i]=P[i];

    double c4=cos(p[3]), s4=sin(p[3]);
    double c5=cos(p[4]), s5=sin(p[4]);
    double c6=cos(p[5]), s6=sin(p[5]);

    double T[4][4]  = {{1,0,0,p[0]},{0,1,0,p[1]},{0,0,1,p[2]},{0,0,0,1}};
    double R1[4][4] = {{1,0,0,0},{0,c4,s4,0},{0,-s4,c4,0},{0,0,0,1}};
    double R2[4][4] = {{c5,0,s5,0},{0,1,0,0},{-s5,0,c5,0},{0,0,0,1}};
    double R3[4][4] = {{c6,s6,0,0},{-s6,c6,0,0},{0,0,1,0},{0,0,0,1}};
    double Z[4][4]  = {{p[6],0,0,0},{0,p[7],0,0},{0,0,p[8],0},{0,0,0,1}};
    double S[4][4]  = {{1,p[9],p[10],0},{0,1,p[11],0},{0,0,1,0},{0,0,0,1}};

    double R[4][4], TR[4][4], TRZ[4][4];
    mat44_mul(R1,R2,R);    /* R = R1*R2 */
    mat44_mul(R,R3,R);     /* R = R1*R2*R3 */
    mat44_mul(T,R,TR);     /* T*R */
    mat44_mul(TR,Z,TRZ);   /* T*R*Z */
    mat44_mul(TRZ,S,A);    /* T*R*Z*S */
}

int mat44_inv(const double A[4][4], double Inv[4][4]) {
    double a[4][8];
    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) { a[i][j]=A[i][j]; a[i][4+j]=(i==j)?1.0:0.0; }
    }
    for (int col=0; col<4; col++) {
        int piv=col; double best=fabs(a[col][col]);
        for (int r=col+1;r<4;r++) if (fabs(a[r][col])>best){best=fabs(a[r][col]);piv=r;}
        if (best==0.0) return 1;
        if (piv!=col) for (int j=0;j<8;j++){double t=a[col][j];a[col][j]=a[piv][j];a[piv][j]=t;}
        double d=a[col][col];
        for (int j=0;j<8;j++) a[col][j]/=d;
        for (int r=0;r<4;r++) {
            if (r==col) continue;
            double f=a[r][col];
            for (int j=0;j<8;j++) a[r][j]-=f*a[col][j];
        }
    }
    for (int i=0;i<4;i++) for (int j=0;j<4;j++) Inv[i][j]=a[i][4+j];
    return 0;
}

#ifdef MATRIX_TEST
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
int main(int argc, char **argv) {
    const char *gf = argc>1 ? argv[1] : "golden/matrix_golden.txt";
    FILE *f=fopen(gf,"r");
    if(!f){fprintf(stderr,"no golden %s\n",gf);return 2;}
    int n; double worst=0.0, worst_inv=0.0; int cases=0;
    while (fscanf(f,"%d",&n)==1) {
        double P[12]; for(int i=0;i<n;i++) if(fscanf(f,"%lf",&P[i])!=1)return 3;
        double G[4][4];
        for(int i=0;i<4;i++)for(int j=0;j<4;j++) if(fscanf(f,"%lf",&G[i][j])!=1)return 3;
        double A[4][4]; spm_matrix(P,n,A);
        for(int i=0;i<4;i++)for(int j=0;j<4;j++){double d=fabs(A[i][j]-G[i][j]); if(d>worst)worst=d;}
        /* inv round-trip */
        double Ai[4][4],I[4][4]; assert(mat44_inv(A,Ai)==0); mat44_mul(A,Ai,I);
        for(int i=0;i<4;i++)for(int j=0;j<4;j++){double e=fabs(I[i][j]-(i==j?1.0:0.0)); if(e>worst_inv)worst_inv=e;}
        cases++;
    }
    fclose(f);
    printf("M1 cases=%d  max|spm_matrix - golden|=%.3e  max|A*inv(A)-I|=%.3e\n",cases,worst,worst_inv);
    assert(worst<1e-12); assert(worst_inv<1e-12);
    printf("M1 GATE PASS\n");
    return 0;
}
#endif
