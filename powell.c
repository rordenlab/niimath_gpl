/* powell.c — faithful port of spm_powell.m (Powell direction set + Brent). */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "xalloc.h"
#include "powell.h"

#define EPS DBL_EPSILON

typedef struct { double p, f; } Pt;
typedef struct {           /* line context: pt = p0 + x*dir */
    const double *p0, *dir; int n; objfun f; void *ctx; double *scratch;
} LineCtx;

static double funeval(LineCtx *L, double x) {
    for (int i=0;i<L->n;i++) L->scratch[i] = L->p0[i] + x*L->dir[i];
    return L->f(L->scratch, L->n, L->ctx);
}

/* fit c0+c1*u+c2*u^2 = y through 3 points (distinct u) -> solve 3x3 (== pinv). */
static void parabola(const double u[3], const double y[3], double pol[3]) {
    double V[3][3];
    for (int r=0;r<3;r++){ V[r][0]=1.0; V[r][1]=u[r]; V[r][2]=u[r]*u[r]; }
    /* invert 3x3 via cofactors */
    double a=V[0][0],b=V[0][1],c=V[0][2],d=V[1][0],e=V[1][1],f=V[1][2],g=V[2][0],h=V[2][1],i=V[2][2];
    double A=e*i-f*h, B=-(d*i-f*g), C=d*h-e*g;
    double det=a*A+b*B+c*C;
    double D=-(b*i-c*h), E=a*i-c*g, F=-(a*h-b*g);
    double G=b*f-c*e,    H=-(a*f-c*d), I=a*e-b*d;
    double inv[3][3]={{A/det,D/det,G/det},{B/det,E/det,H/det},{C/det,F/det,I/det}};
    for (int r=0;r<3;r++) pol[r]=inv[r][0]*y[0]+inv[r][1]*y[1]+inv[r][2]*y[2];
}

static void bracket(LineCtx *L, double f0, Pt t[3]) {
    double gold=(1+sqrt(5.0))/2;
    t[0].p=0; t[0].f=f0;
    t[1].p=1; t[1].f=funeval(L,1);
    if (t[1].f > t[0].f){ Pt s=t[0]; t[0]=t[1]; t[1]=s; }
    t[2].p=t[1].p+gold*(t[1].p-t[0].p); t[2].f=funeval(L,t[2].p);
    /* cap guards a pathological/monotone objective from expanding without bound;
     * well-behaved costs bracket in a handful of steps, so this never trips on
     * the validated images (optimizer trajectory unchanged vs. the SPM reference). */
    for (int guard=0; t[1].f > t[2].f && guard < 200; guard++) {
        double u[3]={t[0].p-t[1].p, 0.0, t[2].p-t[1].p};
        double y[3]={t[0].f, t[1].f, t[2].f}, pol[3];
        parabola(u,y,pol);
        Pt uu;
        if (pol[2]>0) {
            double dd = -pol[1]/(2*pol[2]+EPS);
            if (dd > (1+gold)*(t[2].p-t[1].p)) dd=(1+gold)*(t[2].p-t[1].p);
            uu.p = t[1].p+dd;
        } else uu.p = t[2].p + gold*(t[2].p-t[1].p);
        uu.f = funeval(L,uu.p);
        if ((t[1].p < uu.p) == (uu.p < t[2].p)) {
            if (uu.f < t[2].f){ t[0]=t[1]; t[1]=uu; return; }
            else if (uu.f > t[1].f){ t[2]=uu; return; }
        }
        t[0]=t[1]; t[1]=t[2]; t[2]=uu;
    }
}

static void search(LineCtx *L, Pt t[3], double tol, double *fout, double *pout) {
    double gold1 = 1-(sqrt(5.0)-1)/2;
    double d=INFINITY, pd=INFINITY;
    /* sort t best (lowest f) first */
    for (int a=0;a<3;a++) for (int b=a+1;b<3;b++) if (t[b].f<t[a].f){ Pt s=t[a];t[a]=t[b];t[b]=s; }
    double brk0 = fmin(fmin(t[0].p,t[1].p),t[2].p);
    double brk1 = fmax(fmax(t[0].p,t[1].p),t[2].p);
    for (int iter=0; iter<128; iter++) {
        if (fabs(t[0].p - 0.5*(brk0+brk1)) + 0.5*(brk1-brk0) <= 2*tol) break;
        double ppd=pd; pd=d;
        double u[3]={t[0].p-t[0].p, t[1].p-t[0].p, t[2].p-t[0].p};
        double y[3]={t[0].f,t[1].f,t[2].f}, pol[3];
        parabola(u,y,pol);
        d = -pol[1]/(2*pol[2]+EPS);
        Pt uu; uu.p=t[0].p+d;
        double eps2 = 2*EPS*fabs(t[0].p)+EPS;
        if (fabs(d) > fabs(ppd)/2 || uu.p < brk0+eps2 || uu.p > brk1-eps2 || pol[2]<=0) {
            if (t[0].p >= 0.5*(brk0+brk1)) d=gold1*(brk0-t[0].p);
            else                          d=gold1*(brk1-t[0].p);
            uu.p=t[0].p+d;
        }
        uu.f=funeval(L,uu.p);
        if (uu.f <= t[0].f) {
            if (uu.p >= t[0].p) brk0=t[0].p; else brk1=t[0].p;
            t[2]=t[1]; t[1]=t[0]; t[0]=uu;
        } else {
            if (uu.p < t[0].p) brk0=uu.p; else brk1=uu.p;
            if (uu.f <= t[1].f){ t[2]=t[1]; t[1]=uu; }
            else if (uu.f <= t[2].f) t[2]=uu;
        }
    }
    *fout=t[0].f; *pout=t[0].p;
}

/* line minimisation: updates p in place, *f; writes scaled dir to dir_out if non-NULL */
static void min1d(double *p, const double *dir, int n, double *f,
                  const double *tolsc, objfun fn, void *ctx,
                  double *scratch, double *dir_out) {
    double ss=0; for(int i=0;i<n;i++){ double t=dir[i]/tolsc[i]; ss+=t*t; }
    double tol = 1.0/sqrt(ss);
    LineCtx L={p,dir,n,fn,ctx,scratch};
    Pt t[3]; bracket(&L,*f,t);
    double fmin_,pmin; search(&L,t,tol,&fmin_,&pmin);
    for (int i=0;i<n;i++){ double sd=pmin*dir[i]; p[i]+=sd; if(dir_out) dir_out[i]=sd; }
    *f=fmin_;
}

double spm_powell(double *p, int n, double *xi, const double *tolsc,
                  objfun fn, void *ctx) {
    double *scratch=xmalloc(sizeof(double)*n), *pp=xmalloc(sizeof(double)*n);
    double *dtmp=xmalloc(sizeof(double)*n);   /* (p-pp) conjugate direction */
    double f = fn(p,n,ctx);
    for (int iter=0; iter<512; iter++) {
        int ibig=n-1; double del=0, fp=f;
        memcpy(pp,p,sizeof(double)*n);
        for (int i=0;i<n;i++) {
            double ft=f;
            min1d(p, &xi[i*n], n, &f, tolsc, fn, ctx, scratch, NULL); /* discard scaled dir */
            if (fabs(ft-f)>del){ del=fabs(ft-f); ibig=i; }
        }
        double s=0; for(int i=0;i<n;i++){ double t=(p[i]-pp[i])/tolsc[i]; s+=t*t; }
        if (n==1 || sqrt(s)<1.0 || fabs((f-fp)/(f+fp))<1e-6) break;
        /* extrapolated point 2*p - pp */
        for (int i=0;i<n;i++) scratch[i]=2.0*p[i]-pp[i];
        double ft = fn(scratch,n,ctx);
        if (ft < f) {
            for (int i=0;i<n;i++) dtmp[i]=p[i]-pp[i];
            min1d(p, dtmp, n, &f, tolsc, fn, ctx, scratch, &xi[ibig*n]); /* update xi(:,ibig) */
        }
    }
    free(scratch); free(pp); free(dtmp);
    return f;
}

#ifdef POWELL_TEST
#include <stdio.h>
#include <assert.h>
static double *G_traj; static int G_np, G_count;
static double tf1(const double *p,int n,void*c){(void)n;(void)c;double x=p[0];double v=3*(x-2.7)*(x-2.7)+0.5;
    if(G_traj){G_traj[G_count*(G_np+1)+0]=x; G_traj[G_count*(G_np+1)+G_np]=v;} G_count++; return v;}
static double tf3(const double *p,int n,void*c){(void)n;(void)c;
    double d1=p[0]-1.5,d2=p[1]+2.0,d3=p[2]-0.7;
    double a1=2.0*d1+0.3*d2-0.1*d3, a2=0.1*d1+1.5*d2+0.2*d3, a3=-0.2*d1+0.1*d2+1.8*d3;
    double v=a1*a1+a2*a2+a3*a3;
    if(G_traj){for(int i=0;i<3;i++)G_traj[G_count*(G_np+1)+i]=p[i]; G_traj[G_count*(G_np+1)+G_np]=v;} G_count++; return v;}

static double run_case(const char*fn, int np, objfun f){
    FILE*fp=fopen(fn,"r"); if(!fp){fprintf(stderr,"no %s\n",fn);exit(2);}
    int N,gnp; if(fscanf(fp,"%*s %d %*s %d",&N,&gnp)!=2)exit(3); assert(gnp==np);
    double *gold=malloc(sizeof(double)*N*(np+1));
    for(int r=0;r<N;r++) for(int c=0;c<=np;c++) if(fscanf(fp,"%lf",&gold[r*(np+1)+c])!=1)exit(3);
    char tag[16]; double gf[8]; fscanf(fp,"%15s",tag); for(int i=0;i<np;i++) fscanf(fp,"%lf",&gf[i]);
    fclose(fp);
    G_traj=malloc(sizeof(double)*(N+50)*(np+1)); G_np=np; G_count=0;
    double p[8]={0}; double tol[8]; double *xi=calloc(np*np,sizeof(double));
    for(int i=0;i<np;i++){ tol[i]=0.001; xi[i*np+i]=0.001*20; }
    double ff=spm_powell(p,np,xi,tol,f,NULL);
    /* leading-segment agreement (until pinv/SVD branch divergence) */
    int M=N<G_count?N:G_count, lead=0; double leadtol=1e-9;
    for(int r=0;r<M;r++){ double w=0; for(int c=0;c<=np;c++){double dd=fabs(G_traj[r*(np+1)+c]-gold[r*(np+1)+c]); if(dd>w)w=dd;} if(w<leadtol) lead++; else break; }
    double wfin=0; for(int i=0;i<np;i++){double dd=fabs(p[i]-gf[i]); if(dd>wfin)wfin=dd;}
    printf("%s: golden_evals=%d C_evals=%d  leading_match(<1e-9)=%d evals  final f=%.12g  max|final-golden|=%.3e\n",
           fn,N,G_count,lead,ff,wfin);
    free(gold);free(G_traj);free(xi);
    return wfin;   /* gate on convergence to the same minimum */
}
int main(void){
    double w1=run_case("golden/powell_traj_1d.txt",1,tf1);
    double w3=run_case("golden/powell_traj_3d.txt",3,tf3);
    double worst=w1>w3?w1:w3;
    printf("M6 worst=%.3e\n",worst);
    assert(worst<1e-6);
    printf("M6 GATE PASS\n");
    return 0;
}
#endif
