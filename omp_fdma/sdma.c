#include <stdlib.h>
#include <math.h>

#define isZero fabs(Delta)<EPS
#define EPS 1e-20

void sdma_finalize(double*,double*,double*,double*);

int sdma(double *a, double *b, double *c, double *d, double *e, double *f,double *g,double *h,double *y, int dim)
{
    if (fabs(d[0])<EPS || fabs(d[0])<(fabs(e[0])+fabs(f[0])+fabs(g[0])))
    {
        return 1;
    }
    int malloc_size = dim*sizeof(double);
    double *alpha = (double*)malloc(malloc_size);
    double *beta = (double*)malloc(malloc_size);
    double *gamma = (double*)malloc(malloc_size);
    double *delta=(double*)malloc((dim+1)*sizeof(double));
    alpha[1] = e[0]/d[0];
    beta[1]=f[0]/d[0];
    gamma[1]=g[0]/d[0];
    delta[1]=h[0]/d[0];
    int i;
    int N=dim-1;
    int N_1=N-1;
    int N_2=N_1-1;
    int N_3=N_2-1;
    double Delta;

    Delta=d[1]-c[1]*alpha[1];
    if (isZero || fabs(d[1])<(fabs(c[1])+fabs(e[1])+fabs(f[1])+fabs(g[1])))
    {
        sdma_finalize(alpha,beta,gamma,delta);
        return 2;
    }
    alpha[2] = (e[1]-c[1]*beta[1])/Delta;
    beta[2]=(f[1]-c[1]*gamma[1])/Delta;
    gamma[2]=g[1]/Delta;
    delta[2]=(h[1]+c[1]*delta[1])/Delta;

    Delta=d[2]-c[2]*alpha[2]+b[2]*(alpha[1]*alpha[2]-beta[1]);
    if (isZero || fabs(d[2])<(fabs(b[2])+fabs(c[2])+fabs(e[2])+fabs(f[2])+fabs(g[2])))
    {
        sdma_finalize(alpha,beta,gamma,delta);
        return 3;
    }
    alpha[3]=(e[2]-c[2]*beta[2]-b[2]*(gamma[1]-alpha[1]*beta[2]))/Delta;
    beta[3]=(f[2]-c[2]*gamma[2]+b[2]*alpha[1]*gamma[2])/Delta;
    gamma[3]=g[2]/Delta;
    delta[3]=(h[2]+c[2]*delta[2]-b[2]*(alpha[1]*delta[2]+delta[1]))/Delta;

    for (i=3;i<N;i++)
    {
        Delta = d[i]-a[i]*(alpha[i-2]*(alpha[i-1]*alpha[i]-beta[i-1])-beta[i-2]*alpha[i]+gamma[i-2])+b[i]*(alpha[i-1]*alpha[i]-beta[i-1])-c[i]*alpha[i];
        if (isZero || fabs(d[i])<(fabs(a[i])+fabs(b[i])+fabs(c[i])+fabs(e[i])+fabs(f[i])+fabs(g[i])))
        {
            sdma_finalize(alpha,beta,gamma,delta);
            return i+1;
        }
        Delta=d[i]-a[i]*((alpha[i-1]*alpha[i]-beta[i-1])*alpha[i-2]-alpha[i]*beta[i-2]+gamma[i-2])+b[i]*(alpha[i-1]*alpha[i]-beta[i-1])-c[i]*alpha[i];
        if (isZero || fabs(d[i])<(fabs(a[i])+fabs(b[i])+fabs(c[i])+fabs(e[i])+fabs(f[i])+fabs(g[i])))
        {
            sdma_finalize(alpha,beta,gamma,delta);
            return i+1;
        }
        alpha[i+1]=(e[i]-a[i]*(alpha[i-2]*(alpha[i-1]*beta[i]-gamma[i-1])-beta[i]*beta[i-2])+b[i]*(alpha[i-1]*beta[i]-gamma[i-1])-c[i]*beta[i])/Delta;
        beta[i+1]=(f[i]+b[i]*alpha[i-1]*gamma[i]-c[i]*gamma[i]-a[i]*gamma[i]*(alpha[i-2]*alpha[i-1]-beta[i-2]))/Delta;
        gamma[i+1]=g[i]/Delta;
        delta[i+1]=(h[i]+c[i]*delta[i]+a[i]*(alpha[i-2]*(alpha[i-1]*delta[i]+delta[i-1])+delta[i-2]-beta[i-2]*delta[i])-b[i]*(alpha[i-1]*delta[i]+delta[i-1]))/Delta;
    }
    Delta=d[N]-a[N]*((alpha[N_1]*alpha[N]-beta[N_1])*alpha[N_2]-alpha[N]*beta[N_2]+gamma[N_2])+b[N]*(alpha[N_1]*alpha[N]-beta[N_1])-c[N]*alpha[N];
    delta[dim]=(h[N]+c[N]*delta[N]+a[N]*(alpha[N_2]*(alpha[N_1]*delta[N]+delta[N_1])+delta[N_2]-beta[N_2]*delta[N])-b[N]*(alpha[N_1]*delta[N]+delta[N_1]))/Delta;
    y[N]=delta[dim];
    y[N_1]=alpha[N]*y[N]+delta[N];
    y[N_2]=alpha[N_1]*y[N_1]-beta[N_1]*y[N]+delta[N_1];
    for (i=N_3;i>-1;i--)
    {
        y[i]=alpha[i+1]*y[i+1]-beta[i+1]*y[i+2]+gamma[i+1]*y[i+3]+delta[i+1];
    }
    sdma_finalize(alpha,beta,gamma,delta);
    return 0;
}

void sdma_finalize(double *a, double *b, double *c,double *d)
{
    free(a);
    free(b);
    free(c);
    free(d);
}
