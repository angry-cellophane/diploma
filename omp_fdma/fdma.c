#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define isZero fabs(delta)<EPS
#define EPS 1e-20

void fdma_finalize(double*,double*,double*);

int fdma(double *a, double *b, double *c, double *d, double *e, double *f,double *y, int dim)
{	
    if (fabs(c[0])< EPS || fabs(c[0])<(fabs(d[0])+fabs(e[0])))
    {
        return 1;
    }
    int malloc_size = dim*sizeof(double);
    double *alpha = (double*)malloc(malloc_size);
    double *beta = (double*)malloc(malloc_size);
    double *gamma = (double*)malloc((dim+1)*sizeof(double));
    alpha[1] = d[0] / c[0];
    beta[1] = e[0]/c[0];
    gamma[1] = f[0]/c[0];
    double delta = c[1] - b[1]*alpha[1];
    if (isZero || fabs(c[1])<(fabs(d[1])+fabs(e[1])+fabs(b[1])))
    {
        fdma_finalize(alpha,beta,gamma);
        return 2;
    }
    alpha[2]=(d[1]-b[1]*beta[1])/delta;
    beta[2]=e[1]/delta;
    gamma[2] = (f[1]+b[1]*gamma[1])/delta;
    int i;
    int last = dim-1;
    for (i=2;i<last;i++)
    {
        delta = c[i]-a[i]*beta[i-1]+alpha[i]*(a[i]*alpha[i-1]-b[i]);
        if (isZero || fabs(c[i])<(fabs(a[i])+fabs(b[i]+fabs(d[i])+fabs(e[i]))))
        {
            fdma_finalize(alpha,beta,gamma);
            return i+1;
        }
        alpha[i+1] = (d[i]+beta[i]*(a[i]*alpha[i-1]-b[i]))/delta;
        beta[i+1] = e[i]/delta;
        gamma[i+1] = (f[i]-a[i]*gamma[i-1]-gamma[i]*(a[i]*alpha[i-1]-b[i]))/delta;
    }
    int prelast = last-1;
    delta = c[last] - a[last] * beta[prelast] + alpha[last] * (a[last]*alpha[prelast]-b[last]);
    gamma[last+1] = (f[last]-a[last]*gamma[prelast]-gamma[last]*(a[last]*alpha[prelast]-b[last]))/delta;
    y[last] = gamma[last+1];
    y[prelast] = alpha[last]*y[last]+gamma[last];
    for (i=last-2;i>-1;i--)
    {
        y[i]=alpha[i+1]*y[i+1]-beta[i+1]*y[i+2]+gamma[i+1];
    }
    fdma_finalize(alpha,beta,gamma);
    return 0;
}

void fdma_finalize(double *a, double *b, double *c)
{
    free(a);
    free(b);
    free(c);
}
