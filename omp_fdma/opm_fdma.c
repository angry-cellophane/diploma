#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

extern int sdma(double*,double*,double*,double*,double*,double*,double*,double*,double*,int);

extern int fdma(double *a, double *b, double *c, double *d, double *e, double *f,double *y, int dim);

int delta(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,int);

int union_sdma(double*A,double*B,double*C,double*D,double*E,double*F,double*a,double*b,double*c,double*d,double*e,double*f,double*g,double*h,double*u,double*v,double*p,double*q,double*w,double *y,int*,int,int);



int omp_fdma(double*a,double*b,double*c,double*d,double*e,double*f,double**y,int dim)
{
	int err=0;
	int *pr_el;
	int *disp;
	
	double *u_all;
	double *v_all;
	double *p_all;
	double *q_all;
	double *w_all;
	
	double *Y;
	
	int *offsets;
	int *length;
	
	#pragma omp parallel shared(a,b,c,d,e,f,u_all,v_all,p_all,q_all,w_all,Y,offsets,length,pr_el,disp) 
	
	{
		int pr_count = omp_get_num_threads();
		int thread_num = omp_get_thread_num();
		
		int i;
		
		#pragma omp single
		{
			pr_el = (int*) calloc(pr_count,sizeof(int));
			int row_counts = (dim+pr_count-2) / pr_count+1;
			for (i=0;i<pr_count;i++)
			{
				pr_el[i]=row_counts;
			}
			int rest_rows = (dim+pr_count-2) % pr_count;
			i=0;
			while (rest_rows)
			{
				pr_el[i]++;
				i++;
				rest_rows--;
			}
			
			disp=(int*)malloc(pr_count*sizeof(int));
			disp[0]=0;
			for (i=1;i<pr_count;i++)
			{
				disp[i]=pr_el[i-1]+disp[i-1]-2;
			}
		}	
		
		int lc_pr_count = pr_el[thread_num];
		
		double *_a = &a[disp[thread_num]];
		double *_b = &b[disp[thread_num]];
		double *_c = &c[disp[thread_num]];
		double *_d = &d[disp[thread_num]];
		double *_e = &e[disp[thread_num]];
		double *_f = &f[disp[thread_num]];
		
		double *u=calloc(lc_pr_count,sizeof(double));
		u[0]=1;
		u[1]=0;
		u[lc_pr_count-2]=0;
		u[lc_pr_count-1]=0;
		double *v=calloc(lc_pr_count,sizeof(double));
		v[0]=0;
		v[1]=1.0;
		v[lc_pr_count-2]=0;
		v[lc_pr_count-1]=0;
		double *p=calloc(lc_pr_count,sizeof(double));
		p[0]=0;
		p[1]=0;
		p[lc_pr_count-2]=1;
		p[lc_pr_count-1]=0;
		double *q=calloc(lc_pr_count,sizeof(double));
		q[0]=0;
		q[1]=0;
		q[lc_pr_count-2]=0;
		q[lc_pr_count-1]=1;
		double *w=calloc(lc_pr_count,sizeof(double));
		w[0]=0;
		w[1]=0;
		w[lc_pr_count-2]=0;
		w[lc_pr_count-1]=0;		
		
		err = delta(&_a[2],&_b[2],&_c[2],&_d[2],&_e[2],&_f[2],&u[2],&v[2],&p[2],&q[2],&w[2],lc_pr_count-4);		
		
		
		#pragma omp single
		{
			u_all = (double*)calloc(dim,sizeof(double));
			v_all = (double*)calloc(dim,sizeof(double));
			p_all = (double*)calloc(dim,sizeof(double));
			q_all = (double*)calloc(dim,sizeof(double));
			w_all = (double*)calloc(dim,sizeof(double));
		}
		
		for (i=0;i<lc_pr_count;i++)
		{
			u_all[disp[thread_num] + i] = u[i];
			v_all[disp[thread_num] + i] = v[i];
			p_all[disp[thread_num] + i] = p[i];
			q_all[disp[thread_num] + i] = q[i];
			w_all[disp[thread_num] + i] = w[i];
		} 
		
		#pragma omp single
		{
			int size = (pr_count<<1)+2;
			double *A = (double*)malloc(size*sizeof(double));
			double *B = (double*)malloc(size*sizeof(double));
			double *C = (double*)malloc(size*sizeof(double));
			double *D = (double*)malloc(size*sizeof(double));
			double *E = (double*)malloc(size*sizeof(double));
			double *F = (double*)malloc(size*sizeof(double));
			double *G = (double*)malloc(size*sizeof(double));
			double *H = (double*)malloc(size*sizeof(double));
			Y = (double*)malloc(dim*sizeof(double));
			err = union_sdma(a,b,c,d,e,f,A,B,C,D,E,F,G,H,u_all,v_all,p_all,q_all,w_all,Y,pr_el,size,pr_count);

			free(A);
			free(B);
			free(C);
			free(D);
			free(E);
			free(F);
			free(G);
			free(H);
			
			offsets=(int*)calloc(pr_count,sizeof(int));
			length=(int*)calloc(pr_count,sizeof(int));
			offsets[0]=0;
			length[0]=pr_el[0];
			for (i=1;i<pr_count;i++)
			{
				offsets[i]=offsets[i-1]+pr_el[i-1]-2;
				length[i]=pr_el[i];
			}
		}
		
		double *_y=calloc(lc_pr_count,sizeof(double));
		
		for (i=0;i<lc_pr_count;i++)
		{
			_y[i] = Y[offsets[thread_num]+i];
		}

		int k_m = lc_pr_count-2;
		int k_m_1 = lc_pr_count-1;
		for (i=2;i<k_m;i++)
		{
			_y[i]=_y[0]*u[i]+_y[1]*v[i]+_y[k_m]*p[i]+_y[k_m_1]*q[i]+w[i];
		}

		#pragma omp single
		{
			(*y)=(double*)calloc(dim,sizeof(double));
		}
		
		for (i=0;i<lc_pr_count;i++)
		{
			*( (*y) + disp[thread_num]+i) = _y[i];
		}
		
		#pragma omp master
		{
			free(Y);
			free(disp);
			free(length);
			free(offsets);
			free(pr_el);
		}
		free(u);
		free(v);
		free(p);
		free(q);
		free(w);
	}
	return err;
}

int delta(double *a,double *b,double *c,double *d,double *e,double *f,double *u,double *v,double *p,double *q,double *w,int dim)
{

    int err = fdma(a,b,c,d,e,f,w,dim);

    double *help_array=(double*)malloc(dim*sizeof(double));
    int i;
    for (i=0;i<dim;++i)
    {
        help_array[i]=f[i];
        f[i]=0;
    }

    f[0]=-a[0];
    err = fdma(a,b,c,d,e,f,u,dim);
    if (err!=0)
    {
        return 1;
    }
    f[0]=0;

    f[0] =b[0];
    f[1] =-a[1];
    err = fdma(a,b,c,d,e,f,v,dim);
    if (err!=0)
    {
        return 2;
    }
    f[0]=0;
    f[1]=0;

    int last = dim-1;
    int prelast=last-1;

    f[prelast]=-e[prelast];
    f[last]=d[last];
    err = fdma(a,b,c,d,e,f,p,dim);
    if (err!=0)
    {
        return 3;
    }
    f[prelast]=0;
    f[last]=0;

    f[last]=-e[last];
    err = fdma(a,b,c,d,e,f,q,dim);
    if (err!=0)
    {
        return 4;
    }
    f[last]=0;

    for (i=0;i<dim;i++)
    {
        f[i]=help_array[i];
    }
    free(help_array);
    return 0;
}



int union_sdma(
        double *A,
        double *B,
        double *C,
        double *D,
        double *E,
        double *F,
        double *a,
        double *b,
        double *c,
        double *d,
        double *e,
        double *f,
        double *g,
        double *h,
        double* u,
        double* v,
        double* p,
        double* q,
        double* w,
        double *y,
        int* array_m,
        int dim,
        int pr_count)
{
    int node_size = (pr_count<<1)+2;
    double *node_a = (double*) calloc(node_size,sizeof(double));
    double *node_b = (double*) calloc(node_size,sizeof(double));
    double *node_c = (double*) calloc(node_size,sizeof(double));
    double *node_d = (double*) calloc(node_size,sizeof(double));
    double *node_e = (double*) calloc(node_size,sizeof(double));
    double *node_f = (double*) calloc(node_size,sizeof(double));
    double *node_g = (double*) calloc(node_size,sizeof(double));
    double *node_h = (double*) calloc(node_size,sizeof(double));
    double *node_y = (double*) calloc(node_size,sizeof(double));

    //    row #0
    a[0]=0.0;
    node_a[0]=0.0;
    b[0]=0.0;
    node_b[0]=0.0;
    c[0]=0.0;
    node_c[0]=0.0;
    d[0]=C[0]+E[0]*u[2];
    node_d[0]=d[0];
    e[0]=D[0]-E[0]*v[2];
    node_e[0]=e[0];
    f[0]=E[0]*p[2];
    node_f[0]=f[0];
    g[0]=-E[0]*q[2];
    node_g[0]=g[0];
    h[0]=F[0]-E[0]*w[2];
    node_h[0]=h[0];

    //    row #1
    a[1]=0.0;
    b[1]=0.0;
    c[1]=B[1]+D[1]*u[2]-E[1]*u[3];
    d[1]=C[1]-D[1]*v[2]+E[1]*v[3];
    e[1]=D[1]*p[2]-E[1]*p[3];
    f[1]=E[1]*q[3]-D[1]*q[2];
    g[1]=0;
    h[1]=F[1]+D[1]*w[2]-E[1]*w[3];

    node_a[1]=a[1];
    node_b[1]=b[1];
    node_c[1]=c[1];
    node_d[1]=d[1];
    node_e[1]=e[1];
    node_f[1]=f[1];
    node_g[1]=g[1];
    node_h[1]=h[1];

    //    row #k,k+1
    int k=0;
    int i=0;
    int i_1;
    int last=dim-1;
    int prelast=last-1;
    int pr_countX2 = (pr_count<<1);
    for (i=2;i<pr_countX2;i+=2)
    {
//        printf("array_m[%d]=%d\n",(i>>1),array_m[i>>1]);
        k += array_m[((i-2)>>1)]-2;
//        printf("k=%d\n",k);
        a[i]=0;
        b[i]=A[k]*u[k-2]-B[k]*u[k-1];
        c[i]=B[k]*v[k-1]-A[k]*v[k-2];
        d[i]=A[k]*p[k-2]+C[k]+E[k]*u[k+2]-B[k]*p[k-1];
        e[i]=B[k]*q[k-1]-A[k]*q[k-2]+D[k]-E[k]*v[k+2];
        f[i]=E[k]*p[k+2];
        g[i]=-E[k]*q[k+2];
        h[i]=F[k]-A[k]*w[k-2]+B[k]*w[k-1]-E[k]*w[k+2];

        node_a[i]=a[i];
        node_b[i]=b[i];
        node_c[i]=c[i];
        node_d[i]=d[i];
        node_e[i]=e[i];
        node_f[i]=f[i];
        node_g[i]=g[i];
        node_h[i]=h[i];

//        k+1
        i_1=i+1;
        a[i_1]=-A[k+1]*u[k-1];
        b[i_1]=A[k+1]*v[k-1];
        c[i_1]=B[k+1]+D[k+1]*u[k+2]-A[k+1]*p[k-1]-E[k+1]*u[k+3];
        d[i_1]=A[k+1]*q[k-1]+C[k+1]+E[k+1]*v[k+3]-D[k+1]*v[k+2];
        e[i_1]=D[k+1]*p[k+2]-E[k+1]*p[k+3];
        f[i_1]=E[k+1]*q[k+3]-D[k+1]*q[k+2];
        g[i_1]=0;
        h[i_1]=F[k+1]-A[k+1]*w[k-1]+D[k+1]*w[k+2]-E[k+1]*w[k+3];

        node_a[i_1]=a[i_1];
        node_b[i_1]=b[i_1];
        node_c[i_1]=c[i_1];
        node_d[i_1]=d[i_1];
        node_e[i_1]=e[i_1];
        node_f[i_1]=f[i_1];
        node_g[i_1]=g[i_1];
        node_h[i_1]=h[i_1];
    }

    //    row #N-1
    k += array_m[((prelast-2)>>1)]-2;
    a[prelast]=0.0;
    b[prelast]=A[k]*u[k-2]-B[k]*u[k-1];
//    c[prelast]=B[k]*v[k]-A[k]*v[k-2]; // было
    c[prelast]=B[k]*v[k-1]-A[k]*v[k-2];
    d[prelast]=A[k]*p[k-2]-B[k]*p[k-1]+C[k];
    e[prelast]=B[k]*q[k-1]+D[k]-A[k]*q[k-2];
    f[prelast]=0.0;
    g[prelast]=0.0;
    h[prelast]=F[k]+B[k]*w[k-1]-A[k]*w[k-2];

    node_a[prelast]=a[prelast];
    node_b[prelast]=b[prelast];
    node_c[prelast]=c[prelast];
    node_d[prelast]=d[prelast];
    node_e[prelast]=e[prelast];
    node_f[prelast]=f[prelast];
    node_g[prelast]=g[prelast];
    node_h[prelast]=h[prelast];

    //    row #N
    k++;
//    printf("last=%d\nk[N]=%d\n",last,k);
    a[last]=-A[k]*u[k-2];
    b[last]=A[k]*v[k-2];
    c[last]=B[k]-A[k]*p[k-2];
    d[last]=C[k]+A[k]*q[k-2];
    e[last]=0.0;
    f[last]=0.0;
    g[last]=0.0;
    h[last]=F[k]-A[k]*w[k-2];

    node_a[last]=a[last];
    node_b[last]=b[last];
    node_c[last]=c[last];
    node_d[last]=d[last];
    node_e[last]=e[last];
    node_f[last]=f[last];
    node_g[last]=g[last];
    node_h[last]=h[last];

    int err = sdma(node_a,node_b,node_c,node_d,node_e,node_f,node_g,node_h,node_y,node_size);

    if (err!=0)
    {
        return err;
    }

    y[0]=node_y[0];
    y[1]=node_y[1];
    k=0;
    i=0;
    last=dim-1;
    prelast=last-1;
    pr_countX2 = (pr_count<<1);
    for (i=2;i<pr_countX2;i+=2)
    {

        k += array_m[((i-2)>>1)]-2;
//        printf("k=%d\n",k);
        y[k]=node_y[i];
        y[k+1]=node_y[i+1];

    }

    //    row #N-1
    k += array_m[((prelast-2)>>1)]-2;
    y[k]=node_y[prelast];
    y[k+1]=node_y[last];

    free(node_a);
    free(node_b);
    free(node_c);
    free(node_d);
    free(node_e);
    free(node_f);
    free(node_g);
    free(node_h);
    free(node_y);

    return 0;
}
