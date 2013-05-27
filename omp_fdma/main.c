#include <stdio.h>
#include <stdlib.h>
#include <time.h>

extern int omp_fdma(double*a,double*b,double*c,double*d,double*e,double*f,double**y,int dim);
extern int fdma(double *a, double *b, double *c, double *d, double *e, double *f,double *y, int dim);

int reverse(double*,double*,int);
int get_from_file(double **a,double **b,double **c,double **d,double **e,double **f,int* dim);
int write_to_file(double*,int);

int main(int argc,char *argv[])
{
	int dim = 0;
	//int d_size = sizeof(double);
	double *a,*b,*c,*d,*e,*f,*y;
	
	get_from_file(&a,&b,&c,&d,&e,&f,&dim);
	reverse(b,d,dim);
	
	clock_t start,end;
	start = omp_get_wtime();
	omp_fdma(a,b,c,d,e,f,&y,dim);
	/*
	y = (double*)calloc(dim,sizeof(double));
	fdma(a,b,c,d,e,f,y,dim);
	//*/
	end = omp_get_wtime();
	//printf("%f\n",(double)(end-start)/CLOCKS_PER_SEC);
	//fdma(a,b,c,d,e,f,y,dim);
	
	/*
	int i;
	for (i=0;i<dim;i++)
	{
		printf(" y[%d] = %.4f \n",i,y[i]);
	}
	//*/
	
	write_to_file(y,dim);
	FILE* file = fopen("time.log","a");
    fprintf(file,"dim=%d time=%lf\n",dim,((double)(end-start)));
    fclose(file);
    free(y);
    
    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);
	
	return 0;
}


int get_from_file(double **a,double **b,double **c,double **d,double **e,double **f,int* dim)
{
    FILE *fp;
    fp = fopen("input.txt", "r");
    if(fp == NULL)
    {
        fprintf(stderr,"Cannot open file for reading");
        exit(EXIT_FAILURE);
    }
    fscanf(fp,"%d",dim);
    int i;
    (*a) = (double*)malloc((*dim)*sizeof(double));
    for (i=0;i<(*dim);i++)
    {
        fscanf(fp,"%lf",&(*a)[i]);
    }
    (*b) = (double*)malloc((*dim)*sizeof(double));
    for (i=0;i<(*dim);i++)
    {
        fscanf(fp,"%lf",&(*b)[i]);
    }
    (*c) = (double*)malloc((*dim)*sizeof(double));
    for (i=0;i<(*dim);i++)
    {
        fscanf(fp,"%lf",&(*c)[i]);
    }
    (*d) = (double*)malloc((*dim)*sizeof(double));
    for (i=0;i<(*dim);i++)
    {
        fscanf(fp,"%lf",&(*d)[i]);
    }
    (*e) = (double*)malloc((*dim)*sizeof(double));
    for (i=0;i<(*dim);i++)
    {
        fscanf(fp,"%lf",&(*e)[i]);
    }
    (*f) = (double*)malloc((*dim)*sizeof(double));
    for (i=0;i<(*dim);i++)
    {
        fscanf(fp,"%lf",&(*f)[i]);
    }
    if(fclose(fp) != 0)
    {
        fprintf(stderr,"Could not close file properly!\n");
        exit(EXIT_FAILURE);
    }
    return 0;
}

int reverse(double *b,double *d,int dim)
{
    int i;
    for (i=0;i<dim;++i)
    {
        b[i]=-b[i];
        d[i]=-d[i];
    }
    return 0;
}

int write_to_file(double*y,int dim)
{
    FILE *file = fopen("output.txt","w");
    if (file== NULL) {
      printf("Error in opening a file..");
      return 1;
    }
    int i=0;
    fprintf(file,"%d\n",dim);
    for (i=0;i<dim;i++)
    {
        fprintf(file,"%lf\n",y[i]);
    }
    fclose(file);
    return 0;
}
