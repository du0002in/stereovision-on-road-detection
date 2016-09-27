#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <stdlib.h>
#include <omp.h>

static int colnum = 0;

double** matrix(int rows,int cols){
int k;
double **m;
m = (double **)mxMalloc(rows * sizeof(double *));

for (k=0; k<rows; k++){
    m[k] = (double *)mxMalloc(cols * sizeof(double));
}
return m;
}

int cmp (const void *pa, const void *pb) {
	
	double** x=(double**) pa;
	double** y=(double**) pb;

	double xval, yval;
	xval = *(*(x)+colnum);
	yval = *(*(y)+colnum);

	if (xval<yval) return -1;
	if (xval>yval) return +1;
	return 0;
	/*if (*(double*)pa<*(double*)pb) return -1;
	if (*(double*)pa>*(double*)pb) return +1;
	return 0;*/
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *valuesPtr;
	double *sortedPtr, *indexPtr;
	int nc,n,k;

	valuesPtr=mxGetPr(prhs[0]);
	n=mxGetNumberOfElements(prhs[0]);
	double **z;
	z=matrix(n,2);

	#pragma omp parallel for num_threads(4)
	for (k=0;k<n;k++) 
	{
		z[k][0]=*(valuesPtr+k);
		z[k][1]=k+1;
	}

	plhs[0]=mxCreateDoubleMatrix(n, 1, mxREAL);
	sortedPtr=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(n, 1, mxREAL);
	indexPtr=mxGetPr(plhs[1]);


	qsort(z,n,sizeof(*valuesPtr),cmp);

	#pragma omp parallel for num_threads(4)
	for (k=0;k<n;k++) 
	{
		*(sortedPtr+k)=z[k][0];
		*(indexPtr+k)=z[k][1];
	}
	mxFree(z);
}