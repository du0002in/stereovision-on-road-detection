// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
#include <omp.h>
#include <math.h>
#include <time.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *cumsumPtr, *indPtr;
	mwSize ncols;
	int kk;
	cumsumPtr=mxGetPr(prhs[0]);
	ncols=mxGetN(prhs[0]);
	plhs[0]=mxCreateDoubleMatrix(1, ncols, mxREAL);
	indPtr=mxGetPr(plhs[0]);

	srand(time(NULL));
	#pragma omp parallel for num_threads(4)
	for (kk=0; kk<ncols; kk++)
	{
		double r = ((double) rand() / (RAND_MAX));
		int i;
		for (i=0; i<ncols; i++)
		{
			if (r<=*(cumsumPtr+i))
				break;
		}
		*(indPtr+kk)=(i+1);
	}
}