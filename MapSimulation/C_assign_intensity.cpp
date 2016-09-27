// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
#include <omp.h>
#include <math.h>
#include <time.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *cumsumPtr, *indPtr, *p_numPtr;
	mwSize ncols;
	int kk, p_num;
	cumsumPtr=mxGetPr(prhs[0]);
	ncols=mxGetN(prhs[0]);
	p_numPtr=mxGetPr(prhs[1]);
	p_num=(int)*(p_numPtr);

	plhs[0]=mxCreateDoubleMatrix(1, p_num, mxREAL);
	indPtr=mxGetPr(plhs[0]);

	srand(time(NULL));
	#pragma omp parallel for num_threads(4)
	for (kk=0; kk<p_num; kk++)
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