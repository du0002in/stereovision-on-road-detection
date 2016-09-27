// intersectAB(*A, *B)=> return an binary array corresponding to B element.
// e.g. A=[1 2 3 4 5], B=[3 7 1] => output=[1 0 1]
#include "mex.h"
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *APtr, *BPtr;
	double *indPtr;
	int i;

	APtr=mxGetPr(prhs[0]);
	BPtr=mxGetPr(prhs[1]);
	mwSize Anum=mxGetNumberOfElements(prhs[0]);
	mwSize Bnum=mxGetNumberOfElements(prhs[1]);

	plhs[0]=mxCreateDoubleMatrix(1, Bnum, mxREAL);
	indPtr=mxGetPr(plhs[0]);

	#pragma omp parallel for num_threads(4)
    for(i=0;i<Bnum;i++) {
		for (int j=0; j<Anum; j++) {
			if (*(BPtr+i)==*(APtr+j)) {
				*(indPtr+i)=1;
				break;
			}
		}
	}
}