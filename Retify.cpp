// mex file to calculate gradient
// input:   2 mxn structure tensor field matrix U and V
// output:  mxn divergence matrix
#include "mex.h"
#include <omp.h>

/*
double *spawn_threads(double x1[], double x2[], int cols)
{
	double *sum;
	int i;
	#pragma omp parallel for 
	//shared(x1, x2, cols) private(i)
	for(i=0; i<cols; i++)
		sum[i]=x1[i]+x2[i];
	return sum;
}*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *imgPtr, *a1Ptr, *a2Ptr, *a3Ptr, *a4Ptr, *ind1Ptr, *ind2Ptr, *ind3Ptr, *ind4Ptr, *indnewPtr;
    double *rectimgPtr;
    mwSize mrows;
    mwSize ncols, len;
    int i;
    
	imgPtr=mxGetPr(prhs[0]);
	a1Ptr=mxGetPr(prhs[1]);
	a2Ptr=mxGetPr(prhs[2]);
	a3Ptr=mxGetPr(prhs[3]);
	a4Ptr=mxGetPr(prhs[4]);
	ind1Ptr=mxGetPr(prhs[5]);
	ind2Ptr=mxGetPr(prhs[6]);
	ind3Ptr=mxGetPr(prhs[7]);
	ind4Ptr=mxGetPr(prhs[8]);
	indnewPtr=mxGetPr(prhs[9]);
	len=mxGetN(prhs[1]);
	mrows=mxGetM(prhs[0]);
	ncols=mxGetN(prhs[0]);
	plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	//plhs[0]=mxDuplicateArray(prhs[0]);
	rectimgPtr=mxGetPr(plhs[0]);

	#pragma omp parallel for num_threads(4)
	for (i=0; i<len; i++) {
		*(rectimgPtr+(int)(*(indnewPtr+i)))=*(a1Ptr+i)*(*(imgPtr+(int)(*(ind1Ptr+i))))+*(a2Ptr+i)*(*(imgPtr+(int)(*(ind2Ptr+i))))+*(a3Ptr+i)*(*(imgPtr+(int)(*(ind3Ptr+i))))+*(a4Ptr+i)*(*(imgPtr+(int)(*(ind4Ptr+i))));
	}
}