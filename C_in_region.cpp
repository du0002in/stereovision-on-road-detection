// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
#include <omp.h>
#include <math.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ImgPtr, *outPtr, *xli_lbPtr, *xli_rbPtr, *xri_lbPtr, *xri_rbPtr, *yfirstPtr, *ylastPtr;
    mwSize mrows;
    mwSize ncols;
	mwSize y_length;
    int kk, yfirst, ylast;
	//mxArray *U,*V;
    //double *UPtr, *VPtr;

    ImgPtr=mxGetPr(prhs[0]);
    mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	xli_lbPtr=mxGetPr(prhs[1]);
	xli_rbPtr=mxGetPr(prhs[2]);
	xri_lbPtr=mxGetPr(prhs[3]);
	xri_rbPtr=mxGetPr(prhs[4]);
	yfirstPtr=mxGetPr(prhs[5]);
	yfirst=(int)(*yfirstPtr);
	ylastPtr=mxGetPr(prhs[6]);
	ylast=(int)(*ylastPtr);
    plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	plhs[0]=mxDuplicateArray(prhs[0]);
	outPtr=mxGetPr(plhs[0]);
	
	int li_lb_len=mxGetNumberOfElements(prhs[1]);
	int li_rb_len=mxGetNumberOfElements(prhs[2]);
	int ri_lb_len=mxGetNumberOfElements(prhs[3]);
	int ri_rb_len=mxGetNumberOfElements(prhs[4]);

	int i;
	#pragma omp parallel for num_threads(4)
	for (i=0; i<ncols; i++)
	{
		int li_lb, li_rb, ri_lb, ri_rb;
		for (int kk=yfirst; kk<=ylast; kk++)
		{
			int b=kk-yfirst;
			int a=i*mrows+kk;
			if (*(outPtr+a)!=0)
			{
				li_lb=((b>=li_lb_len)?1:xli_lbPtr[b]);
				li_rb=((b>=li_rb_len)?640:xli_rbPtr[b]);
				ri_lb=((b>=ri_lb_len)?671:xri_lbPtr[b]);
				ri_rb=((b>=ri_rb_len)?1310:xri_rbPtr[b]);
				if (!((i>=li_lb && i<=li_rb) || (i>=ri_lb && i<=ri_rb)))
					*(outPtr+a)=0;
			}
		}
	}
}