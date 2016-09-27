// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
#include <omp.h>
#include <math.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ImgPtr, *outPtr, *xlPtr, *xrPtr, *yfirstPtr, *ylastPtr;
    mwSize mrows;
    mwSize ncols;
	mwSize y_length;
    int kk, yfirst, ylast;
	//mxArray *U,*V;
    //double *UPtr, *VPtr;

    ImgPtr=mxGetPr(prhs[0]);
    mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	xlPtr=mxGetPr(prhs[1]);
	xrPtr=mxGetPr(prhs[2]);
	yfirstPtr=mxGetPr(prhs[3]);
	yfirst=(int)(*yfirstPtr);
	ylastPtr=mxGetPr(prhs[4]);
	ylast=(int)(*ylastPtr);
    plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	plhs[0]=mxDuplicateArray(prhs[0]);
	outPtr=mxGetPr(plhs[0]);

	#pragma omp parallel for num_threads(4)
	for (kk=yfirst; kk<=ylast; kk++)
	//for (kk=*yPtr;kk<=*(yPtr+y_length);kk++)
	{
		for (int i=0; i<ncols; i++)
		{
			int a=i*mrows+kk;
			int lilb=*(xlPtr+kk-yfirst)-100;
			int lirb=(lilb+200<670)?(lilb+200):670;
			int rirb=*(xrPtr+kk-yfirst)+100;
			int rilb=(rirb-200>640)?(rirb-200):640;
			if (!((i>=lilb && i<=lirb) || (i>=rilb && i<=rirb)))
				*(outPtr+a)=0;
		}
	}
}