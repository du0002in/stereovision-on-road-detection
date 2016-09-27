/*% uniform convolusion
	ImgPtr is the image to be convoluted
	hPtr and wPtr is the convolusion window size height and width
*/
#include "mex.h"
#include <omp.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *ImgPtr, *hPtr, *wPtr, *indPtr;
	double *outPtr, *col_sum, *counter_sum;
	mwSize mrows, indlen;
    mwSize ncols;
	int i;

	ImgPtr=mxGetPr(prhs[0]);
	mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	indPtr=mxGetPr(prhs[1]);
	indlen=mxGetM(prhs[1]);
	hPtr=mxGetPr(prhs[2]);
	wPtr=mxGetPr(prhs[3]);
//	int h=(int)(*hPtr);
//	int w=(int)(*wPtr);
	int half_h=(int)((*hPtr)/2);
	int half_w=(int)((*wPtr)/2);
	plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	outPtr=mxGetPr(plhs[0]);

	#pragma omp parallel for num_threads(4)
	for (i=0; i<indlen; i++)
	{
		int ind=(int)(*(indPtr+i));
		int row=ind%mrows;
		int col=ind/mrows;
		int sta_col, end_col, sta_row, end_row;
		sta_col=((col-half_w)<0?0:(col-half_w));
		end_col=((col+half_w)<ncols?(col+half_w):ncols);
		sta_row=((row-half_h)<0?0:(row-half_h));
		end_row=((row+half_h)<mrows?(row+half_h):mrows);
		for (int p=sta_col; p<end_col; p++) {
			int rowind=p*mrows;
			for (int q=sta_row; q<end_row; q++)
				*(outPtr+q+rowind)=1;
		}
	}
}