/*% uniform convolusion
	ImgPtr is the image to be convoluted
	hPtr and wPtr is the convolusion window size height and width
*/
#include "mex.h"
#include <omp.h>
#include <math.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *rowPtr, *colPtr, *ImgPtr, *tempPtr, *win_sizePtr;
	double *out_row_min,*out_col_min, *out_diff;
	int i;
	mwSize mrows, ncols, win_r, win_c;
	mwSize temp_len, rr_cc_len;

	rowPtr=mxGetPr(prhs[0]);
	colPtr=mxGetPr(prhs[1]);
	ImgPtr=mxGetPr(prhs[2]);
	tempPtr=mxGetPr(prhs[3]);
	win_sizePtr=mxGetPr(prhs[4]);

	rr_cc_len=mxGetNumberOfElements(prhs[0]);
	mrows=mxGetM(prhs[2]);
    ncols=mxGetN(prhs[2]);
	temp_len=mxGetM(prhs[3]);
	win_r=*win_sizePtr;
	win_c=*(win_sizePtr+1);
	int win_half_r=win_r/2;
	int win_half_c=win_c/2;
	
	plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[1]=mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2]=mxCreateDoubleMatrix(1, 1, mxREAL);
	out_row_min=mxGetPr(plhs[0]);
	out_col_min=mxGetPr(plhs[1]);
	out_diff=mxGetPr(plhs[2]);

	double *logdiff=new double[rr_cc_len];
	#pragma omp parallel for num_threads(4)
	for (i=0; i<rr_cc_len; i++)
	{
		int r=(int)(*(rowPtr+i))-1;
		int c=(int)(*(colPtr+i))-1;
		double diff_count=0, pixel, pixel_count=0;
		for (int j=0; j<win_c; j++)
		{
			for (int k=0; k<win_r; k++)
			{
				pixel=*(ImgPtr+(r-win_half_r+k)+(c-win_half_c+j)*mrows);
				if (pixel>=0)
				{
					diff_count=diff_count+((pixel==*(tempPtr+k+j*win_r+i*temp_len))?0:1);
					pixel_count++;
				}
			}
		}
		logdiff[i]=diff_count;
	}

	*out_row_min=*rowPtr;
	*out_col_min=*colPtr;
	*out_diff=logdiff[0];

	for (i=1; i<rr_cc_len; i++)
	{
		if (logdiff[i]<logdiff[i-1])
		{
			*out_row_min=rowPtr[i];
			*out_col_min=colPtr[i];
			*out_diff=logdiff[i];
		}
	}

	delete [] logdiff;
}