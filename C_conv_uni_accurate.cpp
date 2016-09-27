/*% uniform convolusion
	ImgPtr is the image to be convoluted
	hPtr and wPtr is the convolusion window size height and width
*/
#include "mex.h"
#include <omp.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *ImgPtr, *hPtr, *wPtr;
	double *outPtr, *col_sum, *counter_sum;
	mwSize mrows;
    mwSize ncols;
	int i;

	ImgPtr=mxGetPr(prhs[0]);
	mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	hPtr=mxGetPr(prhs[1]);
	wPtr=mxGetPr(prhs[2]);
//	int h=(int)(*hPtr);
//	int w=(int)(*wPtr);
	int half_h=(int)((*hPtr)/2);
	int half_w=(int)((*wPtr)/2);
	plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	outPtr=mxGetPr(plhs[0]);
	mxArray *col_sumPtr, *counterPtr;
	col_sumPtr=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	col_sum=mxGetPr(col_sumPtr);
	counterPtr=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	counter_sum=mxGetPr(counterPtr);
	//row_sum=mxCreateDoubleMatrix(mrows, ncols, mxREAL);

	#pragma omp parallel for num_threads(4)
	for (i=0; i<ncols; i++)
	{
		int i_mrows=i*mrows;
		double *col_sum_i_mrows=col_sum+i_mrows;
		double *ImgPtr_i_mrows=ImgPtr+i_mrows;
		double *ImgPtr_i_mrows_p_half_h=ImgPtr_i_mrows+half_h;
		double *ImgPtr_i_mrows_m_half_h=ImgPtr_i_mrows-half_h;
		double sum_tmp=0;
		for (int j=0; j<=(half_h); j++)
			sum_tmp=*(ImgPtr_i_mrows+j)+sum_tmp;
		*(col_sum_i_mrows)=sum_tmp;
		for (int j=1; j<=(half_h); j++)
			*(col_sum_i_mrows+j)=*(col_sum_i_mrows+j-1)+*(ImgPtr_i_mrows_p_half_h+j);
		for (int j=(half_h+1); j<(mrows-half_h); j++)
			*(col_sum_i_mrows+j)=*(col_sum_i_mrows+j-1)+*(ImgPtr_i_mrows_p_half_h+j)-*(ImgPtr_i_mrows_m_half_h+j-1);
		for (int j=(mrows-half_h); j<mrows; j++)
			*(col_sum_i_mrows+j)=*(col_sum_i_mrows+j-1)-*(ImgPtr_i_mrows_m_half_h+j-1);
	}

	#pragma omp parallel for num_threads(4)
	for (i=0; i<mrows; i++)
	{
		double *outPtr_i=outPtr+i;
		double *col_sum_i=col_sum+i;
		double sum_tmp=0;
		for (int j=0; j<=(half_w); j++)
			sum_tmp=sum_tmp+*(col_sum_i+j*mrows);
		*(outPtr_i)=sum_tmp;
		for (int j=1; j<=(half_w); j++)
			*(outPtr_i+j*mrows)=*(outPtr_i+(j-1)*mrows)+*(col_sum_i+(j+half_w)*mrows);
		for (int j=(half_w+1); j<(ncols-half_w); j++)
			*(outPtr_i+j*mrows)=*(outPtr_i+(j-1)*mrows)+*(col_sum_i+(j+half_w)*mrows)-*(col_sum_i+(j-1-half_w)*mrows);
		for (int j=(ncols-half_w); j<ncols; j++)
			*(outPtr_i+j*mrows)=*(outPtr_i+(j-1)*mrows)-*(col_sum_i+(j-1-half_w)*mrows);
	}
	mxDestroyArray(col_sumPtr);
	/*
	int total=ncols*mrows;
	#pragma omp parallel for num_threads(4)
	for (i=0; i<total; i++)
	{
		int row=i%mrows;
		int col=i/mrows;
		double counter=0;
		int sta_col, end_col, sta_row, end_row;
		sta_col=((col-half_w)<0?0:(col-half_w));
		end_col=((col+half_w)<ncols?(col+half_w):ncols);
		sta_row=((row-half_h)<0?0:(row-half_h));
		end_row=((row+half_h)<mrows?(row+half_h):mrows);
		for (int p=sta_col; p<end_col; p++){
			for (int q=sta_row; q<end_row; q++) {
				if (*(negImgPtr+i)!=-1)
					counter++;
			}
		}
		*(outPtr+i)=(counter==0?0:(*(outPtr+i)/counter));
	}
	*/
}