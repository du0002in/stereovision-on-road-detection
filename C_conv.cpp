/*% uniform convolusion
	ImgPtr is the image to be convoluted
	hPtr and wPtr is the convolusion window size height and width
*/
#include "mex.h"
#include <omp.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *ImgPtr, *indPtr, *tempPtr;
	double *outPtr;
	mwSize mrows, indlen, trows_h, trows;
    mwSize ncols, tcols_h, tcols, temp_num;
	const mxArray *templatePtr;
	int i;

	/*
	ImgPtr=mxGetPr(prhs[0]);
	mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	indPtr=mxGetPr(prhs[1]);
	indlen=mxGetM(prhs[1]);
	tempPtr=mxGetPr(prhs[2]);
	trows=mxGetM(prhs[2]);
	tcols=mxGetN(prhs[2]);
	trows_h=trows/2;
	tcols_h=tcols/2;

	plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	outPtr=mxGetPr(plhs[0]);
	
	#pragma omp parallel for num_threads(4)
	for (i=0; i<indlen; i++)
	{
		int ind=(int)(*(indPtr+i));
		int row=ind%mrows;
		int col=ind/mrows;
		int sta_col, end_col, sta_row, end_row;
		double counter=0;
		sta_col=((col-tcols_h)<0?0:(col-tcols_h));
		end_col=((col+tcols_h)<ncols?(col+tcols_h):ncols);
		sta_row=((row-trows_h)<0?0:(row-trows_h));
		end_row=((row+trows_h)<mrows?(row+trows_h):mrows);
		for (int p=sta_col; p<end_col; p++)
			for (int q=sta_row; q<end_row; q++)
				counter=counter+((*(ImgPtr+q+p*mrows)==1 && *(tempPtr+(p-col+tcols_h)*trows_h+(q-row+trows_h))==1)?1:0);
		*(outPtr+ind)=counter;
	}
	*/
	ImgPtr=mxGetPr(prhs[0]);
	mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	indPtr=mxGetPr(prhs[1]);
	indlen=mxGetM(prhs[1]);
	templatePtr=prhs[2];
	temp_num=mxGetNumberOfElements(prhs[2]);
	mwSize dims[3];
	dims[0]=mrows; dims[1]=ncols; dims[2]=temp_num;

	plhs[0]=mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	outPtr=mxGetPr(plhs[0]);

	mxArray *lhs[13], *rhs[3];
	
	rhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	double *rhs_p0;
	rhs_p0=mxGetPr(rhs[0]);
	#pragma omp parallel for num_threads(4)
	for (i=0; i<mrows*ncols; i++)
	{
		rhs_p0[i]=ImgPtr[i];
	}
	
	//#pragma omp parallel for num_threads(4)
	for (int i=0; i<temp_num; i++) {
		const mxArray *cellArray;
		double *p;
		cellArray=mxGetCell(templatePtr,i);
		p=mxGetPr(cellArray);
		int temp_rows=mxGetM(cellArray);
		int temp_cols=mxGetN(cellArray);
		rhs[1]=mxCreateDoubleMatrix(temp_rows, temp_cols, mxREAL);
		double *rhs_p1;	
		rhs_p1=mxGetPr(rhs[1]);
		for (int j=0; j<temp_cols*temp_rows; j++)
			rhs_p1[j]=p[j];
		rhs[2]=mxCreateString("same");
		lhs[i]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
		mexCallMATLAB(1, &(lhs[i]), 3, rhs, "conv2");
		double *lhs_p;
		lhs_p=mxGetPr(lhs[i]);
		#pragma omp parallel for num_threads(4)
		for (int j=0; j<mrows*ncols; j++)
			outPtr[j+i*mrows*ncols]=lhs_p[j];
	}

}