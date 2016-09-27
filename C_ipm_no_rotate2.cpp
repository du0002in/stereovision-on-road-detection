/*% inverse perspective mapping
% img is the retified and half resoluted image, 
% p is view angle Phi, 
% H is cam height, 
% focal is cam focal length, 
% xxc and yyc are cam coordinate grid defined as following:
%	[ximg,yimg]=meshgrid(1:640,1:2:480);
%	yyc=-yimg+cc_y;
%	xxc=-ximg+cc_x_left;
% theta is vehicle moving direction w.r.t road tangent
% output 1: is the iverse perspective image
% output 2: is the rotated inverse perspective image
*/
#include "mex.h"
#include <omp.h>
#include <math.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ImgPtr, *p, *H, *focal, *theta, *td, *cc_x, *cc_y;
	double *ipm_imgPtr, *final_imgPtr;
	mwSize mrows;
    mwSize ncols;
	int i;
	int row=450, col=300;
	int total=row*col;

	ImgPtr=mxGetPr(prhs[0]);
    mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	p=mxGetPr(prhs[1]);
	H=mxGetPr(prhs[2]);
	focal=mxGetPr(prhs[3]);
	td=mxGetPr(prhs[4]);
	cc_x=mxGetPr(prhs[5]);
	cc_y=mxGetPr(prhs[6]);

	plhs[0]=mxCreateDoubleMatrix(row, col, mxREAL);
	ipm_imgPtr=mxGetPr(plhs[0]);

	plhs[1]=mxCreateDoubleMatrix(row, col, mxREAL);
	final_imgPtr=mxGetPr(plhs[1]);

	double H_sinp=*H*(sin(*p));
	double cosp=cos(*p);
	double f_H_cosp=*focal*(*H)*cosp;
	double f_sinp=*focal*(sin(*p));

	// start of inverse persepective mapping
	#pragma omp parallel for num_threads(4)
	for (i=0; i<total; i++)
	{
		int r=i%row;
		int c=i/row;
		double yc=(50*f_H_cosp+(r-500)*f_sinp)/((r-500)*cosp-50*H_sinp);
		int yimg=(int)((-yc+*cc_y)/2.0);
		if (yimg>=0 && yimg<mrows)
		{
			double Zw=(yc*H_sinp+f_H_cosp)/(f_sinp-yc*cosp);
			double xc=-*focal*(c-150.0)/50.0/(H_sinp+Zw*cosp);
			int ximg=-(int)xc+*cc_x;
			if (ximg>=0 && ximg<ncols)
			{
				*(ipm_imgPtr+i)=*(ImgPtr+yimg+ximg*mrows);
				if (*(ipm_imgPtr+i)>*td)
					*(final_imgPtr+i)=1;
				else *(final_imgPtr+i)=0;
			}
			else
			{	
				*(ipm_imgPtr+i)=0;
				*(final_imgPtr+i)=-1;
			}
		}
		else
		{
			*(ipm_imgPtr+i)=0;
			*(final_imgPtr+i)=-1;
		}
	}
}