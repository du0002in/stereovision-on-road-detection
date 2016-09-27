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
    double *ImgPtr, *p, *H, *focal, *xc, *yc, *theta, *td;
	double *ipm_imgPtr, *final_imgPtr;
	mwSize mrows;
    mwSize ncols;
	int i;
	int row=450, col=300;

	ImgPtr=mxGetPr(prhs[0]);
    mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	p=mxGetPr(prhs[1]);
	H=mxGetPr(prhs[2]);
	focal=mxGetPr(prhs[3]);
	xc=mxGetPr(prhs[4]);
	yc=mxGetPr(prhs[5]);
	td=mxGetPr(prhs[6]);

	plhs[0]=mxCreateDoubleMatrix(row, col, mxREAL);
	ipm_imgPtr=mxGetPr(plhs[0]);

	double H_sinp=*H*(sin(*p));
	double cosp=cos(*p);
	double f_H_cosp=*focal*(*H)*cosp;
	double f_sinp=*focal*(sin(*p));

	// start of inverse persepective mapping
	#pragma omp parallel for num_threads(4)
	for (i=0; i<mrows*ncols; i++)
	{
		double Zw=(*(yc+i)*H_sinp+f_H_cosp)/(f_sinp-*(yc+i)*cosp);
		if (Zw>1.0 && Zw<10.0) {
			double Xw=*(xc+i)*(H_sinp+Zw*cosp)/(*focal);
			if (Xw<3.0 && Xw>-3.0) {
				//int row_ind=(int)(-(Zw-1)*25+225);
				//int col_ind=(int)(-(Xw+3)*50+300);
				int row_ind=(int)(-Zw*50+500);
				int col_ind=(int)(-Xw*50+150);
				*(ipm_imgPtr+row_ind+col_ind*row)=*(ImgPtr+i);
			}
		}
	}
	
	double *non_zero_row=new double[row];
	#pragma omp parallel for num_threads(4)
	for (i=0; i<row; i++)
	{
		for (int j=10; j<300; j=j+30) {
			if (*(ipm_imgPtr+i+j*row)!=0) {
				non_zero_row[i]=1;
				break;
			}
			else {
				non_zero_row[i]=0;
			}
		}
	}
	
	int j=1;
	if (non_zero_row[0]==0) {
		while (non_zero_row[j]==0 && j<(row-1))
			j++;
		for (int i=0; i<j; i++)
			for (int k=0; k<300; k++)
				*(ipm_imgPtr+i+k*row)=*(ipm_imgPtr+j+k*row);
		if (j>=(row-1)) {
			delete [] non_zero_row;
			return;
		}
	}
	else j=0;

	int start_row=j;
	int end_row=j+1;
	for (int i=j; i<(row-1); i++) {
		if (non_zero_row[i]!=0)
			start_row=i;
		else {
			while (non_zero_row[i]==0 && i<(row-1)) {
				i=i+1;
				end_row=i;
			}
			i=i-1;
			if (i<(row-2))
				for (int k=start_row+1; k<end_row; k++)
					for (int p=0; p<300; p++)
						*(ipm_imgPtr+k+p*row)=0.5*(*(ipm_imgPtr+start_row+p*row)+*(ipm_imgPtr+end_row+p*row));
		}
	}
	
	delete [] non_zero_row;

	int total=row*col;	
	plhs[1]=mxCreateDoubleMatrix(row, col, mxREAL);
	final_imgPtr=mxGetPr(plhs[1]);
	
	#pragma omp parallel for num_threads(4)
	for (int i=0; i<total; i++)
	{
		double a=*(ipm_imgPtr+i);
		if (a>*td)
			*(final_imgPtr+i)=1;
		else if (a>0)
			*(final_imgPtr+i)=0;
		else *(final_imgPtr+i)=-1;
		//*(final_imgPtr+i)=(*(ipm_imgPtr+i)>*td?1:0);
	}
}