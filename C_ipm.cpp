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

void max_value(int *first, int *last, int *largest)
{
	if (first==last) {
		*largest=*first;
		return;
	}
	*largest=*first;
	while (++first!=last)
		if (*largest<*first)
			*largest=*first;
}

void min_value(int *first, int *last, int *smallest)
{
	if (first==last) {
		*smallest=*first;
		return;
	}
	*smallest=*first;
	while (++first!=last)
		if (*smallest>*first)
			*smallest=*first;
}

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
	theta=mxGetPr(prhs[6]);
	td=mxGetPr(prhs[7]);

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

	// start of coordinate rotation
	double R[4]; int total=row*col;
	R[0]=cos(*theta); R[1]=-sin(*theta); R[2]=-R[1]; R[3]=R[0];
	double *R0x=new double[col];
	double *R1y=new double[row];
	double *R2x=new double[col];
	double *R3y=new double[row];
	int *outPtr=new int[total*2];

	R0x[0]=R[0]; R1y[0]=R[1];
	R2x[0]=R[2]; R3y[0]=R[3];

	//#pragma omp parallel for num_threads(4)
	for (int i=1; i<col; i++)
	{
		R0x[i]=R0x[i-1]+R[0];	R1y[i]=R1y[i-1]+R[1];
		R2x[i]=R2x[i-1]+R[2];	R3y[i]=R3y[i-1]+R[3];
	}
	for (int i=col; i<row; i++)
	{
		R1y[i]=R1y[i-1]+R[1];
		R3y[i]=R3y[i-1]+R[3];
	}

	#pragma omp parallel for num_threads(4)
	for (int i=0; i<col; i++)
	{
		int rowind=i*row;
		int rowind2=i*row+total;
		for (int j=0; j<row; j++)
		{
			*(outPtr+j+rowind)=(int)(R0x[i]+R1y[j]);
			*(outPtr+j+rowind2)=(int)(R2x[i]+R3y[j]);
		}
	}

	delete [] R0x;
	delete [] R1y;
	delete [] R2x;
	delete [] R3y;
	

	// start of mapping ipm_img to the final image
	int max_x, min_x, max_y, min_y;
	max_value(outPtr, outPtr+total, &max_x);
	min_value(outPtr, outPtr+total, &min_x);
	max_value(outPtr+total, outPtr+2*total, &max_y);
	min_value(outPtr+total, outPtr+2*total, &min_y);

	int move_x=0; int move_y=0;

	if (min_x<=0) {
		move_x=-min_x+1;
		#pragma omp parallel for num_threads(4)
		for (int i=0; i<total; i++)
		{
			*(outPtr+i)=*(outPtr+i)+move_x;
		}
	}
	if (min_y<=0) {
		move_y=-min_y+1;
		#pragma omp parallel for num_threads(4)
		for (int i=0; i<total; i++)
		{
			*(outPtr+i+total)=*(outPtr+i+total)+move_y;
		}
	}

	max_x=max_x+move_x+1;
	max_y=max_y+move_y+1;

	plhs[1]=mxCreateDoubleMatrix(max_y, max_x, mxREAL);
	final_imgPtr=mxGetPr(plhs[1]);
	plhs[2]=mxCreateDoubleMatrix(1, 2, mxREAL);
	double *move_x_yPtr;
	move_x_yPtr=mxGetPr(plhs[2]);
	*move_x_yPtr=move_x; *(move_x_yPtr+1)=move_y;
	#pragma omp parallel for num_threads(4)
	for (int i=0; i<total; i++)
	{
		*(final_imgPtr+(*(outPtr+i+total))+(*(outPtr+i))*max_y)=
			(*(ipm_imgPtr+i)>*td?1:0);
	}
	delete [] outPtr;
}