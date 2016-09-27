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
/*
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
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ImgPtr, *p, *H, *focal, *xc, *yc, *theta, *td;
	double *final_imgPtr;
	mwSize row;
    mwSize col;
	int i;

	ImgPtr=mxGetPr(prhs[0]);
    row=mxGetM(prhs[0]);
    col=mxGetN(prhs[0]);
	theta=mxGetPr(prhs[1]);

	// start of coordinate rotation
	double R[4]; 
	int total=row*col;
	R[0]=cos(*theta); R[1]=-sin(*theta); R[2]=-R[1]; R[3]=R[0];
	//double *R0x=new double[col];
	//double *R1y=new double[row];
	//double *R2x=new double[col];
	//double *R3y=new double[row];
	/*
	R0x[0]=R[0]; R1y[0]=R[1];
	R2x[0]=R[2]; R3y[0]=R[3];
	*/
	int max_x, min_x, max_y, min_y;
	if (*theta>0) {
		min_x=(int)((row-1)*R[1]);	max_x=(int)((col-1)*R[0]);
		min_y=0;					max_y=(int)((col-1)*R[2]+(row-1)*R[3]);
	}
	else {
		min_x=0;					max_x=(int)((col-1)*R[0]+(row-1)*R[1]);
		min_y=(int)((col-1)*R[2]);	max_y=(int)((row-1)*R[3]);
	}
	int move_x=((min_x<0)?(-min_x):0);
	int move_y=((min_y<0)?(-min_y):0);
	max_x=max_x+move_x+1;
	max_y=max_y+move_y+1;

	plhs[0]=mxCreateDoubleMatrix(max_y, max_x, mxREAL);
	final_imgPtr=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(1, 2, mxREAL);
	double *move_x_yPtr;
	move_x_yPtr=mxGetPr(plhs[1]);
	*move_x_yPtr=move_x; *(move_x_yPtr+1)=move_y;
	/*
	R0x[0]=0; R1y[0]=0;
	R2x[0]=0; R3y[0]=0;

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
		for (int j=0; j<row; j++)
		{
			//*(outPtr+j+rowind)=(int)(R0x[i]+R1y[j])+move_x;
			//*(outPtr+j+rowind2)=(int)(R2x[i]+R3y[j])+move_y;
			int x=(int)(R0x[i]+R1y[j])+move_x;
			int y=(int)(R2x[i]+R3y[j])+move_y;
			*(final_imgPtr+y+x*max_y)=*(ImgPtr+rowind+j)+2;
		}
	}	*/
	/*
	#pragma omp parallel for num_threads(4)
	for (int i=0; i<col; i++)
	{
		double R0i=R[0]*i, R2i=R[2]*i;
		int rowind=i*row;
		for (int j=0; j<row; j++)
		{
			int x=(int)(R0i+R[1]*j)+move_x;
			int y=(int)(R2i+R[3]*j)+move_y;
			*(final_imgPtr+y+x*max_y)=*(ImgPtr+rowind+j)+2;
		}
	}*/

	#pragma omp parallel for num_threads(4)
	for (i=0; i<total; i++)
	{
		int r=i%row;
		int c=i/row;
		int x=(int)(R[0]*c+R[1]*r)+move_x;
		int y=(int)(R[2]*c+R[3]*r)+move_y;
		*(final_imgPtr+y+x*max_y)=*(ImgPtr+i)+2;
	}

	total=max_y*max_x-1;
	#pragma omp parallel for num_threads(4)
	for (int i=0; i<total; i++)
	{
		if (*(final_imgPtr+i)==0)
			*(final_imgPtr+i)=*(final_imgPtr+i+1);
	}
	/*
	delete [] R0x;
	delete [] R1y;
	delete [] R2x;
	delete [] R3y;*/
}