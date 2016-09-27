// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
#include <omp.h>
#include <math.h>
#include <time.h>

void GetLine(double *x, double *y, double *abc)
{
	if (y[0]!=y[1]) {
		abc[0]=1.0;
		abc[1]=(x[1]-x[0])/(y[0]-y[1]);
		abc[2]=-x[0]-y[0]*abc[1];
	}
	else { 
		abc[0]=0.0; abc[1]=1.0; abc[2]=-y[0];
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *rowPtr, *colPtr, *slope_radPtr;
	double *outPtr, *lineParaPtr;
	double min_slope, max_slope;
	mwSize no_ele;
	int i, no_ite=200, tcon_pre=0;

	rowPtr=mxGetPr(prhs[0]);
	colPtr=mxGetPr(prhs[1]);
	no_ele=mxGetNumberOfElements(prhs[0]);
	slope_radPtr=mxGetPr(prhs[2]);
	min_slope=slope_radPtr[0];
	max_slope=slope_radPtr[1];

	plhs[0]=mxCreateDoubleMatrix(no_ele, 1, mxREAL);
	outPtr=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(3, 1, mxREAL);
	lineParaPtr=mxGetPr(plhs[1]);


	srand(time(NULL));
	#pragma omp parallel for num_threads(4)
	for (i=0; i<no_ite; i++)
	{
		int tcon=0;
		double u[2], v[2], abc[3];
		int rand1=rand() % no_ele;
		int rand2=rand() % no_ele;
		while (rand2==rand1)	rand2=rand() % no_ele;
		v[0]=rowPtr[rand1];
		u[0]=colPtr[rand1];
		v[1]=rowPtr[rand2];
		u[1]=colPtr[rand2];
		GetLine(u,v,abc); // au+bv+c=0 or ax+by+c=0
		double a=abc[0], b=abc[1], c=abc[2];
		double *line_ptr=new double[no_ele];

		double slope=((b==0)?-3.414/2:atan(-a/b));
		if (slope>=min_slope && slope<=max_slope)
		{
			double deno=pow(a*a+b*b, 0.5);
			double dist, x, y;			
			for (int j=0; j<no_ele; j++) {
				if (tcon+no_ele-j>tcon_pre) {
					x=colPtr[j]; y=rowPtr[j];
					dist=abs(a*x+b*y+c)/deno;
					if (dist<1.5) {
						tcon++;
						line_ptr[j]=1;
					}
					else line_ptr[j]=0;
				}
				else break;
			}
		}
		#pragma omp critical
		{
			if (tcon>tcon_pre) {
				tcon_pre=tcon;
				lineParaPtr[0]=a; lineParaPtr[1]=b; lineParaPtr[2]=c; 
				for (int j=0; j<no_ele; j++) 
					outPtr[j]=line_ptr[j];
			}
		}
		delete [] line_ptr;
	}
}