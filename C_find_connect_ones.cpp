// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
//#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *arrayPtr, *conditionPtr;
	double *outPtr;
	mwSize no_ele;
	int i,j,s=0,e=0;
	double count_f=0;


	arrayPtr=mxGetPr(prhs[0]);
	no_ele=mxGetNumberOfElements(prhs[0]);
	conditionPtr=mxGetPr(prhs[1]);
	plhs[0]=mxCreateDoubleMatrix(1, 2, mxREAL);
	outPtr=mxGetPr(plhs[0]);

	for (i=0; i<no_ele; i++)
	{
		if (arrayPtr[i]==1)
		{
			double count_1=1;
			for (j=i+1; j<no_ele; j++)
			{
				if (arrayPtr[j]==0) break;
				else count_1++;
			}
			if (count_1>=*conditionPtr)
			{
				if (count_1>count_f)
				{
					count_f=count_1;
					s=i; e=j-1;
				}
			}
			i=j;
		}
	}
	outPtr[0]=(double)s+1;
	outPtr[1]=(double)e+1;
}