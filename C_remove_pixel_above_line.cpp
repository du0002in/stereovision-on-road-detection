/*% uniform convolusion
	ImgPtr is the image to be convoluted
	hPtr and wPtr is the convolusion window size height and width
*/
#include "mex.h"
#include <omp.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *linePtr, *ImgPtr;
	mwSize ncols, mrows;
	int i, counter=0;

	ImgPtr=mxGetPr(prhs[0]);
	linePtr=mxGetPr(prhs[1]);
	
	mrows=mxGetM(prhs[0]);
	ncols=mxGetN(prhs[0]);

	#pragma omp parallel for num_threads(4)
	for (i=0; i<ncols; i++)
	{
		for (int j=0; j<mrows; j++)
		{
			if (j<(*(linePtr+i)))
				*(ImgPtr+i*mrows+j)=0;
			else break;
		}
	}
}