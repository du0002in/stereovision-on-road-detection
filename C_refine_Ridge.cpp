// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
#include <omp.h>
#include <math.h>
//#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ImgPtr, *OffsetPtr;
	//double *outPtr;
	mwSize mrows;
    mwSize ncols;

	ImgPtr=mxGetPr(prhs[0]);
    mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	OffsetPtr=mxGetPr(prhs[1]);

	//plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	//outPtr=mxGetPr(plhs[0]);

	int Offset=(int)*OffsetPtr;

	#pragma omp parallel for num_threads(4)
	for (int i=Offset/2; i<(mrows-5); i++)
	{
		for (int j=0; j<ncols; j++) {
			int x=0; int counter=0;
			if (*(ImgPtr+i+j*mrows)==1) {
				*(ImgPtr+i+j*mrows)=0;
				x=x+j; counter++; j++;
				while (*(ImgPtr+i+j*mrows)==1 || *(ImgPtr+i+(j+1)*mrows)==1) {
					if (*(ImgPtr+i+j*mrows)==1) {
						*(ImgPtr+i+j*mrows)=0;
						x=x+j; counter++; j++;
					}
					else {
						*(ImgPtr+i+(j+1)*mrows)=0;
						x=x+j+1; counter++; j=j+2;
					}
				}
				x=x/counter;
				*(ImgPtr+i+x*mrows)=1;
			}
		}
	}
}