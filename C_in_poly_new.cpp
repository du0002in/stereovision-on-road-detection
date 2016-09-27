/*% uniform convolusion
	ImgPtr is the image to be convoluted
	hPtr and wPtr is the convolusion window size height and width
*/
#include "mex.h"
#include <omp.h>
#include <math.h>

int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *ImgPtr, *test_xPtr, *test_yPtr, *vert_xPtr, *vert_yPtr;
	mwSize testlen, mrows;
	int i;

	ImgPtr=mxGetPr(prhs[0]);
	test_xPtr=mxGetPr(prhs[1]);
	test_yPtr=mxGetPr(prhs[2]);
	vert_xPtr=mxGetPr(prhs[3]);
	vert_yPtr=mxGetPr(prhs[4]);

	mrows=mxGetM(prhs[0]);
	testlen=mxGetM(prhs[1]);

	#pragma omp parallel for num_threads(4)
	for (i=0; i<testlen; i++)
	{
		int c;
		c=pnpoly(4, vert_xPtr, vert_yPtr, test_xPtr[i], test_yPtr[i]);
		if (c==1)
			*(ImgPtr+((int)(test_xPtr[i]))*mrows+(int)(test_yPtr[i]))=0;
	}

}