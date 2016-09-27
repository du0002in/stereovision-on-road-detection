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
	double *outPtr, *test_xPtr, *test_yPtr, *vert_xPtr, *vert_yPtr;
	mwSize testlen, mrows;
	int i, counter=0;

	test_xPtr=mxGetPr(prhs[0]);
	test_yPtr=mxGetPr(prhs[1]);
	vert_xPtr=mxGetPr(prhs[2]);
	vert_yPtr=mxGetPr(prhs[3]);

	testlen=mxGetNumberOfElements(prhs[0]);

	plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
	outPtr=mxGetPr(plhs[0]);


	for (i=0; i<testlen; i++)
	{
		int c;
		c=pnpoly(4, vert_xPtr, vert_yPtr, test_xPtr[i], test_yPtr[i]);
		if (c==1)	counter++;
		if (counter>=10) {
			*outPtr=1;
			return;
		}
	}
	*outPtr=0;

}