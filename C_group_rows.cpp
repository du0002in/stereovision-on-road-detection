/*% uniform convolusion
	ImgPtr is the image to be convoluted
	hPtr and wPtr is the convolusion window size height and width
*/
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *rowArrayPtr, *distPtr;
	double *outPtr;
	mwSize no_ele;
	double dist;

	rowArrayPtr=mxGetPr(prhs[0]);
	no_ele=mxGetNumberOfElements(prhs[0]);
	distPtr=mxGetPr(prhs[1]);
	dist=(*distPtr);

	int total_row;
	total_row=(int)(rowArrayPtr[no_ele-1]-rowArrayPtr[0]+1);
	plhs[0]=mxCreateDoubleMatrix(1, total_row, mxREAL);
	outPtr=mxGetPr(plhs[0]);

	int row_counter=0;
	int start_row=(int)rowArrayPtr[0];
	for (int i=1; i<no_ele; i++)
	{
		double row_diff=rowArrayPtr[i]-rowArrayPtr[i-1];
		if (row_diff<=dist) {
			for (int j=(int)rowArrayPtr[i-1]; j<=(int)rowArrayPtr[i]; j++)
				outPtr[j-start_row]=j;
		}
	}
	double s, e;
	for (int i=0; i<total_row; i++)
	{
		if (outPtr[i]!=0) {
			s=i; e=i;
			int j;
			for (j=i+1; j<total_row; j++) {
				if (outPtr[j]!=0) e=j;
				else break;
			}
			i=j;
			if (e-s+1<9) {
				for (int k=s; k<=e; k++)
					outPtr[k]=0;
			}
		}
	}


}