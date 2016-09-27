/*% L1 norm
	ImgPtr is the image to be worked on
	indPtr is the index of ImgPtr where L1 norm is performed
	templatePtr is the comparison based to calculate L1 norm
	outPtr is 2x1 array [max_L1 norm; index in Img]
*/
#include "mex.h"
#include <omp.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *ImgPtr, *indPtr, *templatePtr;
	double *outPtr;
	mwSize mrows, temprows, indlen;
    mwSize ncols, tempcols;
	int i;

	ImgPtr=mxGetPr(prhs[0]);
	mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	indPtr=mxGetPr(prhs[1]);
	indlen=mxGetM(prhs[1]);
	templatePtr=mxGetPr(prhs[2]);
	temprows=mxGetM(prhs[2]);
    tempcols=mxGetN(prhs[2]);

	int temprows_h=temprows/2;
	int tempcols_h=tempcols/2;

	plhs[0]=mxCreateDoubleMatrix(indlen,1, mxREAL);
	outPtr=mxGetPr(plhs[0]);

	#pragma omp parallel for num_threads(4)
	for (i=0; i<indlen; i++)
	{
		int ind=(int)(*(indPtr+i))-1;
		int row=ind%mrows;
		int col=ind/mrows;
		double sum_diff=0; double sum_counter=0; int counter=0;
		int start_row=row-temprows_h; int start_col=col-tempcols_h;
		double *ImgPtr_k_mrows;
		for (int k=start_col; k<start_col+tempcols; k++) {
			if (k>=0 && k<ncols) {
				ImgPtr_k_mrows=ImgPtr+k*mrows;
				for (int j=start_row; j<start_row+temprows; j++) {
					if (j>=0 && j<mrows && *(ImgPtr_k_mrows+j)!=-1 && *(templatePtr+counter)!=-1) {
//						if (*(ImgPtr_k_mrows+j)!=-1) {
						sum_diff=sum_diff+(((int)(*(ImgPtr_k_mrows+j))==(int)(*(templatePtr+counter)))?0:1);
						sum_counter++;
//						}
					}
					counter++;
				}
			}
			else counter=counter+temprows;
		}
		*(outPtr+i)=(sum_counter>883)?(sum_diff/sum_counter):1.0;
	}
}