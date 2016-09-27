// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
#include <omp.h>
#include <math.h>
//#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *ImgPtr, *tdxPtr, *widthPtr, *validRowPtr;
	double *outPtr;
	mwSize mrows;
	mwSize ncols;
	int i;

	ImgPtr=mxGetPr(prhs[0]);
	mrows=mxGetM(prhs[0]);
	ncols=mxGetN(prhs[0]);
	widthPtr=mxGetPr(prhs[1]);
	int width=*widthPtr;
	int s=(int)(1.5*width); 
	int e=ncols-s-1;
	validRowPtr=mxGetPr(prhs[2]);

	plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	outPtr=mxGetPr(plhs[0]);



	#pragma omp parallel for num_threads(4)
	for (i=0; i<mrows; i++)
	{
		if (validRowPtr[i]!=0)
		{
			for (int j=s; j<e; j++)
			{
				int rg_len=2*width+1;
				double *rg=new double[rg_len];
				double max=*(ImgPtr+i+j*mrows);
				double max_ind=width;

				for (int k=0; k<rg_len; k++) {
					rg[k]=*(ImgPtr+i+(j-width+k)*mrows);
					if (rg[k]>max) {
						max=rg[k];
						max_ind=k;
						break;
					}
				}
				if (max_ind==width) {
					double *ll=new double[width];
					double *rr=new double[width];
					double min_l=rg[width]; double min_r=rg[width];
					int ct_ll=0, ct_rr=0;
					for (int h=0; h<width; h++) {
						ll[h]=rg[h+1]-rg[h];
						rr[h]=rg[width+h]-rg[width+h+1];
						if (rg[h]<min_l) min_l=rg[h];
						if (rg[width+h+1]<min_r) min_r=rg[width+h+1];
						if (ll[h]<=-3) {
							ct_ll++;
							break;
						}
						if (rr[h]<=-3) {
							ct_rr++;
							break;
						}
					}
					if (ct_ll<=0 && ct_rr<=0 && (max-(min_l>min_r?min_l:min_r))>20)
					{
						*(outPtr+i+j*mrows)=1;
						j=j+width;
					}
					else *(outPtr+i+j*mrows)=0;

					delete [] ll;
					delete [] rr;
				}
				else *(outPtr+i+j*mrows)=0;
				delete [] rg;
			}
		}
	}
}