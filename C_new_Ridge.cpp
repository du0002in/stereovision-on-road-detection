// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
#include <omp.h>
#include <math.h>
//#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ImgPtr, *tdxPtr, *OffsetPtr;
	double *outPtr;
	mwSize mrows;
    mwSize ncols;
	int i;

	ImgPtr=mxGetPr(prhs[0]);
    mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	tdxPtr=mxGetPr(prhs[1]);
	OffsetPtr=mxGetPr(prhs[2]);

	plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	outPtr=mxGetPr(plhs[0]);

	int Offset=(int)*OffsetPtr;
	int dumy_region_s=ncols/2-15;
	int dumy_region_e=ncols/2+15;

	//#pragma omp parallel for num_threads(4)
	for (i=Offset/2; i<(mrows-5); i++)
	{
		int width=2*(*(tdxPtr+(2*i-Offset)));
		int s=(int)(1.5*width); 
		int e=ncols-s-1;
		int corupt_region_s=dumy_region_s-s;
		int corupt_region_e=dumy_region_e+s;

		#pragma omp parallel for num_threads(4)
		for (int j=s; j<e; j++)
		{
			if (j>=corupt_region_s && j<=corupt_region_e)
				*(outPtr+i+j*mrows)=0;
			else {
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
					//if (abs(rg[0]-rg[2*width])<=25.0) {
						double *ll=new double[width];
						double *rr=new double[width];
						double min_l=rg[width]; double min_r=rg[width];
						int ct_ll=0, ct_rr=0;
						for (int h=0; h<width; h++) {
							ll[h]=rg[h+1]-rg[h];
							rr[h]=rg[width+h]-rg[width+h+1];
							//if (abs(ll[h])<3) ll[h]=0;
							//if (abs(rr[h])<3) rr[h]=0;
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
						if (ct_ll<=0 && ct_rr<=0 && (max-(min_l>min_r?min_l:min_r))>30)
							*(outPtr+i+j*mrows)=1;
						else *(outPtr+i+j*mrows)=0;

						delete [] ll;
						delete [] rr;
					//}
					//else *(outPtr+i+j*mrows)=0;
				}
				else *(outPtr+i+j*mrows)=0;

				delete [] rg;
			}
		}
	}

}