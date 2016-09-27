/*% uniform convolusion
	ImgPtr is the image to be convoluted
	hPtr and wPtr is the convolusion window size height and width
*/
#include "mex.h"
#include <omp.h>
#include <math.h>

void dotproduct(double *temp, double *eigVec, double *coe)
{
	for (int i=0; i<20; i++)
	{
		coe[i]=0;
		for (int j=0; j<5166; j++)
			coe[i] += temp[j]*eigVec[j+i*5166];
	}
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *ImgPtr, *win_sizePtr, *eigPtr, *coePtr, *pca_meanPtr;
	double *outPtr,*outPtr2;
	mwSize mrows, win_r;
    mwSize ncols, win_c;
	int eig_len, pca_len;
	const mxArray *templatePtr;
	int i;


	ImgPtr=mxGetPr(prhs[0]);
	mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	win_sizePtr=mxGetPr(prhs[1]);
	win_r=*win_sizePtr;
	win_c=*(win_sizePtr+1);
	eigPtr=mxGetPr(prhs[2]);
	eig_len=win_r*win_c;
	pca_len=mxGetN(prhs[2]);
	coePtr=mxGetPr(prhs[3]);
	pca_meanPtr=mxGetPr(prhs[4]);

	int win_half_r=win_r/2+1;
	int win_half_c=win_c/2+1;


	plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	outPtr=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	outPtr2=mxGetPr(plhs[1]);

	//#pragma omp parallel for num_threads(4)
	for (i=0; i<mrows*ncols; i++)
	{
		int r=i%mrows;
		int c=i/mrows;
		if (r>=win_half_r && r<=(mrows-win_half_r) && c>=win_half_c && c<=(ncols-win_half_c)
			&& r%2==0 && c%2==0)
		{
			double temp[5166], atmp[20], win_counter=0, cent_r=0, cent_c=0;
			for (int j=0; j<win_c; j++)
			{
				for (int k=0; k<win_r; k++)
				{
					temp[k+j*win_r]=*(ImgPtr+(r-win_half_r+k)+(c-win_half_c+j)*mrows)-*(pca_meanPtr+k+j*win_r);
					win_counter=win_counter+*(ImgPtr+(r-win_half_r+k)+(c-win_half_c+j)*mrows);
					cent_r=cent_r+k*(*(ImgPtr+(r-win_half_r+k)+(c-win_half_c+j)*mrows));
					cent_c=cent_c+j*(*(ImgPtr+(r-win_half_r+k)+(c-win_half_c+j)*mrows));
				}
			}
			cent_r=cent_r/win_counter;
			cent_c=cent_c/win_counter;

			if (win_counter>280) {
				if (abs(34-cent_r)<=10 || abs(41-cent_c)<=10)
				{
					dotproduct(temp, eigPtr, atmp);
					double diff_sq[91];
					double min_diff=5166,min_diff_ind;
					for (int j=0; j<91; j++)
					{
						diff_sq[j]=0;
						for (int k=0; k<20; k++)
						{
							diff_sq[j] = diff_sq[j]+pow((atmp[k]-*(coePtr+k+20*j)),2);
						}
						if (diff_sq[j]<=min_diff) {
							min_diff=diff_sq[j];
							min_diff_ind=j;
						}
					}
					*(outPtr+i)=min_diff;
					*(outPtr2+i)=min_diff_ind;
				}
				else *(outPtr+i)=290.0;
			}
			else *(outPtr+i)=280.0;
		}
		else *(outPtr+i)=206.0;
	}
}