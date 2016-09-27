// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#include "mex.h"
#include <omp.h>
#include <math.h>
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ImgPtr, *divPtr;
    mwSize mrows;
    mwSize ncols;
    int kk;
	mxArray *U,*V;
    double *UPtr, *VPtr;

    ImgPtr=mxGetPr(prhs[0]);
    mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
    plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	divPtr=mxGetPr(plhs[0]);

	U=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
    UPtr=mxGetPr(U);
	V=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
    VPtr=mxGetPr(V);

	#pragma omp parallel for num_threads(4)
    for (kk=mrows+1;kk<(ncols*mrows-mrows-1) ;kk++)
    {
            double Grad_x=(*(ImgPtr+kk+mrows)-*(ImgPtr+kk-mrows));
			double Grad_y=(*(ImgPtr+kk+1)-*(ImgPtr+kk-1));
			double ssd11=Grad_x*Grad_x;
			double ssd12=Grad_x*Grad_y;
			double ssd22=Grad_y*Grad_y;
			double Ssdi12=ssd12*.6306+(ssd11+ssd22)*0.0838;
			double Ssdi11_22=(ssd11-ssd22)*.608;
			double EigenValue_22D12=0.5*(Ssdi11_22+sqrt(Ssdi11_22*Ssdi11_22+4*Ssdi12*Ssdi12))/Ssdi12;
			if (Ssdi12==0)
				EigenValue_22D12=0;
			double wpsdiV=sqrt(1/(1+EigenValue_22D12*EigenValue_22D12));
			double wpsdiU=wpsdiV*EigenValue_22D12;
			double psdi=wpsdiU*Grad_x+wpsdiV*Grad_y;
			if (psdi>0) {
				*(UPtr+kk)=wpsdiU;
				*(VPtr+kk)=wpsdiV;
			}
			else if (psdi==0) {
				*(UPtr+kk)=0;
				*(VPtr+kk)=0;
			}
			else {
				*(UPtr+kk)=-wpsdiU;
				*(VPtr+kk)=-wpsdiV;
			}
    }

	#pragma omp parallel for num_threads(4)
    for (kk=mrows+1;kk<(ncols*mrows-mrows-1) ;kk++)
	{
		*(divPtr+kk)=-0.5*(*(UPtr+kk+mrows)-*(UPtr+kk-mrows)+*(VPtr+kk+1)-*(VPtr+kk-1));
	}
	mxDestroyArray(U);
	mxDestroyArray(V);



	mxArray *tmp=mxDuplicateArray(plhs[0]);
	double *arr=mxGetPr(tmp);
	int n=ncols*mrows;
	double *kPtr=mxGetPr(prhs[1]);
	unsigned k=(unsigned)(*kPtr);
	plhs[1]=mxCreateDoubleMatrix(1, 1, mxREAL);
	double *out=mxGetPr(plhs[1]);

	unsigned long i,ir,j,l,mid;
	double a,temp;

  l=0;
  ir=n-1;
  for(;;) {
    if (ir <= l+1) { 
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir]);
      }
      *out=*(arr+k);
	  break;
    }
    else {
      mid=(l+ir) >> 1; 
      SWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1; 
      j=ir;
      a=arr[l+1]; 
      for (;;) { 
	do i++; while (arr[i] < a); 
	do j--; while (arr[j] > a); 
	if (j < i) break; 
	SWAP(arr[i],arr[j]);
      } 
      arr[l+1]=arr[j]; 
      arr[j]=a;
      if (j >= k) ir=j-1; 
      if (j <= k) l=i;
    }
  }
  mxDestroyArray(tmp);
  
  #pragma omp parallel for num_threads(4)
    for (kk=mrows+1;kk<(ncols*mrows-mrows-1) ;kk++)
	{
		//double ridge=*(divPtr+kk);
		//*(divPtr+kk)=(ridge>=*out)?1:0;
		if (*(divPtr+kk)<*out)
			*(divPtr+kk)=0;
	}
}