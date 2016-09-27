// mex file to calculate gradient
// input:   mxn gray image
// output:  gradient along x (mxn) and y (mxn)
// http://www.stat.cmu.edu/~ryantibs/median/quickselect.c

#include "mex.h"
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *arr, *kPtr, *out;
	int n;
	unsigned k;
	mxArray *tmp;

	tmp=mxDuplicateArray(prhs[0]);
	arr=mxGetPr(tmp);
	n=mxGetNumberOfElements(prhs[0]);
	kPtr=mxGetPr(prhs[1]);
	k=(unsigned)(*kPtr);
	plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
	out=mxGetPr(plhs[0]);

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
	  return;
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

}