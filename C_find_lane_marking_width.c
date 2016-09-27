// mex file to calculate zero mean normalized cross correlation

#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *IntensityPtr, *lbPtr, *rbPtr,*tdx1Ptr, *tdx2Ptr, *tdx3Ptr, *widthPtr, *xPtr;
	mwSize mrows, ncols;
	double f, sum_f, sum_fsqr, ave_f, ave_fsqr, var_f;
	double g1, sum_g1, sum_gsqr1, ave_g1, sum_fg1, var_g1, cov_fg1;
	double g2, sum_g2, sum_gsqr2, ave_g2, sum_fg2, var_g2, cov_fg2;
	double g3, sum_g3, sum_gsqr3, ave_g3, sum_fg3, var_g3, cov_fg3;	
	double c, cpre=0;
	int i,j,k,p,len_r;

	IntensityPtr=mxGetPr(prhs[0]);
	lbPtr=mxGetPr(prhs[1]);
	rbPtr=mxGetPr(prhs[2]);
	tdx1Ptr=mxGetPr(prhs[3]);
	tdx2Ptr=mxGetPr(prhs[4]);
	tdx3Ptr=mxGetPr(prhs[5]);
	xPtr=mxGetPr(prhs[6]);
	mrows=mxGetM(prhs[0]);
	ncols=mxGetN(prhs[0]);

	plhs[0]=mxCreateDoubleMatrix(mrows, 1, mxREAL);
	widthPtr=mxGetPr(plhs[0]);

	for (i=0;i<mrows;i++)
	{
		sum_f=0; sum_fsqr=0; cpre=0;
		for (j=(*(lbPtr+i)); j<=(*(rbPtr+i)); j++) {
			f=*(IntensityPtr+mrows*j+i);
			sum_f=sum_f+f;
			sum_fsqr=sum_fsqr+f*f;
		}
		len_r=(*(rbPtr+i)-*(lbPtr+i)+1);
		var_f=pow(sum_fsqr-sum_f*sum_f/len_r,0.5);

		sum_g1=0; sum_gsqr1=0; sum_fg1=0;
		sum_g2=0; sum_gsqr2=0; sum_fg2=0;
		sum_g3=0; sum_gsqr3=0; sum_fg3=0;
		for (p=0;p<len_r;p++)
		{
			f=*(IntensityPtr+mrows*((int)((*(lbPtr+i))+p))+i);
			if ((*(lbPtr+i)+p)<(*(xPtr+i)-*(tdx1Ptr+i)) || (*(lbPtr+i)+p)>(*(xPtr+i)+*(tdx1Ptr+i)))
				g1=0; // template intensity level, out of the boundary is 0, in the boundary is 255
			else g1=255;
			sum_g1=sum_g1+g1;
			sum_gsqr1=sum_gsqr1+g1*g1;
			sum_fg1=sum_fg1+f*g1;

			if ((*(lbPtr+i)+p)<(*(xPtr+i)-*(tdx2Ptr+i)) || (*(lbPtr+i)+p)>(*(xPtr+i)+*(tdx2Ptr+i)))
				g2=0;
			else g2=255;
			sum_g2=sum_g2+g2;
			sum_gsqr2=sum_gsqr2+g2*g2;
			sum_fg2=sum_fg2+f*g2;

			if ((*(lbPtr+i)+p)<(*(xPtr+i)-*(tdx3Ptr+i)) || (*(lbPtr+i)+p)>(*(xPtr+i)+*(tdx3Ptr+i)))
				g3=0;
			else g3=255;
			sum_g3=sum_g3+g3;
			sum_gsqr3=sum_gsqr3+g3*g3;
			sum_fg3=sum_fg3+f*g3;
		}
		var_g1=pow(sum_gsqr1-sum_g1*sum_g1/len_r,0.5);
		cov_fg1=sum_fg1-sum_f*sum_g1/len_r;
		var_g2=pow(sum_gsqr2-sum_g2*sum_g2/len_r,0.5);
		cov_fg2=sum_fg2-sum_f*sum_g2/len_r;
		var_g3=pow(sum_gsqr3-sum_g3*sum_g3/len_r,0.5);
		cov_fg3=sum_fg3-sum_f*sum_g3/len_r;
		c=cov_fg1/var_f/var_g1;
		*(widthPtr+i)=0.10;
		if (cov_fg2/var_f/var_g2>c)
		{
			c=cov_fg2/var_f/var_g2;
			*(widthPtr+i)=0.15;
		}
		if (cov_fg3/var_f/var_g3>c)
		{
			c=cov_fg3/var_f/var_g3;
			*(widthPtr+i)=0.20;
		}
	}
}