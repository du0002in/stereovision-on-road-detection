// mex file to calculate zero mean normalized cross correlation

#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //input pointer
	double *IntensityRPtr, *CorespLPtr;
	double *li_lbPtr, *li_rbPtr, *ri_lbPtr, *ri_rbPtr;
    //output pointer
	double *disparityPtr;
    mwSize mrows, curcols, precols, nextcols;
    mwSize ncols;
    int i, j, k, p, len_l, len_r;
	double f, sum_f, sum_fsqr, ave_f, ave_fsqr, var_f;
	double g, sum_g, sum_gsqr, ave_g, sum_fg, var_g, cov_fg;
	double c, cpre=0;
    //read in input
	ri_lbPtr=mxGetPr(prhs[0]);
	ri_rbPtr=mxGetPr(prhs[1]);
	li_lbPtr=mxGetPr(prhs[2]);
	li_rbPtr=mxGetPr(prhs[3]);
	IntensityRPtr=mxGetPr(prhs[4]);
	CorespLPtr=mxGetPr(prhs[5]);
	mrows=mxGetM(prhs[5]);
	ncols=mxGetN(prhs[5]);
	//create output matrix
	plhs[0]=mxCreateDoubleMatrix(mrows, 1, mxREAL);
	disparityPtr=mxGetPr(plhs[0]);

	for (i=0; i<mrows; i++) {
		sum_f=0; sum_fsqr=0; cpre=0;
		for (j=(*(ri_lbPtr+i)); j<=(*(ri_rbPtr+i)); j++) {
			f=*(IntensityRPtr+mrows*j+i);
			sum_f=sum_f+f;
			sum_fsqr=sum_fsqr+f*f;
		}
		len_r=(*(ri_rbPtr+i)-*(ri_lbPtr+i)+1);
		var_f=pow(sum_fsqr-sum_f*sum_f/len_r,0.5);

		len_l=(*(li_rbPtr+i)-*(li_lbPtr+i)+1);
		for (k=0; k<len_l-len_r; k++) {
			sum_g=0; sum_gsqr=0; sum_fg=0;
			for (p=0; p<len_r; p++) {
				g=*(CorespLPtr+mrows*(((int)(*(li_lbPtr+i)))+k+p)+i);
				f=*(IntensityRPtr+mrows*(((int)(*(ri_lbPtr+i))+p))+i);
				sum_g=sum_g+g;
				sum_gsqr=sum_gsqr+g*g;
				sum_fg=sum_fg+f*g;
			}
			var_g=pow(sum_gsqr-sum_g*sum_g/len_r,0.5);
			cov_fg=sum_fg-sum_f*sum_g/len_r;
			c=cov_fg/var_f/var_g;
			if (c>cpre) {
				cpre=c;
				*(disparityPtr+i)=k;
			}
		}
	}
}