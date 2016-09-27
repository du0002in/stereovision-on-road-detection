#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *edgePtr, *nodePtr, *cnPtr, *onPtr, *xPtr, *yPtr, *tolPtr;
	double *out_cnPtr, *out_onPtr;
	int nc,n,k;

	edgePtr=mxGetPr(prhs[0]);
	nc=mxGetM(prhs[0]);
	nodePtr=mxGetPr(prhs[1]);
	cnPtr=mxGetPr(prhs[2]);
	onPtr=mxGetPr(prhs[3]);
	xPtr=mxGetPr(prhs[4]);
	yPtr=mxGetPr(prhs[5]);
	n=mxGetM(prhs[5]);
	tolPtr=mxGetPr(prhs[6]);

	plhs[0]=mxCreateDoubleMatrix(n, 1, mxREAL);
	plhs[1]=mxCreateDoubleMatrix(n, 1, mxREAL);
	out_cnPtr=mxGetPr(plhs[0]);
	out_onPtr=mxGetPr(plhs[1]);

	for (k=0;k<n;k++) {
		*(out_cnPtr+k)=*(cnPtr+k);
		*(out_onPtr+k)=*(onPtr+k);
	}

	for (k=0;k<nc;k++)
	{
		int n1=(int)(*(edgePtr+k))-1;
		int n2=(int)(*(edgePtr+nc+k))-1;

		int y1=(int)(*(nodePtr+n1+nc));
		int y2=(int)(*(nodePtr+n2+nc));
		int x1,x2,xmin,xmax;
		if (y1<y2) {
			x1=(int)(*(nodePtr+n1));
			x2=(int)(*(nodePtr+n2));
		}
		else {
			int yt=y1;
			y1=y2;
			y2=yt;
			x1=(int)(*(nodePtr+n2));
			x2=(int)(*(nodePtr+n1));
		}
		if (x1>x2) {
			xmin=x2;
			xmax=x1;
		}
		else {
			xmin=x1;
			xmax=x2;
		}

		int start;
		if (*yPtr>=y1)
			start=0;
		else if (*(yPtr+n-1)<y1)
			start=n;
		else {
			int lower=0;
			int upper=n-1;
			for (int j=0; j<n; j++) {
				start= (int)(0.5*(lower+upper+0.5));
				if (*(yPtr+start)<y1)
					lower=start;
				else if (*(yPtr+start-1)<y1)
					break;
				else
					upper=start;
			}
		}

		for (int j=start;j<n;j++) {
			int Y=*(yPtr+j);
			if (Y<y2) {
				int X=*(xPtr+j);
				if (X>=xmin) {
					if (X<=xmax) {
						*(out_onPtr+j)=*(out_onPtr+j) || (abs((y2-Y)*(x1-X)-(y1-Y)*(x2-X))<*tolPtr);
						if ((Y<y2) && ((y2-y1)*(X-x1)<(Y-y1)*(x2-x1)))
							*(out_cnPtr+j)=!*(out_cnPtr+j);
					}
				}
				else if (Y<y2) {
					*(out_cnPtr+j)=!*(out_cnPtr+j);
				}
			}
			else
				break;
		}
	}


}