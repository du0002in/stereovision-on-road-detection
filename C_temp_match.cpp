#include "mex.h"
#include <omp.h>
#include <math.h>

void temp_match(double *Img, int Ir,int Ic,int totalIr, int totalIc,
	double *temp,int tr, int tc, int totaltr, int totaltc, double *mis_rate)
{
	double sum_count=0, miss_count=0;
	for (int i=0; i<totaltc; i++) {
		if ((Ic-tc+i)<0) continue;
		else if ((Ic-tc+i)>=totalIc) break;
		for (int j=0; j<totaltr; j++) {
			if ((int)(*(temp+i*totaltr+j))==-1 || 
				((int)(*(Img+(Ic-tc+i)*totalIr+(Ir-tr+j))))==-1) continue;
			miss_count=miss_count+(((int)(*(Img+(Ic-tc+i)*totalIr+(Ir-tr+j)))==
				(int)(*(temp+i*totaltr+j)))?0:1);
			sum_count++;
		}
	}
	*mis_rate=miss_count/sum_count;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *ImgPtr, *indPtr, *ref_rPtr, *ref_cPtr, *char_widthPtr, *thetaPtr;
	const mxArray *templatePtr;
	double *outPtr;
	mwSize mrows, temp_num, indlen;
    mwSize ncols, tempcols;
	int i;

	ImgPtr=mxGetPr(prhs[0]);
	mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
	indPtr=mxGetPr(prhs[1]);
	indlen=mxGetM(prhs[1]);
	templatePtr=prhs[2];
	temp_num=mxGetNumberOfElements(prhs[2]);
	ref_rPtr=mxGetPr(prhs[3]);
	ref_cPtr=mxGetPr(prhs[4]);
	char_widthPtr=mxGetPr(prhs[5]);
	thetaPtr=mxGetPr(prhs[6]);

	plhs[0]=mxCreateDoubleMatrix(1,2, mxREAL);
	outPtr=mxGetPr(plhs[0]); *outPtr=1;

	//#pragma omp parallel for num_threads(4)
	for (i=0; i<indlen; i++)
	{
		int ind=(int)(*(indPtr+i));
		int row=ind%mrows;
		int col=ind/mrows;
		if (row==50 && col==66)
			int a=0;
		for (int j=0; j<temp_num; j++)
		{
			if (j==0 || j==3 || j==8 || j==12 || j==15 || j==16 || j==19
				 || j==21 || j==22 || j==23 || j==24)
				 continue;
			if (j==13)
				double a=0;
			const mxArray *cellArray;
			double *p;
			cellArray=mxGetCell(templatePtr,j);
			p=mxGetPr(cellArray);
			int temp_rows=mxGetM(cellArray);
			int temp_cols=mxGetN(cellArray);
			double mis_rate;
			temp_match(ImgPtr,row,col,mrows,ncols,p,
				(int)(*(ref_rPtr+j)),(int)(*(ref_cPtr+j)),temp_rows,temp_cols,&mis_rate);
			if (mis_rate<0.3) {
				switch (j)
				{
				case 13:
					double overall_miss_rate[4];
					overall_miss_rate[1]=mis_rate;
					int ori_col=col, ori_row=row;
					for (int k=14; k<16; k++) {
						double char_wid=*(char_widthPtr+k-1);
						int seed_c=(int)(col+char_wid*cos(*thetaPtr)/2.0);
						int seed_r=(int)(row-char_wid*sin(*thetaPtr)/5.0);
						const mxArray *cellArray;
						double *p;
						cellArray=mxGetCell(templatePtr,k);
						p=mxGetPr(cellArray);
						int temp_rows=mxGetM(cellArray);
						int temp_cols=mxGetN(cellArray);
						double mis_rate=1.0;
						for (int m=3; m<7; m++) {
							int col_tmp=seed_c+m;
							for (int n=-1; n<2; n++) {
								int row_tmp=seed_r+n;
								double current_mis_rate;
								temp_match(ImgPtr,row_tmp,col_tmp,mrows,ncols,p,
									(int)(*(ref_rPtr+k)),(int)(*(ref_cPtr+k)),temp_rows,temp_cols,&current_mis_rate);
								if (current_mis_rate<mis_rate) {
									mis_rate=current_mis_rate;
									row=row_tmp; col=col_tmp;
								}
								if (mis_rate<0.3)
									break;
							}
							if (mis_rate<0.3)
								break;
						}
						overall_miss_rate[k-12]=mis_rate;
					}
					col=ori_col; row=ori_row;
					for (int k=12; k<13; k++) {
						double char_wid=-(*(char_widthPtr+k));
						int seed_c=(int)(col+char_wid*cos(*thetaPtr)/2.0);
						int seed_r=(int)(row-char_wid*sin(*thetaPtr)/5.0);
						const mxArray *cellArray;
						double *p;
						cellArray=mxGetCell(templatePtr,k);
						p=mxGetPr(cellArray);
						int temp_rows=mxGetM(cellArray);
						int temp_cols=mxGetN(cellArray);
						double mis_rate=1.0;
						for (int m=-3; m>-7; m--) {
							int col_tmp=seed_c+m;
							for (int n=-1; n<2; n++) {
								int row_tmp=seed_r+n;
								double current_mis_rate;
								temp_match(ImgPtr,row_tmp,col_tmp,mrows,ncols,p,
									(int)(*(ref_rPtr+k)),(int)(*(ref_cPtr+k)),temp_rows,temp_cols,&current_mis_rate);
								if (current_mis_rate<mis_rate) {
									mis_rate=current_mis_rate;
									row=row_tmp; col=col_tmp;
								}
								if (mis_rate<0.3)
									break;
							}
							if (mis_rate<0.3)
								break;
						}
						overall_miss_rate[k-12]=mis_rate;
					}
					double ave_mis_rate=0;
					for (int m=0; m<4; m++)
						ave_mis_rate=ave_mis_rate+overall_miss_rate[m];
					ave_mis_rate=ave_mis_rate/4;
					*outPtr=ave_mis_rate;
					if (ave_mis_rate<0.3) {
						*outPtr=ave_mis_rate;
						*(outPtr+1)=4;
					}


				}


			}



		}



	}
}