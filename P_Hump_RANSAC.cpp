// mex file to run RANSAC to improve speed
// input: 2xn matrix with 1st row being x and 2nd row being y. n refers to the number of points
// ouput: 1xn matrix (or row vector) to indicate whether the corresponding points are inliers (1) or not (0)
// fitting left and right img simoutanously
// Cordinate must be in camera coordinate form instead of img coordinate!! Follow exactly as notes in 3D Vision.
#include "mex.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <vector>

void GetModel(double *v,double *u, double *ab)
{
 ab[0]=(u[0]-u[1])/(v[0]-v[1]);
 ab[1]=u[0]-ab[0]*v[0];
}


void Random_2_Unique(int *arr, int interval)
{
	arr[0]=rand() % interval;
	arr[1]=rand() % interval;
	while (arr[1]==arr[0])
		arr[1]=rand() % interval;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *Datapointl, *Datapointr, *FinalLinel, *FinalLiner, *FinalLinel_index,
		*SubDatapointl, *SubDatapointr, *FinalTheta,*Finalxc, *Phi_pre_ptr, *H_pre_ptr, *num_ite_ptr;
    double *focalPtr, *dist_left_rightPtr, *cc_yPtr, *cc_xPtr, *y_int_ptr, *theta_pre_ptr;
	size_t ncolsl, ncolsr, sub_ncolsl, sub_ncolsr;
    //size_t Randpoint[4];
    int j, num_ite;
	int tcon_pre=0;
	int line_num=13; //this has to be an odd number

	srand(time(NULL));
    Datapointl=mxGetPr(prhs[0]);
	SubDatapointl=mxGetPr(prhs[1]);
    ncolsl=mxGetN(prhs[0]);
	sub_ncolsl=mxGetN(prhs[1]);
	
	Phi_pre_ptr=mxGetPr(prhs[2]);
	H_pre_ptr=mxGetPr(prhs[3]);

	num_ite_ptr=mxGetPr(prhs[4]);
	num_ite=int(*num_ite_ptr);

	focalPtr=mxGetPr(prhs[5]);
	theta_pre_ptr=mxGetPr(prhs[6]);

	cc_yPtr=mxGetPr(prhs[7]);
	cc_xPtr=mxGetPr(prhs[8]);
	double top_row=*cc_yPtr;
	double btm_row=-480+*cc_yPtr;
	double right_clm=*cc_xPtr;
	double left_clm=-640+*cc_xPtr;

    plhs[0]=mxCreateDoubleMatrix(1,ncolsl,mxREAL);
    FinalLinel=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(1,ncolsl,mxREAL);
    FinalLinel_index=mxGetPr(plhs[1]);
	plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
    FinalTheta=mxGetPr(plhs[2]);
	plhs[3]=mxCreateDoubleMatrix(1,line_num,mxREAL);
    Finalxc=mxGetPr(plhs[3]);

	
	for (int i=0;i<ncolsl;i++) {
		*(FinalLinel+i)=0;
		*(FinalLinel_index+i)=line_num+1;
	}
	*FinalTheta=0;
	for (int i=0;i<line_num;i++)
		*(Finalxc+i)=0;


	//center of theta lb and ub is theta_pre-pi/4 coz the hump lines are 45 slanted wrt road tangent
	double theta_lb=*theta_pre_ptr-3.1415926/8.0*3.0;
	double theta_ub=*theta_pre_ptr-3.1415926/8.0*1.0;
	double sin_phi=sin(*Phi_pre_ptr);
	double cos_phi=cos(*Phi_pre_ptr);
	double sin_phi_to_H=sin_phi/(*H_pre_ptr);
	double cos_phi_to_H=cos_phi/(*H_pre_ptr);


	//#pragma omp parallel for num_threads(4)
    for(j=0;j<num_ite;j++)
    {
		int tcon1=0;
		//int rand1=rand() % sub_ncolsl;
		// choose 2 pt from left img
		double vl[2],ul[2];
		int rand_num[2];
		vl[0]=0;vl[1]=0;
		while (vl[0]==vl[1]) {
			Random_2_Unique(rand_num, sub_ncolsl);
			//rand_num[0]=280; rand_num[1]=377;
			*vl=*(SubDatapointl+2*rand_num[0]);
			*ul=*(SubDatapointl+2*rand_num[0]+1);
			*(vl+1)=*(SubDatapointl+2*rand_num[1]);
			*(ul+1)=*(SubDatapointl+2*rand_num[1]+1);
		}
		double ab[2],theta,xc[13];
		GetModel(vl,ul,ab);
		theta=atan(sin_phi*ab[0]+cos_phi*ab[1]/(*focalPtr));
		double *line_ptr_l1=new double[ncolsl];
		double *line_index=new double[ncolsl];
		if (theta>theta_lb && theta<theta_ub) {
			xc[line_num/2]=(sin_phi*ab[1]/(*focalPtr)-ab[0]*cos_phi)*(*H_pre_ptr)*cos(theta);
			for (int i=0;i<line_num/2; i++) {
				xc[line_num/2-1-i]=xc[line_num/2-i]+0.38;
				xc[line_num/2+1+i]=xc[line_num/2+i]-0.38;
			}
			/*xc[2]=xc[3]+0.38;
			xc[1]=xc[2]+0.38;
			xc[0]=xc[1]+0.38;
			xc[4]=xc[3]-0.38;
			xc[5]=xc[4]-0.38;
			xc[6]=xc[5]-0.38;*/
			double a[13],b[13],deno[13],t_dist=5.0,dist,x,y;
			int valid_counter=0, valid_index_tmp[13];
			for (int i=0;i<line_num;i++) {
				a[i]=sin_phi*tan(theta)-cos_phi_to_H*xc[i]/cos(theta);
				b[i]=(cos_phi*tan(theta)+sin_phi_to_H*xc[i]/cos(theta))*(*focalPtr);
				deno[i]=pow(a[i]*a[i]+1.0,0.5);
				double u_top=a[i]*top_row+b[i];
				if (u_top>left_clm && u_top<right_clm) {
					valid_index_tmp[valid_counter]=i;
					valid_counter++;
					continue;
				}
				double u_btm=a[i]*btm_row+b[i];
				if (u_btm>left_clm && u_btm<right_clm) {
					valid_index_tmp[valid_counter]=i;
					valid_counter++;
					continue;
				}
				double v_left=(left_clm-b[i])/a[i];
				if (v_left<top_row && v_left>btm_row) {
					valid_index_tmp[valid_counter]=i;
					valid_counter++;
					continue;
				}
				double v_right=(right_clm-b[i])/a[i];
				if (v_right<top_row && v_right>btm_row) {
					valid_index_tmp[valid_counter]=i;
					valid_counter++;
					continue;
				}
				a[i]=0; b[i]=0; xc[i]=0;
			}
			
			for (int i=0; i<ncolsl; i++) {
				if (tcon1+ncolsl-i>tcon_pre) {
					x=*(Datapointl+2*i+1);
					y=*(Datapointl+2*i);
					t_dist=5.0-(y-top_row)/480*5;
					for (int k=0; k<valid_counter; k++) {
						//dist=abs(a[valid_index_tmp[k]]*y+b[valid_index_tmp[k]]-x)/deno[valid_index_tmp[k]];
						dist=abs(a[valid_index_tmp[k]]*y+b[valid_index_tmp[k]]-x);
						if (dist<t_dist) {
							tcon1++;
							*(line_ptr_l1+i)=1;
							*(line_index+i)=valid_index_tmp[k];
							break;
						}
						else {
							*(line_ptr_l1+i)=0;
							*(line_index+i)=line_num+1;
						}
					}
				}
				else
					break;
			}
		}

		//#pragma omp critical
		{
			if (tcon1>tcon_pre)
			{			
				tcon_pre=tcon1;
				*FinalTheta=theta;
				for (int i=0;i<line_num;i++)
					*(Finalxc+i)=xc[i];
				for (int i=0;i<ncolsl;i++) {
					*(FinalLinel+i)=*(line_ptr_l1+i);
					*(FinalLinel_index+i)=*(line_index+i);
				}
			}
		}						
		delete [] line_ptr_l1;
		delete [] line_index;	
	}
}