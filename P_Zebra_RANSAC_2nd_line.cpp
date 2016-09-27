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

void inverse3x3(double *a, double *in)
{
	double det;
	det=a[0]*(a[4]*a[8]-a[5]*a[7])-a[3]*(a[1]*a[8]-a[7]*a[2])+a[6]*(a[1]*a[5]-a[4]*a[2]);
	det=((det==0)?0:1/det);
	in[0]=(a[4]*a[8]-a[5]*a[7])*det;
	in[1]=-(a[1]*a[8]-a[7]*a[2])*det;
	in[2]=(a[1]*a[5]-a[2]*a[4])*det;
	in[3]=-(a[3]*a[8]-a[6]*a[5])*det;
	in[4]=(a[0]*a[8]-a[6]*a[2])*det;
	in[5]=-(a[0]*a[5]-a[2]*a[3])*det;
	in[6]=(a[3]*a[7]-a[6]*a[4])*det;
	in[7]=-(a[0]*a[7]-a[1]*a[6])*det;
	in[8]=(a[0]*a[4]-a[1]*a[3])*det;
}

void GetModel(double *v,double *u, double *D_pre, double *BCD_L)
{
 double L[9], dU[3], inv_L[9];
 L[0]=v[0]-v[2];	L[3]=v[2];	L[6]=-1;	dU[0]=u[0]-u[2];
 L[1]=v[1]-v[3];	L[4]=v[3];	L[7]=-1;	dU[1]=u[1]-u[3];
 L[2]=v[0]-v[3];	L[5]=v[3];	L[8]=-1;	dU[2]=u[0]-u[3];
 inverse3x3(L,inv_L);
 
 BCD_L[0]=inv_L[0]*dU[0]+inv_L[3]*dU[1]+inv_L[6]*dU[2];
 BCD_L[1]=inv_L[1]*dU[0]+inv_L[4]*dU[1]+inv_L[7]*dU[2];
 BCD_L[2]=inv_L[2]*dU[0]+inv_L[5]*dU[1]+inv_L[8]*dU[2];
 BCD_L[2]=BCD_L[2]/BCD_L[1];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *Datapointl, *Datapointr, *FinalLinel, *FinalLiner, *SubDatapointl, *SubDatapointr, *hypo_ptr, *D_pre_ptr, *H_pre_ptr, *num_ite_ptr;
    double *focalPtr, *dist_left_rightPtr, *cc_yPtr;
	size_t ncolsl, ncolsr, sub_ncolsl, sub_ncolsr;
    //size_t Randpoint[4];
    int j, num_ite;
	int tcon_pre=0, total_dist_pre=0;
	double top_row=0;

	srand(time(NULL));
    Datapointl=mxGetPr(prhs[0]);
	SubDatapointl=mxGetPr(prhs[1]);
    ncolsl=mxGetN(prhs[0]);
	sub_ncolsl=mxGetN(prhs[1]);
	
    Datapointr=mxGetPr(prhs[2]);
	SubDatapointr=mxGetPr(prhs[3]);
    ncolsr=mxGetN(prhs[2]);
	sub_ncolsr=mxGetN(prhs[3]);

	D_pre_ptr=mxGetPr(prhs[4]);
	H_pre_ptr=mxGetPr(prhs[5]);
	

	num_ite_ptr=mxGetPr(prhs[6]);
	num_ite=int(*num_ite_ptr);

	focalPtr=mxGetPr(prhs[7]);
	dist_left_rightPtr=mxGetPr(prhs[8]);
	cc_yPtr=mxGetPr(prhs[9]);
	double phi=atan(*D_pre_ptr/(*focalPtr));

	double *Rnd_l=mxGetPr(prhs[10]);
	double *Rnd_r=mxGetPr(prhs[11]);
	double *theta_l=mxGetPr(prhs[12]);
	double *theta_r=mxGetPr(prhs[13]);


	top_row=*cc_yPtr-40;

    plhs[0]=mxCreateDoubleMatrix(1,ncolsl,mxREAL);
    FinalLinel=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(1,ncolsr,mxREAL);
    FinalLiner=mxGetPr(plhs[1]);
	plhs[2]=mxCreateDoubleMatrix(5,1,mxREAL);
    hypo_ptr=mxGetPr(plhs[2]);
	plhs[3]=mxCreateDoubleMatrix(1,1,mxREAL);
    double *break_location_ptr=mxGetPr(plhs[3]);
	// C_limit is related to distance btw left and right camera 0.2566. Focal=852. 1.2 is a buffering factor.
	double C_limit=1/sqrt(1+(*D_pre_ptr/(*focalPtr))*(*D_pre_ptr/(*focalPtr)))*(*dist_left_rightPtr)/(*H_pre_ptr)*1.5;
	double ymax=(top_row<*D_pre_ptr?top_row:(*D_pre_ptr-10));
	#pragma omp parallel for num_threads(4)
    for(j=0;j<num_ite;j++)
    {
		int tcon=0; int tcon1=0; double total_dist=0;
		int rand1=rand() % sub_ncolsl;
		double v[4],u[4];
		*v=*(SubDatapointl+2*rand1);
		*u=*(SubDatapointl+2*rand1+1);
		//rand1=rand() % sub_ncolsl;
		*(v+1)=*(Rnd_l+2*j);
		*(u+1)=*(Rnd_l+2*j+1);
		rand1=rand() % sub_ncolsr;
		*(v+2)=*(SubDatapointr+2*rand1);
		*(u+2)=*(SubDatapointr+2*rand1+1);
		//rand1=rand() % sub_ncolsr;
		*(v+3)=*(Rnd_r+2*j);
		*(u+3)=*(Rnd_r+2*j+1);
		
		double bcd_l[3];
		GetModel(v,u,D_pre_ptr,bcd_l);
		//left and right Line
		double c=1, eleft=-bcd_l[0], fpleft=-eleft*v[0]-u[0], fleft=fpleft*2;
		double eright=-(bcd_l[0]-bcd_l[1]), fpright=fpleft-bcd_l[1]*bcd_l[2], fright=fpright*2;
        double denoleft=c*c+eleft*eleft;
		double denoright=c*c+eright*eright;

		double Cp1, Cp2, Cp3, x, y, x1,t_dist,dist;
		double *line_ptr_l1=new double[ncolsl];
		double *line_ptr_r1=new double[ncolsr];

		double theta_l_p=atan((-eleft*sin(phi)+fpleft/(*focalPtr)*cos(phi))/(1-2*cos(phi)*cos(phi)));
		double theta_r_p=atan((-eright*sin(phi)+fpright/(*focalPtr)*cos(phi))/(1-2*cos(phi)*cos(phi)));
		bool valid=(abs(theta_l_p-*theta_l)>0.12) && (abs(theta_r_p-*theta_r)>0.12);

		if (bcd_l[1]<C_limit && (eleft-eright)/c<0 && ((eleft-eright)/c*top_row+(fpleft-fpright)/c)>0 && valid) {
			for (int i=0; i<ncolsl; i++) {
				if (tcon1+ncolsl+ncolsr-i>tcon_pre) {
					x=*(Datapointl+2*i+1);
					y=*(Datapointl+2*i);
					t_dist=1.0-(y-top_row)/440*9;
					dist=pow(c*x+eleft*y+fpleft,2)/denoleft;
					if (dist<t_dist) {
						tcon1++;
						*(line_ptr_l1+i)=1;
						total_dist=total_dist+dist;
					}
					else
						*(line_ptr_l1+i)=0;	
				}
				else
					break;
			}
			for (int i=0; i<ncolsr; i++) {
				if (tcon1+ncolsr-i>tcon_pre) {
					x=*(Datapointr+2*i+1);
					y=*(Datapointr+2*i);
					t_dist=1.0-(y-top_row)/440*9;
					dist=pow(c*x+eright*y+fpright,2)/denoright;
					if (dist<t_dist) {
						tcon1++;
						*(line_ptr_r1+i)=1;
						total_dist=total_dist+dist;
					}
					else
						*(line_ptr_r1+i)=0;	
				}
				else
					break;
			}
		}
		#pragma omp critical
		{
			if ((tcon1>tcon_pre) || (tcon1==tcon_pre && total_dist<total_dist_pre))
			{
				tcon_pre=tcon1;
				total_dist_pre=total_dist;
				*(hypo_ptr+0)=c;	*(hypo_ptr+1)=eleft;	*(hypo_ptr+2)=fleft;
				*(hypo_ptr+3)=eright;	*(hypo_ptr+4)=fright;
				for (int i=0;i<ncolsl;i++)
					*(FinalLinel+i)=*(line_ptr_l1+i);
				for (int i=0;i<ncolsr;i++)
					*(FinalLiner+i)=*(line_ptr_r1+i);
				break_location_ptr[0]=*(Rnd_l+2*j);
			}
		}
		delete [] line_ptr_l1;
		delete [] line_ptr_r1;
	}
}