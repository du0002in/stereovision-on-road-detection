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
/*
void GetModel(double *v,double *u, double *cos_Phai_H, double *D_pre, double *A, double *B, double *C, double *sign, double *L1L2)
{
 double H[4], rhs[2], inv_H[4];
 H[0]=*sign*(*cos_Phai_H*v[0]);	H[2]=*sign*(-*cos_Phai_H*(*D_pre));	rhs[0]=u[0]-*A/(v[0]-*D_pre)-*B*v[0]-*C;
 H[1]=*sign*(*cos_Phai_H*v[1]);	H[3]=*sign*(-*cos_Phai_H*(*D_pre)); rhs[1]=u[1]-*A/(v[1]-*D_pre)-*B*v[1]-*C;
 double det=H[0]*H[3]-H[2]*H[1];
 double inv_dtm=((det==0)?0:1.0/det);
 inv_H[0]=inv_dtm*H[3];		inv_H[2]=-inv_dtm*H[2];
 inv_H[1]=-inv_dtm*H[1];	inv_H[3]=inv_dtm*H[0];

 L1L2[0]=inv_H[0]*rhs[0]+inv_H[2]*rhs[1];
 L1L2[1]=inv_H[1]*rhs[0]+inv_H[3]*rhs[1];
}*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *leftrightlane_ptr, *Datapointl, *Datapointr, *FinalLinel, *FinalLiner, *hypo_ptr, *para_ptr, *D_pre_ptr, *H_pre_ptr, *Phai_ptr;
    double *cc_yPtr;
	size_t ncolsl, ncolsr;
    //size_t Randpoint[4];
    int j;
	int tcon_pre=0;
	double sign=1.0, A, Bli, Cli, Bri, Cri;
	
	srand(time(NULL));
	leftrightlane_ptr=mxGetPr(prhs[0]);
	if (*leftrightlane_ptr==0)
		sign=1.0;
	else sign=-1.0;

    Datapointl=mxGetPr(prhs[1]);
    ncolsl=mxGetN(prhs[1]);
    Datapointr=mxGetPr(prhs[2]);
    ncolsr=mxGetN(prhs[2]);

	para_ptr=mxGetPr(prhs[3]);
	A=*para_ptr; Bli=*(para_ptr+1); Cli=*(para_ptr+2);
	Bri=*(para_ptr+3); Cri=*(para_ptr+4);

	D_pre_ptr=mxGetPr(prhs[4]);
	Phai_ptr=mxGetPr(prhs[5]);
	H_pre_ptr=mxGetPr(prhs[6]);

	cc_yPtr=mxGetPr(prhs[7]);
	double top_row=*cc_yPtr-40;

    plhs[0]=mxCreateDoubleMatrix(1,ncolsl,mxREAL);
    FinalLinel=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(1,ncolsr,mxREAL);
    FinalLiner=mxGetPr(plhs[1]);
	plhs[2]=mxCreateDoubleMatrix(5,2,mxREAL);
    hypo_ptr=mxGetPr(plhs[2]);
	
	// C_limit is related to distance btw left and right camera 0.2566. Focal=852. 1.2 is a buffering factor.
	// double C_limit=1/sqrt(1+(*D_pre_ptr/852)*(*D_pre_ptr/852))*.2566/(*H_pre_ptr)*1.2;
	double cos_Phai_H=cos(*Phai_ptr)/(*H_pre_ptr);
	double cos_Phai_H_D=cos_Phai_H*(*D_pre_ptr);
	double ymax=(top_row<*D_pre_ptr?top_row:(*D_pre_ptr-10));
	double ratiol2r=(double)ncolsl/(ncolsr+ncolsl);
	#pragma omp parallel for num_threads(4)
    for(j=0;j<600;j++)
    {
		int tcon=0; int tcon1=0;
		double v[2], u[2];
		int rand1=rand() % ncolsl;
		*v=*(Datapointl+2*rand1);
		*u=*(Datapointl+2*rand1+1);
		v[1]=v[0];
		while (v[0]==v[1]) {
			rand1=rand() % ncolsr;
			v[1]=*(Datapointr+2*rand1);
			u[1]=*(Datapointr+2*rand1+1);
		}
		double L;
		L=(u[0]-u[1]-(A*(1/(v[0]-*D_pre_ptr)-1/(v[1]-*D_pre_ptr))+Bli*v[0]-Bri*v[1]+Cli-Cri))/(v[0]-v[1])/cos_Phai_H/sign;
		double Bli2=Bli+cos_Phai_H*L*sign;
		double Cli2=u[0]-(A/(v[0]-*D_pre_ptr)+Bli2*v[0]);
		double Eli2=Cli2-Bli2**D_pre_ptr, Fli2=A-Cli2**D_pre_ptr;
		double Bri2=Bri+cos_Phai_H*L*sign;
		double Cri2=u[1]-(A/(v[1]-*D_pre_ptr)+Bri2*v[1]);
		double Eri2=Cri2-Bri2**D_pre_ptr, Fri2=A-Cri2**D_pre_ptr;
		double Cp1, Cp2, Cp3, x, y, x1,t_dist,dist;
		double *line_ptr_l=new double[ncolsl];
		double *line_ptr_r=new double[ncolsr];
		if (L>2.2 && L<4.7) {
			
			for (int i=0; i<ncolsl; i++) {
				x=*(Datapointl+2*i+1);
				y=*(Datapointl+2*i);
				if (tcon+ncolsl+ncolsr-i>tcon_pre) {
					if (y>=*D_pre_ptr) {
						*(line_ptr_l+i)=0;
						continue;
					}
					x1=A/(y-*D_pre_ptr)+Bli2*y+Cli2;
					if (abs(x-x1)>10) {
						*(line_ptr_l+i)=0;
						continue;
					}
					Cp1=(*D_pre_ptr-y)*0.5;
					Cp2=(Eli2-x1)*0.5+Bli2*y;
					Cp3=Fli2+(*D_pre_ptr*x1+Eli2*y)*0.5;
					dist=pow(x*Cp1+y*Cp2+Cp3,2)/4.0/(Cp1*Cp1+Cp2*Cp2);
					t_dist=0.5-(y-top_row)/220;
					if (dist<t_dist){
						tcon++;
						*(line_ptr_l+i)=1;}
					else
						*(line_ptr_l+i)=0;
				}
				else
					break;
			}
			for (int i=0; i<ncolsr; i++) {
				x=*(Datapointr+2*i+1);
				y=*(Datapointr+2*i);
				if (tcon+ncolsr-i>tcon_pre) {
					if (y>=*D_pre_ptr) {
						*(line_ptr_r+i)=0;
						continue;
					}
					x1=A/(y-*D_pre_ptr)+Bri2*y+Cri2;
					if (abs(x-x1)>10) {
						*(line_ptr_r+i)=0;
						continue;
					}
					Cp1=(*D_pre_ptr-y)*0.5;
					Cp2=(Eri2-x1)*0.5+Bri2*y;
					Cp3=Fri2+(*D_pre_ptr*x1+Eri2*y)*0.5;
					dist=pow(x*Cp1+y*Cp2+Cp3,2)/4.0/(Cp1*Cp1+Cp2*Cp2);
					t_dist=0.5-(y-top_row)/220;
					if (dist<t_dist){
						tcon++;
						*(line_ptr_r+i)=1;}
					else
						*(line_ptr_r+i)=0;
				}
				else
					break;
			}
		}
		#pragma omp critical
		{
		if (tcon>tcon_pre)
		{
			tcon_pre=tcon;
			//*total_num_ptr=tcon_pre;
			*(hypo_ptr)=A;	*(hypo_ptr+1)=Bli2;	*(hypo_ptr+2)=Cli2;	*(hypo_ptr+3)=Bri2;	*(hypo_ptr+4)=Cri2;
			*(hypo_ptr+5)=0;	*(hypo_ptr+6)=0;	*(hypo_ptr+7)=0;		*(hypo_ptr+8)=0;		*(hypo_ptr+9)=0;
			for (int i=0;i<ncolsl;i++)
				*(FinalLinel+i)=*(line_ptr_l+i);
			for (int i=0;i<ncolsr;i++)
				*(FinalLiner+i)=*(line_ptr_r+i);
		}
		}
		delete [] line_ptr_l;
		delete [] line_ptr_r;
	}
}