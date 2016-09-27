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

void GetModel(double *v,double *u, double *cos_Phai_H, double *D_pre, double *c, double *e, double *f, double *sign, double *L1L2)
{
 double H[4], rhs[2], inv_H[4];
 H[0]=*sign*(*cos_Phai_H*v[0]);	H[2]=*sign*(-*cos_Phai_H*(*D_pre));	rhs[0]=u[0]+*e*v[0]+*f*0.5;
 H[1]=*sign*(*cos_Phai_H*v[1]);	H[3]=*sign*(-*cos_Phai_H*(*D_pre)); rhs[1]=u[1]+*e*v[1]+*f*0.5;
 double det=H[0]*H[3]-H[2]*H[1];
 double inv_dtm=((det==0)?0:1.0/det);
 inv_H[0]=inv_dtm*H[3];		inv_H[2]=-inv_dtm*H[2];
 inv_H[1]=-inv_dtm*H[1];	inv_H[3]=inv_dtm*H[0];

 L1L2[0]=inv_H[0]*rhs[0]+inv_H[2]*rhs[1];
 L1L2[1]=inv_H[1]*rhs[0]+inv_H[3]*rhs[1];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *leftrightlane_ptr, *Datapointl, *Datapointr, *FinalLinel, *FinalLiner, *hypo_ptr, *para_ptr, *D_pre_ptr, *H_pre_ptr, *Phai_ptr;
    double *cc_yPtr;
	size_t ncolsl, ncolsr;
    //size_t Randpoint[4];
    int j;
	int tcon_pre=0;
	double sign=1.0, c, eli, fli, eri, fri;

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
	c=*(para_ptr+0); eli=*(para_ptr+1); fli=*(para_ptr+2);
	eri=*(para_ptr+3); fri=*(para_ptr+4);

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
    for(j=0;j<300;j++)
    {
		int tcon=0; int tcon1=0;
		double rand0=((double) rand() / (RAND_MAX));
		double v[2], u[2], L1L2[2];;
		if (rand0<ratiol2r) //select points from left image
		{
			int rand1=rand() % ncolsl;
			*v=*(Datapointl+2*rand1);
			*u=*(Datapointl+2*rand1+1);
			rand1=rand() % ncolsl;
			*(v+1)=*(Datapointl+2*rand1);
			*(u+1)=*(Datapointl+2*rand1+1);
			GetModel(v,u,&cos_Phai_H,D_pre_ptr,&c,&eli,&fli,&sign,L1L2);
		}
		else //select points from right image
		{
			int rand1=rand() % ncolsr;
			*v=*(Datapointr+2*rand1);
			*u=*(Datapointr+2*rand1+1);
			rand1=rand() % ncolsr;
			*(v+1)=*(Datapointr+2*rand1);
			*(u+1)=*(Datapointr+2*rand1+1);
			GetModel(v,u,&cos_Phai_H,D_pre_ptr,&c,&eri,&fri,&sign,L1L2);
		}

		double eli2=eli-cos_Phai_H*L1L2[0]*sign, fli2=fli+2.0*cos_Phai_H_D*L1L2[1]*sign;
		double eri2=eri-cos_Phai_H*L1L2[0]*sign, fri2=fri+2.0*cos_Phai_H_D*L1L2[1]*sign;		
        double denoleft=c*c+eli2*eli2;
		double denoright=c*c+eri2*eri2;

		double Cp1, Cp2, Cp3, x, y, x1,t_dist,dist;
		double *line_ptr_l1=new double[ncolsl];
		double *line_ptr_r1=new double[ncolsr];
		if (L1L2[0]>2.2 && L1L2[0]<4.7 && L1L2[1]>2.2 && L1L2[1]<4.7) {
			for (int i=0; i<ncolsl; i++) {
				if (tcon1+ncolsl+ncolsr-i>tcon_pre) {
					x=*(Datapointl+2*i+1);
					y=*(Datapointl+2*i);
					t_dist=1.0-(y-top_row)/440*9;
					dist=pow(c*x+eli2*y+fli2*0.5,2)/denoleft;
					if (dist<t_dist) {
						tcon1++;
						*(line_ptr_l1+i)=1;}
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
					dist=pow(c*x+eri2*y+fri2*0.5,2)/denoright;
					if (dist<t_dist) {
						tcon1++;
						*(line_ptr_r1+i)=1;}
					else
						*(line_ptr_r1+i)=0;	
				}
				else
					break;
			}
		}
		#pragma omp critical
		{
		if (tcon1>tcon_pre)
		{
			tcon_pre=tcon1;
			*(hypo_ptr)=0;		*(hypo_ptr+1)=0;		*(hypo_ptr+2)=0;		*(hypo_ptr+3)=0;		*(hypo_ptr+4)=0;
			*(hypo_ptr+5)=c;	*(hypo_ptr+6)=eli2;	*(hypo_ptr+7)=fli2;	*(hypo_ptr+8)=eri2;	*(hypo_ptr+9)=fri2;
			for (int i=0;i<ncolsl;i++)
				*(FinalLinel+i)=*(line_ptr_l1+i);
			for (int i=0;i<ncolsr;i++)
				*(FinalLiner+i)=*(line_ptr_r1+i);
		}
		}
		delete [] line_ptr_l1;
		delete [] line_ptr_r1;
	}
}