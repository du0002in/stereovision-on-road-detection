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

void GetModel(double *v,double *u, double *D_pre, double *BAC_H, double *BCD_L)
{
 double H[9], L[9], dU[3], inv_H[9], inv_L[9];
 H[0]=v[0]-v[2];	H[3]=1/(v[0]-*D_pre)-1/(v[2]-*D_pre);	H[6]=v[2]-*D_pre;	dU[0]=u[0]-u[2];
 H[1]=v[1]-v[3];	H[4]=1/(v[1]-*D_pre)-1/(v[3]-*D_pre);	H[7]=v[3]-*D_pre;	dU[1]=u[1]-u[3];
 H[2]=v[0]-v[3];	H[5]=1/(v[0]-*D_pre)-1/(v[3]-*D_pre);	H[8]=v[3]-*D_pre;	dU[2]=u[0]-u[3];
 L[0]=H[0];	L[3]=v[2];	L[6]=-1;
 L[1]=H[1];	L[4]=v[3];	L[7]=-1;
 L[2]=H[2];	L[5]=v[3];	L[8]=-1;
 inverse3x3(H,inv_H);
 inverse3x3(L,inv_L);
 BAC_H[0]=inv_H[0]*dU[0]+inv_H[3]*dU[1]+inv_H[6]*dU[2];
 BAC_H[1]=inv_H[1]*dU[0]+inv_H[4]*dU[1]+inv_H[7]*dU[2];
 BAC_H[2]=inv_H[2]*dU[0]+inv_H[5]*dU[1]+inv_H[8]*dU[2];
 
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
	int tcon_pre=0;
	double top_row=0;
	double prob=0.0;


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
	double *pre_kPtr=mxGetPr(prhs[10]);
	double v_btm=-480+*cc_yPtr;
	double dist_left_right_x_focal=*focalPtr*(*dist_left_rightPtr);

	top_row=*cc_yPtr-40;
	double Phai=atan(*D_pre_ptr/(*focalPtr));
	double k_const=4*pow(cos(Phai),3)/(*H_pre_ptr*(*focalPtr)*(*focalPtr));
	//double miu=-0.00032972;
	//double sigma=0.0017556;
	double miu=0.000014854;
	double sigma=0.04;

    plhs[0]=mxCreateDoubleMatrix(1,ncolsl,mxREAL);
    FinalLinel=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(1,ncolsr,mxREAL);
    FinalLiner=mxGetPr(plhs[1]);
	plhs[2]=mxCreateDoubleMatrix(5,2,mxREAL);
    hypo_ptr=mxGetPr(plhs[2]);
	
	// C_limit is related to distance btw left and right camera 0.2566. Focal=852. 1.2 is a buffering factor.
	double C_limit=1/sqrt(1+(*D_pre_ptr/(*focalPtr))*(*D_pre_ptr/(*focalPtr)))*(*dist_left_rightPtr)/(*H_pre_ptr)*1.5;
	double ymax=(top_row<*D_pre_ptr?top_row:(*D_pre_ptr-10));
	#pragma omp parallel for num_threads(4)
    for(j=0;j<num_ite;j++)
    {
		int tcon=0; int tcon1=0;
		int rand1=rand() % sub_ncolsl;
		double v[4],u[4];
		*v=*(SubDatapointl+2*rand1);
		*u=*(SubDatapointl+2*rand1+1);
		rand1=rand() % sub_ncolsl;
		*(v+1)=*(SubDatapointl+2*rand1);
		*(u+1)=*(SubDatapointl+2*rand1+1);
		rand1=rand() % sub_ncolsr;
		*(v+2)=*(SubDatapointr+2*rand1);
		*(u+2)=*(SubDatapointr+2*rand1+1);
		rand1=rand() % sub_ncolsr;
		*(v+3)=*(SubDatapointr+2*rand1);
		*(u+3)=*(SubDatapointr+2*rand1+1);
		
		//*v=-27; *(v+1)=67; *(v+2)=-9; *(v+3)=69;
		//*u=-254; *(u+1)=-83; *(u+2)=-169; *(u+3)=-38;
		
		double BAC_H[3], bcd_l[3];
		GetModel(v,u,D_pre_ptr,BAC_H,bcd_l);
		//left and right Hypobola
		double Bleft=BAC_H[0], A=BAC_H[1], C=BAC_H[2];
		double K_H=A*k_const;
		//double prob_H=exp(-((K_H-*pre_kPtr)-miu)/sigma)/sigma/pow((1+exp(-((K_H-*pre_kPtr)-miu)/sigma)),2);
		double prob_H=1/2.506628/sigma*exp(-pow(((K_H-*pre_kPtr)-miu),2)/2.0/pow(sigma,2));
		prob_H=(prob_H>5.0)?prob_H:5.0;
		double F_theta_cos=u[0]-A/(v[0]-*D_pre_ptr)-Bleft*v[0]+Bleft**D_pre_ptr;
		double Cleft=-Bleft**D_pre_ptr+F_theta_cos; //ul=A/(vl-D)+Bleft*vl+Cleft;
		double Eleft=Cleft-Bleft**D_pre_ptr, Fleft=A-Cleft**D_pre_ptr;
		double Bright=Bleft-C;
		double Cright=-(Bleft-C)**D_pre_ptr+F_theta_cos; //ur=A/(vr-D)+Bright*vl+Cright;
		double Eright=Cright-Bright**D_pre_ptr, Fright=A-Cright**D_pre_ptr;
		double Zc_btm_H=dist_left_right_x_focal/((A/(v_btm-*D_pre_ptr)+Bright*v_btm+Cright)-
			(A/(v_btm-*D_pre_ptr)+Bleft*v_btm+Cleft));

		//left and right Line
		double c=1, eleft=-bcd_l[0], fpleft=-eleft*v[0]-u[0], fleft=fpleft*2;
		double eright=-(bcd_l[0]-bcd_l[1]), fpright=fpleft-bcd_l[1]*bcd_l[2], fright=fpright*2;
		double zc_btm_l=dist_left_right_x_focal/((-fpright-eright*v_btm)/c-(-fpleft-eleft*v_btm)/c);
		double k_l=0.0;
		//double prob_l=exp(-((k_l-*pre_kPtr)-miu)/sigma)/sigma/pow((1+exp(-((k_l-*pre_kPtr)-miu)/sigma)),2);
		double prob_l=1/2.506628/sigma*exp(-pow(((k_l-*pre_kPtr)-miu),2)/2.0/pow(sigma,2));
		prob_l=(prob_l>5.0)?prob_l:5.0;
        double denoleft=c*c+eleft*eleft;
		double denoright=c*c+eright*eright;

		double Cp1, Cp2, Cp3, x, y, x1,t_dist,dist;
		double *line_ptr_l1=new double[ncolsl];
		double *line_ptr_l=new double[ncolsl];
		double *line_ptr_r1=new double[ncolsr];
		double *line_ptr_r=new double[ncolsr];
		if (BAC_H[2]<C_limit && (Bright-Bleft)<0 && ((Bright-Bleft)*ymax+Cright-Cleft)>0
			&& Zc_btm_H>=2.0 && Zc_btm_H<=3.5) 
		{
			for (int i=0; i<ncolsl; i++) {
				x=*(Datapointl+2*i+1);
				y=*(Datapointl+2*i);
				if (tcon+ncolsl+ncolsr-i>tcon_pre) {
					if (y>=*D_pre_ptr) {
						*(line_ptr_l+i)=0;
						continue;
					}
					x1=A/(y-*D_pre_ptr)+Bleft*y+Cleft;
					if (abs(x-x1)>10) {
						*(line_ptr_l+i)=0;
						continue;
					}
					Cp1=(*D_pre_ptr-y)*0.5;
					Cp2=(Eleft-x1)*0.5+Bleft*y;
					Cp3=Fleft+(*D_pre_ptr*x1+Eleft*y)*0.5;
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
					x1=A/(y-*D_pre_ptr)+Bright*y+Cright;
					if (abs(x-x1)>10) {
						*(line_ptr_r+i)=0;
						continue;
					}
					Cp1=(*D_pre_ptr-y)*0.5;
					Cp2=(Eright-x1)*0.5+Bright*y;
					Cp3=Fright+(*D_pre_ptr*x1+Eright*y)*0.5;
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
		if (bcd_l[1]<C_limit && (eleft-eright)/c<0 && ((eleft-eright)/c*top_row+(fpleft-fpright)/c)>0
			&& zc_btm_l>=2.0 && zc_btm_l<=3.5) 
		{
			for (int i=0; i<ncolsl; i++) {
				if (tcon1+ncolsl+ncolsr-i>tcon_pre) {
					x=*(Datapointl+2*i+1);
					y=*(Datapointl+2*i);
					t_dist=1.0-(y-top_row)/440*9;
					dist=pow(c*x+eleft*y+fpleft,2)/denoleft;
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
					dist=pow(c*x+eright*y+fpright,2)/denoright;
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
			if ((tcon*prob_H)>(tcon_pre*prob) || (tcon1*prob_l)>(tcon_pre*prob))
			{
				if ((1.0*tcon*prob_H)>(1.0*tcon1*prob_l)) {
					tcon_pre=tcon;
					prob=prob_H;
					//*total_num_ptr=tcon_pre;
					*(hypo_ptr)=A;	*(hypo_ptr+1)=Bleft;	*(hypo_ptr+2)=Cleft;	*(hypo_ptr+3)=Bright;	*(hypo_ptr+4)=Cright;
					*(hypo_ptr+5)=0;	*(hypo_ptr+6)=0;	*(hypo_ptr+7)=0;		*(hypo_ptr+8)=0;		*(hypo_ptr+9)=0;
					for (int i=0;i<ncolsl;i++)
						*(FinalLinel+i)=*(line_ptr_l+i);
					for (int i=0;i<ncolsr;i++)
						*(FinalLiner+i)=*(line_ptr_r+i);
				}
				else {
					tcon_pre=tcon1;
					prob=prob_l;
					*(hypo_ptr)=0;		*(hypo_ptr+1)=0;		*(hypo_ptr+2)=0;		*(hypo_ptr+3)=0;		*(hypo_ptr+4)=0;
					*(hypo_ptr+5)=c;	*(hypo_ptr+6)=eleft;	*(hypo_ptr+7)=fleft;	*(hypo_ptr+8)=eright;	*(hypo_ptr+9)=fright;
					for (int i=0;i<ncolsl;i++)
						*(FinalLinel+i)=*(line_ptr_l1+i);
					for (int i=0;i<ncolsr;i++)
						*(FinalLiner+i)=*(line_ptr_r1+i);
				}	
			}
		}
		delete [] line_ptr_l1;
		delete [] line_ptr_l;
		delete [] line_ptr_r1;
		delete [] line_ptr_r;
	}
}