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

void GetModel(double *v,double *u, double *BCD_L)
{
 double L[9], dU[3], inv_L[9];
/* H[0]=v[0]-v[2];	H[3]=1/(v[0]-*D_pre)-1/(v[2]-*D_pre);	H[6]=v[2]-*D_pre;	dU[0]=u[0]-u[2];
 H[1]=v[1]-v[3];	H[4]=1/(v[1]-*D_pre)-1/(v[3]-*D_pre);	H[7]=v[3]-*D_pre;	dU[1]=u[1]-u[3];
 H[2]=v[0]-v[3];	H[5]=1/(v[0]-*D_pre)-1/(v[3]-*D_pre);	H[8]=v[3]-*D_pre;	dU[2]=u[0]-u[3];*/
 L[0]=v[0]-v[2];	L[3]=v[2];	L[6]=-1;	dU[0]=u[0]-u[2];
 L[1]=v[1]-v[3];	L[4]=v[3];	L[7]=-1;	dU[1]=u[1]-u[3];
 L[2]=v[0]-v[3];	L[5]=v[3];	L[8]=-1;	dU[2]=u[0]-u[3];
 //inverse3x3(H,inv_H);
 inverse3x3(L,inv_L);
 /*BAC_H[0]=inv_H[0]*dU[0]+inv_H[3]*dU[1]+inv_H[6]*dU[2];
 BAC_H[1]=inv_H[1]*dU[0]+inv_H[4]*dU[1]+inv_H[7]*dU[2];
 BAC_H[2]=inv_H[2]*dU[0]+inv_H[5]*dU[1]+inv_H[8]*dU[2];*/
 
 BCD_L[0]=inv_L[0]*dU[0]+inv_L[3]*dU[1]+inv_L[6]*dU[2];
 BCD_L[1]=inv_L[1]*dU[0]+inv_L[4]*dU[1]+inv_L[7]*dU[2];
 BCD_L[2]=inv_L[2]*dU[0]+inv_L[5]*dU[1]+inv_L[8]*dU[2];
 BCD_L[2]=BCD_L[2]/BCD_L[1];
}

void Sort_4_Array(double *x, double *y) //sort based on y in accending order
{
	int pos_min, n=4; double ytmp, xtmp;
	for (int i=0; i<n-1; i++)
	{
		pos_min=i;
		for (int j=i+1; j<n; j++)
		{
			if (y[j]<y[pos_min])
				pos_min=j;
		}
		if (pos_min != i)
		{
			ytmp=y[i]; xtmp=x[i];
			y[i]=y[pos_min]; y[pos_min]=ytmp;
			x[i]=x[pos_min]; x[pos_min]=xtmp;
		}
	}
}

void Random_4_Unique(int *arr, int interval)
{
	arr[0]=rand() % interval;
	arr[1]=rand() % interval;
	while (arr[1]==arr[0])
		arr[1]=rand() % interval;
	arr[2]=rand() % interval;
	while (arr[2]==arr[0] || arr[2]==arr[1])
		arr[2]=rand() % interval;
	arr[3]=rand() % interval;
	while (arr[3]==arr[0] || arr[3]==arr[1] || arr[3]==arr[2])
		arr[3]=rand() % interval;
}

void valid_line_intersect(double *x,double *y,int *valid, double *y_int, double focal, double phi, double *ab, double *ab_p)
{
	double theta1, theta2;
	double xx[2], yy[2];
	xx[0]=x[0]; xx[1]=x[1];
	yy[0]=y[0]; yy[1]=y[1];
	ab[0]=(xx[0]-xx[1])/(yy[0]-yy[1]);
	ab[1]=xx[0]-ab[0]*yy[0];
	//GetLine(xx,yy,ab);
	theta1=atan((ab[0]*sin(phi)-ab[1]/focal*cos(phi))/(sin(phi)*sin(phi)-cos(phi)*cos(phi)));

	xx[0]=x[2]; xx[1]=x[3];
	yy[0]=y[2]; yy[1]=y[3];
	ab_p[0]=(xx[0]-xx[1])/(yy[0]-yy[1]);
	ab_p[1]=xx[0]-ab_p[0]*yy[0];
	//GetLine(xx,yy,ab);
	theta2=atan((ab_p[0]*sin(phi)-ab_p[1]/focal*cos(phi))/(sin(phi)*sin(phi)-cos(phi)*cos(phi)));

	if (abs(theta1-theta2)<0.1571) {
		valid[0]=0;
		return;
	}


	double a=(x[0]-x[1])*(y[2]-y[3])-(x[2]-x[3])*(y[0]-y[1]);
	double b=(x[0]-x[1])*(y[2]-y[3])*y[0]-(x[2]-x[3])*(y[0]-y[1])*y[3]
				+(x[3]-x[0])*(y[0]-y[1])*(y[2]-y[3]);
	*y_int=b/a;
	if ((*y_int-y[1])>5 && (y[2]-*y_int)>5) //the intersect is between 2nd and 3rd points
		*valid=1;
	else {
		*valid=0;
		return;
	}

//	if (y[0]!=y[1]) *x_int=(x[0]-x[1])*(*y_int-y[0])/(y[0]-y[1])+x[0];
//	else *x_int=(x[2]-x[3])*(*y_int-y[3])/(y[2]-y[3])+x[3];
}

/*void valid_line_intersect(double *x,double *y,int *valid, double *y_int)
{
	double theta1, theta2;
	if (x[0]==x[1]) theta1=3.1415926;
	else theta1=atan((y[0]-y[1])/(x[0]-x[1]));
	if (x[2]==x[3]) theta2=3.1415926;
	else theta2=atan((y[2]-y[3])/(x[2]-x[3]));
	if (abs(theta1-theta2)<0.1571) {
		*valid=0;
		return;
	}

	double a=(x[0]-x[1])*(y[2]-y[3])-(x[2]-x[3])*(y[0]-y[1]);
	double b=(x[0]-x[1])*(y[2]-y[3])*y[0]-(x[2]-x[3])*(y[0]-y[1])*y[3]
				+(x[3]-x[0])*(y[0]-y[1])*(y[2]-y[3]);
	*y_int=b/a;
	if ((*y_int-y[1])>5 && (y[2]-*y_int)>5) //the intersect is between 2nd and 3rd points
		*valid=1;
	else {
		*valid=0;
		return;
	}

//	if (y[0]!=y[1]) *x_int=(x[0]-x[1])*(*y_int-y[0])/(y[0]-y[1])+x[0];
//	else *x_int=(x[2]-x[3])*(*y_int-y[3])/(y[2]-y[3])+x[3];
}*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *Datapointl, *Datapointr, *FinalLinel, *FinalLiner, *SubDatapointl, *SubDatapointr, *hypo_ptr, *D_pre_ptr, *H_pre_ptr, *num_ite_ptr;
    double *focalPtr, *dist_left_rightPtr, *cc_yPtr, *y_int_ptr;
	size_t ncolsl, ncolsr, sub_ncolsl, sub_ncolsr;
    //size_t Randpoint[4];
    int j, num_ite;
	int tcon_pre=0;
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
	num_ite=int(*num_ite_ptr)*4;

	focalPtr=mxGetPr(prhs[7]);
	dist_left_rightPtr=mxGetPr(prhs[8]);
	cc_yPtr=mxGetPr(prhs[9]);
	double phi=atan(*D_pre_ptr/(*focalPtr));
	top_row=*cc_yPtr-40;

    plhs[0]=mxCreateDoubleMatrix(1,ncolsl,mxREAL);
    FinalLinel=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(1,ncolsr,mxREAL);
    FinalLiner=mxGetPr(plhs[1]);
	plhs[2]=mxCreateDoubleMatrix(5,2,mxREAL);
    hypo_ptr=mxGetPr(plhs[2]);
	plhs[3]=mxCreateDoubleMatrix(1,2,mxREAL);
    y_int_ptr=mxGetPr(plhs[3]);
	
	for (int i=0;i<ncolsl;i++)
		*(FinalLinel+i)=0;
	for (int i=0;i<ncolsr;i++)
		*(FinalLiner+i)=0;
	for (int i=0;i<10;i++)
		*(hypo_ptr+i)=0;
	*y_int_ptr=0; *(y_int_ptr+1)=0;

	// C_limit is related to distance btw left and right camera 0.2566. Focal=852. 1.2 is a buffering factor.
	double C_limit=1/sqrt(1+(*D_pre_ptr/(*focalPtr))*(*D_pre_ptr/(*focalPtr)))*(*dist_left_rightPtr)/(*H_pre_ptr)*1.5;
	double ymax=(top_row<*D_pre_ptr?top_row:(*D_pre_ptr-10));
	#pragma omp parallel for num_threads(4)
    for(j=0;j<num_ite;j++)
    {
		int tcon1=0;
		//int rand1=rand() % sub_ncolsl;
		// choose 4 pt from left img
		double vl[4],ul[4];
		int rand_num[4];
		Random_4_Unique(rand_num, sub_ncolsl);
		*vl=*(SubDatapointl+2*rand_num[0]);
		*ul=*(SubDatapointl+2*rand_num[0]+1);
		*(vl+1)=*(SubDatapointl+2*rand_num[1]);
		*(ul+1)=*(SubDatapointl+2*rand_num[1]+1);
		*(vl+2)=*(SubDatapointl+2*rand_num[2]);
		*(ul+2)=*(SubDatapointl+2*rand_num[2]+1);
		*(vl+3)=*(SubDatapointl+2*rand_num[3]);
		*(ul+3)=*(SubDatapointl+2*rand_num[3]+1);
		// choose 4 pt from right img
		double vr[4],ur[4];
		Random_4_Unique(rand_num, sub_ncolsr);
		*vr=*(SubDatapointr+2*rand_num[0]);
		*ur=*(SubDatapointr+2*rand_num[0]+1);
		*(vr+1)=*(SubDatapointr+2*rand_num[1]);
		*(ur+1)=*(SubDatapointr+2*rand_num[1]+1);
		*(vr+2)=*(SubDatapointr+2*rand_num[2]);
		*(ur+2)=*(SubDatapointr+2*rand_num[2]+1);
		*(vr+3)=*(SubDatapointr+2*rand_num[3]);
		*(ur+3)=*(SubDatapointr+2*rand_num[3]+1);
		// sort the two sets based on u value so that the upper two becomes one set and botom two is the other set
		Sort_4_Array(ul,vl);
		Sort_4_Array(ur,vr);
		int valid_l, valid_r; double y_int_l=0, y_int_r=0;
		double ab_dumy[2], ab_p_dumy[2];
		valid_line_intersect(ul,vl,&valid_l,&y_int_l,*focalPtr,phi, ab_dumy, ab_p_dumy);
		valid_line_intersect(ur,vr,&valid_r,&y_int_r,*focalPtr,phi, ab_dumy, ab_p_dumy);


		if ((valid_l==1) && (valid_r)==1) {
		double u[4],v[4],up[4],vp[4];
		u[0]=ul[0]; u[1]=ul[1]; u[2]=ur[0]; u[3]=ur[1];      //btm lines
		v[0]=vl[0]; v[1]=vl[1]; v[2]=vr[0]; v[3]=vr[1];
		up[0]=ul[2]; up[1]=ul[3]; up[2]=ur[2]; up[3]=ur[3];
		vp[0]=vl[2]; vp[1]=vl[3]; vp[2]=vr[2]; vp[3]=vr[3];  //top lines
		
		double bcd_l[3], bcd_lp[3]; //l -- line, lp -- line prime
		GetModel(v,u,bcd_l);
		GetModel(vp,up,bcd_lp);
		//left and right Line
		double c=1,		eleft=-bcd_l[0],	fpleft=-eleft*v[0]-u[0],		fleft=fpleft*2;
		double eright=-(bcd_l[0]-bcd_l[1]),			fpright=fpleft-bcd_l[1]*bcd_l[2],		fright=fpright*2;
		double cp=1,	eleftp=-bcd_lp[0],	fpleftp=-eleftp*vp[0]-up[0],	fleftp=fpleftp*2;
		double erightp=-(bcd_lp[0]-bcd_lp[1]),		fprightp=fpleftp-bcd_lp[1]*bcd_lp[2],	frightp=fprightp*2;

        double denoleft=c*c+eleft*eleft;
		double denoright=c*c+eright*eright;
		double denoleftp=cp*cp+eleftp*eleftp;
		double denorightp=cp*cp+erightp*erightp;

		double Cp1, Cp2, Cp3, x, y, x1,t_dist,dist,distp;
		double *line_ptr_l1=new double[ncolsl];
		//double *line_ptr_l=new double[ncolsl];
		double *line_ptr_r1=new double[ncolsr];
		//double *line_ptr_r=new double[ncolsr];
		
		if ((bcd_l[1]<C_limit && (eleft-eright)/c<0 && ((eleft-eright)/c*top_row+(fpleft-fpright)/c)>0) && 
			(bcd_lp[1]<C_limit && (eleftp-erightp)/cp<0 && ((eleftp-erightp)/cp*top_row+(fpleftp-fprightp)/cp)>0)){
			for (int i=0; i<ncolsl; i++) {
				if (tcon1+ncolsl+ncolsr-i>tcon_pre) {
					x=*(Datapointl+2*i+1);
					y=*(Datapointl+2*i);
					t_dist=1.0-(y-top_row)/440*9;
					if (y<y_int_l)
						dist=pow(c*x+eleft*y+fpleft,2)/denoleft;
					else
						dist=pow(cp*x+eleftp*y+fpleftp,2)/denoleftp;
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
					if (y<y_int_r)
						dist=pow(c*x+eright*y+fpright,2)/denoright;
					else
						dist=pow(cp*x+erightp*y+fprightp,2)/denorightp;
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
			*(hypo_ptr)=c;		*(hypo_ptr+1)=eleft;		*(hypo_ptr+2)=fleft;	*(hypo_ptr+3)=eright;		*(hypo_ptr+4)=fright;   //btm lines
			*(hypo_ptr+5)=cp;	*(hypo_ptr+6)=eleftp;		*(hypo_ptr+7)=fleftp;	*(hypo_ptr+8)=erightp;		*(hypo_ptr+9)=frightp;  //top lines
			*y_int_ptr=y_int_l; *(y_int_ptr+1)=y_int_r;

			for (int i=0;i<ncolsl;i++)
				*(FinalLinel+i)=*(line_ptr_l1+i);
			for (int i=0;i<ncolsr;i++)
				*(FinalLiner+i)=*(line_ptr_r1+i);

		}
		}
		delete [] line_ptr_l1;
		//delete [] line_ptr_l;
		delete [] line_ptr_r1;
		//delete [] line_ptr_r;
		}
	}
}