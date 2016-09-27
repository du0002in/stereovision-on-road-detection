// mex file to run RANSAC to improve speed
// input: 2xn matrix with 1st row being x and 2nd row being y. n refers to the number of points
// ouput: 1xn matrix (or row vector) to indicate whether the corresponding points are inliers (1) or not (0)
#include "mex.h"
#include <math.h>
#include <time.h>
#include <string.h>

void CallInv(double *v, double *u, double *bedf)
{
    mxArray *A, *invA, *B;
    mxArray *ppLhs[1];
    double *ptr, *A_ptr, *B_ptr; int i;
    A=mxCreateDoubleMatrix(4,4,mxREAL);
	A_ptr=mxGetPr(A);
	B=mxCreateDoubleMatrix(4,1,mxREAL);
	B_ptr=mxGetPr(B);
	for (i=0;i<4;i++) {
		*(A_ptr+i)=pow(*(v+i),2);
		*(A_ptr+i+4)=*(v+i);
		*(A_ptr+i+8)=*(u+i);
		*(A_ptr+i+12)=1;
		*(B_ptr+i)=(*(u+i))*(*(v+i));
	}
   // memcpy(mxGetPr(A), Data, sizeof(double)*4*4);
    mexCallMATLAB(1, ppLhs, 1, &A, "inv");
    invA=ppLhs[0];
	ptr=mxGetPr(invA);
    for (i=0; i<4; i++) {
        bedf[i]=(*(ptr+i))*(*B_ptr)+(*(ptr+i+4))*(*(B_ptr+1))+(*(ptr+i+8))*(*(B_ptr+2))+(*(ptr+i+12))*(*(B_ptr+3));
    }
    mxDestroyArray(A);
	mxDestroyArray(B);
    mxDestroyArray(invA);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *Datapoint, *v, *u, *FinalLine, *distsqr_ptr, *line_ptr, *bedf, *SubDatapoint, *hypo_ptr, *line_ptr1, *distsqr_ptr1, *total_num_ptr;
    size_t ncols, sub_ncols;
    size_t Randpoint[4];
    int i,j, tcon=0, tcon_pre=0, *s, tcon1=0, tcon_pre1=0;
    mxArray *line, *distsqr, *V, *U, *BEDF, *hypo, *line1, *distsqr1;
	double A,B,C,D,E,F,Cp1, Cp2, Cp3, x, y, x1, c,e,f,fp,deno,t_dist;
 /*
    v=mxGetPr(prhs[0]);
	u=mxGetPr(prhs[1]);
    plhs[0]=mxCreateDoubleMatrix(4,1,mxREAL);
    FinalLine=mxGetPr(plhs[0]);
    CallInv(v,u, FinalLine);
	if (~mxIsDouble(bedf))
		mexPrintf("no solution");
*/

	srand(time(NULL));
    Datapoint=mxGetPr(prhs[0]);
	SubDatapoint=mxGetPr(prhs[1]);
    ncols=mxGetN(prhs[0]);
	sub_ncols=mxGetN(prhs[1]);
    plhs[0]=mxCreateDoubleMatrix(1,ncols,mxREAL);
    FinalLine=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(7,1,mxREAL);
    hypo_ptr=mxGetPr(plhs[1]);
	plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
	total_num_ptr=mxGetPr(plhs[2]);

    Randpoint[0]=0; Randpoint[1]=0; Randpoint[2]=0; Randpoint[3]=0;
    distsqr=mxCreateDoubleMatrix(1,(int)ncols,mxREAL);
    distsqr_ptr=mxGetPr(distsqr);
    line=mxCreateDoubleMatrix(1,ncols,mxREAL);
    line_ptr=mxGetPr(line);
	distsqr1=mxCreateDoubleMatrix(1,(int)ncols,mxREAL);
    distsqr_ptr1=mxGetPr(distsqr1);
    line1=mxCreateDoubleMatrix(1,ncols,mxREAL);
    line_ptr1=mxGetPr(line1);
	V=mxCreateDoubleMatrix(4,1,mxREAL);
	U=mxCreateDoubleMatrix(4,1,mxREAL);
	v=mxGetPr(V);
	u=mxGetPr(U);
	BEDF=mxCreateDoubleMatrix(4,1,mxREAL);
	bedf=mxGetPr(BEDF);
    for(j=0;j<250;j++)
    {
		tcon=0; tcon1=0;
		for (i=0;i<4;i++) {
			Randpoint[i]=rand() % sub_ncols;
			*(v+i)=*(SubDatapoint+2*Randpoint[i]);
			*(u+i)=*(SubDatapoint+2*Randpoint[i]+1);
		}
		//*v=-213; *(v+1)=-148; *(v+2)=-90; *(v+3)=-179;
		//*u=41; *(u+1)=131; *(u+2)=212; *(u+3)=89;
		CallInv(v,u,bedf);
		B=*bedf; E=*(bedf+1); D=*(bedf+2);F=*(bedf+3);
		C=E+B*D; A=F+C*D;

		if (*u==*(u+3))
        {
            e=0; c=1; f=-2**u; fp=f/2;
        }
        else
        {
            c=-(*v-*(v+3))/(*u-*(u+3));
            e=1;
            f=2*(-c**u-*v);
            fp=f/2;
        }
        deno=c*c+e*e;

		if (D>-120) {
			for (i=0; i<ncols; i++) {
				x=*(Datapoint+2*i+1);
				y=*(Datapoint+2*i);
				if (tcon1+ncols-i>tcon_pre) {
					t_dist=1.0-(y+20.0)/55;
					*(distsqr_ptr1+i)=pow(c*x+e*y+fp,2)/deno;
					if (*(distsqr_ptr1+i)<t_dist){
						tcon1++;
						*(line_ptr1+i)=1;}
					else
						*(line_ptr1+i)=0;
				}
				if (tcon+ncols-i>tcon_pre) {
					if (y>=D) {
						*(line_ptr+i)=0;
						continue;
					}
					x1=A/(y-D)+B*y+C;
					if (abs(x-x1)>10) {
						*(line_ptr+i)=0;
						continue;
					}
					Cp1=(D-y)*0.5;
					Cp2=(E-x1)*0.5+B*y;
					Cp3=F+(D*x1+E*y)*0.5;
					*(distsqr_ptr+i)=pow(x*Cp1+y*Cp2+Cp3,2)/4.0/(Cp1*Cp1+Cp2*Cp2);
					t_dist=0.5-(y+20.0)/110;
					if (*(distsqr_ptr+i)<t_dist){
						tcon++;
						*(line_ptr+i)=1;}
					else
						*(line_ptr+i)=0;
				}
				if (tcon1+ncols-i<=tcon_pre && tcon+ncols-i<=tcon_pre)
					break;
			}

		}
		else {
			for (i=0; i<ncols; i++) {
				if (tcon1+ncols-i>tcon_pre) {
					x=*(Datapoint+2*i+1);
					y=*(Datapoint+2*i);
					t_dist=1.0-(y+20.0)/55;
					*(distsqr_ptr1+i)=pow(c*x+e*y+fp,2)/deno;
					if (*(distsqr_ptr1+i)<t_dist) {
						tcon1++;
						*(line_ptr1+i)=1;}
					else
						*(line_ptr1+i)=0;
					*(line_ptr+i)=0;
				}
				else
					break;
			}
		}

		if (tcon>tcon_pre || tcon1>tcon_pre)
		{
			if (tcon>tcon1) {
			tcon_pre=tcon;
			*total_num_ptr=tcon_pre;
			*(hypo_ptr)=B;*(hypo_ptr+1)=E;*(hypo_ptr+2)=D;*(hypo_ptr+3)=F;
			*(hypo_ptr+4)=0;*(hypo_ptr+5)=0;*(hypo_ptr+6)=0;
			for (i=0;i<ncols;i++)
				*(FinalLine+i)=*(line_ptr+i);
			//if (tcon>0.9*ncols)
			//	break;
			}
			else {
				tcon_pre=tcon1;
				*total_num_ptr=tcon_pre;
			*(hypo_ptr)=0;*(hypo_ptr+1)=0;*(hypo_ptr+2)=0;*(hypo_ptr+3)=0;
			*(hypo_ptr+4)=c;*(hypo_ptr+5)=e;*(hypo_ptr+6)=f;
			for (i=0;i<ncols;i++)
				*(FinalLine+i)=*(line_ptr1+i);
			
			//if (tcon1>0.9*ncols)
			//	break;
			}
		}
	}
	mxDestroyArray(line);
	mxDestroyArray(line1);
    mxDestroyArray(distsqr);
	mxDestroyArray(V);
	mxDestroyArray(U);
	mxDestroyArray(BEDF);
}

/*
        if (u[0]==u[1])
        {
            e=0; c=1; f=-2*u[0]; fp=f/2;
        }
        else
        {
            c=-(v[0]-v[1])/(u[0]-u[1]);
            e=1;
            f=2*(-c*u[0]-v[0]);
            fp=f/2;
        }
        deno=c*c+e*e;

        for (i=0;i<ncols;i++)
        {
            *(distsqr_ptr+i)=pow(c**(Datapoint+2*i)+e*(*(Datapoint+2*i+1))+fp,2)/deno;
            if (*(distsqr_ptr+i)<10){
				tcon++;
				*(line_ptr+i)=1;}
            else
                *(line_ptr+i)=0;
        }
		if (tcon<50)
			continue;
		if (tcon>tcon_pre)
		{
			tcon_pre=tcon;
			for (i=0;i<ncols;i++)
				*(FinalLine+i)=*(line_ptr+i);
			if (tcon>0.9*ncols)
				break;
		}    
    }
	mxDestroyArray(line);
    mxDestroyArray(distsqr);

}
*/