// mex file to run RANSAC to improve speed
// input: 2xn matrix with 1st row being x and 2nd row being y. n refers to the number of points
// ouput: 1xn matrix (or row vector) to indicate whether the corresponding points are inliers (1) or not (0)
// fitting left and right img seperately
#include "mex.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <vector>

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

void Invert2(double *v,double *u, double *bedf)
{
 double mat[16];
 double tmp[12]; /* temp array for pairs */
 double src[16]; /* array of transpose source matrix */
 double dst[16];
 double B[4];
 double det; /* determinant */
 for (int i=0; i<4; i++) {
	 mat[i]=pow(*(v+i),2);
	 mat[i+4]=*(v+i);
	 mat[i+8]=*(u+i);
	 mat[i+12]=1;
	 B[i]=(*(u+i))*(*(v+i));
 }

 /* transpose matrix */
 for (int i = 0; i < 4; i++) {
 src[i] = mat[i*4];
 src[i + 4] = mat[i*4 + 1];
 src[i + 8] = mat[i*4 + 2];
 src[i + 12] = mat[i*4 + 3];
 }
 /* calculate pairs for first 8 elements (cofactors) */
 tmp[0] = src[10] * src[15];
 tmp[1] = src[11] * src[14];
 tmp[2] = src[9] * src[15];
 tmp[3] = src[11] * src[13];
 tmp[4] = src[9] * src[14];
 tmp[5] = src[10] * src[13];
 tmp[6] = src[8] * src[15];
 tmp[7] = src[11] * src[12];
 tmp[8] = src[8] * src[14];
 tmp[9] = src[10] * src[12];
 tmp[10] = src[8] * src[13];
 tmp[11] = src[9] * src[12];
 /* calculate first 8 elements (cofactors) */
 dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
 dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
 dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
 dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
 dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
 dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
 dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
 dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
 dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
 dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
 dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
 dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
 dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
 dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
 dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
 dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
 /* calculate pairs for second 8 elements (cofactors) */
 tmp[0] = src[2]*src[7];
 tmp[1] = src[3]*src[6];
 tmp[2] = src[1]*src[7];
 tmp[3] = src[3]*src[5];
 tmp[4] = src[1]*src[6];
 tmp[5] = src[2]*src[5];
 tmp[6] = src[0]*src[7];
 tmp[7] = src[3]*src[4];
 tmp[8] = src[0]*src[6];
 tmp[9] = src[2]*src[4];
 tmp[10] = src[0]*src[5];
 tmp[11] = src[1]*src[4];
 /* calculate second 8 elements (cofactors) */
 dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
 dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
 dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
 dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
 dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
 dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
 dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
 dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
 dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
 dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
 dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
 dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
 dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
 dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
 dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
 dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
 /* calculate determinant */
 det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];
 /* calculate matrix inverse */
 det = ((det==0)?0:1/det);
 for (int j = 0; j < 16; j++)
 dst[j] *= det;
 for (int i=0; i<4; i++) {
        bedf[i]=(*(dst+i))*(*B)+(*(dst+i+4))*(*(B+1))+(*(dst+i+8))*(*(B+2))+(*(dst+i+12))*(*(B+3));
    }


}
/*
void LRegression(double *v,double *u, double *cefp)
{
	double vbar=0, ubar=0, vubar=0, u2bar=0;
	for (int i=0;i<4;i++) {
		ubar=ubar+*(u+i);
		vbar=vbar+*(v+i);
		vubar=vubar+*(u+i)**(v+i);
		u2bar=u2bar+*(u+i)**(u+i);
	}
	ubar=ubar/4; vbar=vbar/4; vubar=vubar/2; u2bar=u2bar/4;
	*(cefp+1)=1.0;
	*cefp=-(vubar-vbar*ubar)/(u2bar-ubar*ubar);
	*(cefp+2)=-(vbar+*cefp*ubar);
}
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *Datapoint, *FinalLine, *SubDatapoint, *hypo_ptr;
    size_t ncols, sub_ncols;
    //size_t Randpoint[4];
    int j;
	int tcon_pre=0;
    //mxArray *hypo;
	//double Cp1, Cp2, Cp3, x, y, x1, deno,t_dist;
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

	//plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
	//total_num_ptr=mxGetPr(plhs[2]);

	//plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
	//dumy=mxGetPr(plhs[2]);
   // Randpoint[0]=0; Randpoint[1]=0; Randpoint[2]=0; Randpoint[3]=0;
   // distsqr=mxCreateDoubleMatrix(1,(int)ncols,mxREAL);
    //distsqr_ptr=mxGetPr(distsqr);
    //line=mxCreateDoubleMatrix(1,ncols,mxREAL);
    //line_ptr=mxGetPr(line);
	//distsqr1=mxCreateDoubleMatrix(1,(int)ncols,mxREAL);
    //distsqr_ptr1=mxGetPr(distsqr1);
    //line1=mxCreateDoubleMatrix(1,ncols,mxREAL);
    //line_ptr1=mxGetPr(line1);
	//V=mxCreateDoubleMatrix(4,1,mxREAL);
	//U=mxCreateDoubleMatrix(4,1,mxREAL);
	//v=mxGetPr(V);
	//u=mxGetPr(U);
	//BEDF=mxCreateDoubleMatrix(4,1,mxREAL);
	//bedf=mxGetPr(BEDF);

	#pragma omp parallel for num_threads(4)
    for(j=0;j<960;j++)
    {
		int tcon=0; int tcon1=0;
		//for (int i=0;i<4;i++) {
		//	Randpoint[i]=rand() % sub_ncols;
		//	*(v+i)=*(SubDatapoint+2*Randpoint[i]);
		//	*(u+i)=*(SubDatapoint+2*Randpoint[i]+1);
		//}
		//*v=-227; *(v+1)=-156; *(v+2)=-91; *(v+3)=-22;
		//*u=15; *(u+1)=116; *(u+2)=207; *(u+3)=305;
		int rand1=rand() % sub_ncols;
		//mxArray *V=mxCreateDoubleMatrix(4,1,mxREAL);
		//mxArray *U=mxCreateDoubleMatrix(4,1,mxREAL);
		//double *v=mxGetPr(V);
		//double *u=mxGetPr(U);
		double v[4],u[4];
		*v=*(SubDatapoint+2*rand1);
		*u=*(SubDatapoint+2*rand1+1);
		rand1=rand() % sub_ncols;
		*(v+1)=*(SubDatapoint+2*rand1);
		*(u+1)=*(SubDatapoint+2*rand1+1);
		rand1=rand() % sub_ncols;
		*(v+2)=*(SubDatapoint+2*rand1);
		*(u+2)=*(SubDatapoint+2*rand1+1);
		rand1=rand() % sub_ncols;
		*(v+3)=*(SubDatapoint+2*rand1);
		*(u+3)=*(SubDatapoint+2*rand1+1);
		
		//mxArray *BEDF=mxCreateDoubleMatrix(4,1,mxREAL);
		//double *bedf=mxGetPr(BEDF);
		double bedf[4];
		//CallInv(v,u,bedf);
		Invert2(v,u,bedf);

		double B=*bedf;double E=*(bedf+1); double D=*(bedf+2);double F=*(bedf+3);
		double C=E+B*D; double A=F+C*D;
		double e,c,f,fp;
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
		//double cefp[3];
		//LRegression(v,u,cefp);
		//c=cefp[0]; e=cefp[1]; fp=cefp[2];
        double deno=c*c+e*e;

		double Cp1, Cp2, Cp3, x, y, x1,t_dist,dist;
		//mxArray *line1=mxCreateDoubleMatrix(1,ncols,mxREAL);
		//double *line_ptr1=mxGetPr(line1);
		double *line_ptr1=new double[ncols];
		//mxArray *line=mxCreateDoubleMatrix(1,ncols,mxREAL);
		//double *line_ptr=mxGetPr(line);
		//double line_ptr[ncols];
		double *line_ptr=new double[ncols];
		//D=-300;
		if (D>-120) {
			for (int i=0; i<ncols; i++) {
				x=*(Datapoint+2*i+1);
				y=*(Datapoint+2*i);
				if (tcon1+ncols-i>tcon_pre) {
					t_dist=1.0-(y+20.0)/55;
					dist=pow(c*x+e*y+fp,2)/deno;
					if (dist<t_dist){
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
					dist=pow(x*Cp1+y*Cp2+Cp3,2)/4.0/(Cp1*Cp1+Cp2*Cp2);
					t_dist=0.5-(y+20.0)/110;
					if (dist<t_dist){
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
			for (int i=0; i<ncols; i++) {
				if (tcon1+ncols-i>tcon_pre) {
					x=*(Datapoint+2*i+1);
					y=*(Datapoint+2*i);
					t_dist=1.0-(y+20.0)/55;
					dist=pow(c*x+e*y+fp,2)/deno;
					if (dist<t_dist) {
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
		#pragma omp critical
		{
		if (tcon>tcon_pre || tcon1>tcon_pre)
		{
			if (tcon>tcon1) {
			tcon_pre=tcon;
			//*total_num_ptr=tcon_pre;
			*(hypo_ptr)=B;*(hypo_ptr+1)=E;*(hypo_ptr+2)=D;*(hypo_ptr+3)=F;
			*(hypo_ptr+4)=0;*(hypo_ptr+5)=0;*(hypo_ptr+6)=0;
			for (int i=0;i<ncols;i++)
				*(FinalLine+i)=*(line_ptr+i);
			//if (tcon>0.9*ncols)
			//	j=250;
			}
			else {
				tcon_pre=tcon1;
				//*total_num_ptr=tcon_pre;
			*(hypo_ptr)=0;*(hypo_ptr+1)=0;*(hypo_ptr+2)=0;*(hypo_ptr+3)=0;
			*(hypo_ptr+4)=c;*(hypo_ptr+5)=e;*(hypo_ptr+6)=f;
			for (int i=0;i<ncols;i++)
				*(FinalLine+i)=*(line_ptr1+i);
			//if (tcon1>0.9*ncols)
			//	j=250;
			}
			
		}
		}
		delete [] line_ptr1;
		delete [] line_ptr;
	}
	//mxDestroyArray(line);
    //mxDestroyArray(distsqr);
	//mxDestroyArray(V);
	//mxDestroyArray(U);
	//mxDestroyArray(BEDF);
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