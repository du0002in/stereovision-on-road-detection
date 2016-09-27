#include <string.h>
#include "mex.h"
#include <omp.h>
#include <math.h>
#include <time.h>

void DisplayMatrix(char *Name, double *Data, int M, int N)
{
int m, n;
mexPrintf("%s = nn", Name);
for(m = 0; m < M; m++, mexPrintf("nn"))
for(n = 0; n < N; n++)
mexPrintf("%8.4f ", Data[m + M*n]);
}

void CallQR(double *Data, int M, int N)
{ /* Perform QR factorization by calling the MATLAB function */
mxArray *Q, *R, *A;
mxArray *ppLhs[1];
double *ptr;
//DisplayMatrix("Input", Data, M, N);
A = mxCreateDoubleMatrix(4, 4, mxREAL); /* Put input in an mxArray */
memcpy(mxGetPr(A), Data, sizeof(double)*4*4);
mexCallMATLAB(1, ppLhs, 1, &A, "inv"); /* Call MATLAB's qr function */
Q = ppLhs[0];
ptr=mxGetPr(Q);
//R = ppLhs[1];
DisplayMatrix("Q", mxGetPr(Q), M, N);
//DisplayMatrix("R", mxGetPr(R), M, N);
mxDestroyArray(Q);
mxDestroyArray(A);
}
void Invert2(double *mat, double *dst)
{
 double tmp[12]; /* temp array for pairs */
 double src[16]; /* array of transpose source matrix */
 double det; /* determinant */
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
 det = 1/det;
 for (int j = 0; j < 16; j++)
 dst[j] *= det;
}

double * inverse3x3(double *a)
{
	double det;static double in[9];
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
	return in;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//#define M_IN prhs[0]
//CallQR(mxGetPr(M_IN), mxGetM(M_IN), mxGetN(M_IN));
	double *APtr, *BPtr; int i; double c=0; int sub_ncols=10;
	//mxArray *A;
	//A=mxDuplicateArray(prhs[0]);
	APtr=mxGetPr(prhs[0]);
	 plhs[0]=mxCreateDoubleMatrix(3,3,mxREAL);
	BPtr=mxGetPr(plhs[0]);
	//PIII_Inverse_4x4((float)APtr);
	//Invert2(APtr, BPtr);
	/*
	srand(time(NULL));
	#pragma omp parallel for num_threads(4)
	for (i=0;i<1000;i++) {
		int r=rand()%651;
		*(BPtr+i)=r;
	}*/
	BPtr=inverse3x3(APtr);

}
