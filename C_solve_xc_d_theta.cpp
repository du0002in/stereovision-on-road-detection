// mex file to calculate divergence
// input:   mxn smoothed gray image
// output:  mxn matrix, value at (i,j) refers to divergence at that pixel
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include "mex.h"
#include <omp.h>
//#include <math.h>
//#include "solver.h"

static int isZero(double x)
{
return x > -0.000001 && x < 0.000001;
}

void solveCubic(double c[4], double s[3], double *num_root)
{
int	i, num;
double	sub,
	A, B, C,
	sq_A, p, q,
	cb_p, D;

// normalize the equation:x ^ 3 + Ax ^ 2 + Bx  + C = 0
A = c[2] / c[3];
B = c[1] / c[3];
C = c[0] / c[3];

// substitute x = y - A / 3 to eliminate the quadric term: x^3 + px + q = 0

sq_A = A * A;
p = 1.0/3.0 * (-1.0/3.0 * sq_A + B);
q = 1.0/2.0 * (2.0/27.0 * A *sq_A - 1.0/3.0 * A * B + C);

// use Cardano's formula

cb_p = p * p * p;
D = q * q + cb_p;

if (isZero(D))
    {
    if (isZero(q))
	{
	// one triple solution
	s[0] = 0.0;
	num = 1;
	}
    else
	{
	// one single and one double solution
	double u = pow(-q,1.0/3.0);
	s[0] = 2.0 * u;
	s[1] = - u;
	num = 2;
	}
    }
else
    if (D < 0.0)
	{
	// casus irreductibilis: three real solutions
	double phi = 1.0/3.0 * acos(-q / sqrt(-cb_p));
	double t = 2.0 * sqrt(-p);
	s[0] = t * cos(phi);
	s[1] = -t * cos(phi + M_PI / 3.0);
	s[2] = -t * cos(phi - M_PI / 3.0);
	num = 3;
	}
    else
	{
	// one real solution
	double sqrt_D = sqrt(D);
	double u = pow(sqrt_D + fabs(q),1.0/3.0);
	if (q > 0.0)
	    s[0] = - u + p / u ;
	else
	    s[0] = u - p / u;
	num = 1;
	}

// resubstitute
sub = 1.0 / 3.0 * A;
for (i = 0; i < num; i++)
    s[i] -= sub;
*num_root=num;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *polyPtr, *ParticlePtr, *xcPtr, *d_thetaPtr;
	mwSize P_No;
	int kk;
	polyPtr=mxGetPr(prhs[0]);
	double a=*polyPtr; 
	double b=*(polyPtr+1); 
	double c=*(polyPtr+2);

	ParticlePtr=mxGetPr(prhs[1]);
	P_No=mxGetN(prhs[1]);

	plhs[0]=mxCreateDoubleMatrix(1, P_No, mxREAL);
	plhs[1]=mxCreateDoubleMatrix(1, P_No, mxREAL);
	xcPtr=mxGetPr(plhs[0]);
	d_thetaPtr=mxGetPr(plhs[1]);

	if (a<-0.0000001 || a>0.0000001) {
	#pragma omp parallel for num_threads(4)
	for (kk=0; kk<P_No; kk++)
	{
		double co[4], s[3], num_root;
		double x0=*(ParticlePtr+kk*3),
			   z0=*(ParticlePtr+kk*3+1),
			   theta0=*(ParticlePtr+kk*3+2);
		co[0]=b*c-b*x0-z0;
		co[1]=b*b+2*a*(c-x0)+1;
		co[2]=3*a*b;
		co[3]=2*a*a;
		solveCubic(co,s,&num_root);
		if (num_root>=2)
		{
			double d_sqr_pre_inv=0;
			for (int i=0; i<num_root; i++)
			{
				double d_sqr_inv=1.0/(pow(a*s[i]*s[i]+b*s[i]+c-x0,2)+pow(s[i]-z0,2));
				if (d_sqr_inv>d_sqr_pre_inv)
				{
						s[0]=s[i];
						d_sqr_pre_inv=d_sqr_inv;
				}
			}
		}
		double z_tangent=s[0];
		double x_tangent=a*z_tangent*z_tangent+b*z_tangent+c;
		double d_theta=atan(1/(2*a*z_tangent+b));
		if (d_theta<0)
		    d_theta=d_theta+M_PI;
		*(xcPtr+kk)=-sin(d_theta)*(x0-x_tangent)+cos(d_theta)*(z0-z_tangent);
		*(d_thetaPtr+kk)=(-d_theta+theta0);
	}
	}
	else {
	#pragma omp parallel for num_threads(4)
	for (kk=0; kk<P_No; kk++)
	{
		double co[2], s[1], num_root;
		double x0=*(ParticlePtr+kk*3),
			   z0=*(ParticlePtr+kk*3+1),
			   theta0=*(ParticlePtr+kk*3+2);
		co[0]=b*c-b*x0-z0;
		co[1]=b*b+1;
		//co[2]=3*a*b;
		//co[3]=2*a*a;
		s[0]=-co[0]/co[1];
		double z_tangent=s[0];
		double x_tangent=a*z_tangent*z_tangent+b*z_tangent+c;
		double d_theta=atan(1/(2*a*z_tangent+b));
		if (d_theta<0)
		    d_theta=d_theta+M_PI;
		*(xcPtr+kk)=-sin(d_theta)*(x0-x_tangent)+cos(d_theta)*(z0-z_tangent);
		*(d_thetaPtr+kk)=(-d_theta+theta0);
	}
	}

}