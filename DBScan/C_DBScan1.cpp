#include "mex.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <omp.h>
//#include "type.hpp"

void GetDist(double *x, double *y, int noPts, double *distPtr)
{
	int i;
	#pragma omp parallel for num_threads(4)
	for (i=0; i<noPts; i++)
	{
		for (int j=0; j<i; j++)
		{
			distPtr[i+j*noPts]=pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
			distPtr[j+i*noPts]=distPtr[i+j*noPts];
		}
	}
}

std::vector<int> regionQuery(int ind0, double eps, int noPts, double *distPtr)
{
	double dist;
	std::vector<int> retKeys;
	ind0=ind0*noPts;
	for(int i = 0; i< noPts; i++)
	{
		//dist = sqrt(pow((x0-x[i]),2)+pow((y0-y[i]),2));
		dist=distPtr[i+ind0];
		if(dist <= eps && dist != 0.0f)
		{
			retKeys.push_back(i);
		}
	}
	return retKeys;
}

std::vector<int> DBScan(double *x, double *y, double eps, int minPts, int noPts)
{
	std::vector<int> cluster;
	std::vector<int> noise;
	std::vector<bool> visited;
	std::vector<int> neighborPts;
	std::vector<int> neighborPts_;
	int c=0;

	for(int k = 0; k < noPts; k++)
	{
		cluster.push_back(0);
		visited.push_back(false);
	}

	mxArray *distMx;
	double *distPtr;
	distMx=mxCreateDoubleMatrix(noPts,noPts,mxREAL);
	distPtr=mxGetPr(distMx);
	GetDist(x,y,noPts,distPtr);
	eps=eps*eps;

	for(int i=0; i<noPts; i++)
	{
		if(!visited[i])
		{
			visited[i] = true;
			neighborPts = regionQuery(i,eps, noPts, distPtr);
			if(neighborPts.size() < minPts)
				noise.push_back(i);
			else
			{
				c++;
				cluster[i]=c;
				for(int j = 0; j < neighborPts.size(); j++)
				{
					if(!visited[neighborPts[j]])
					{
						visited[neighborPts[j]] = true;
						neighborPts_ = regionQuery(neighborPts[j],eps, noPts, distPtr);
						if(neighborPts_.size() >= minPts)
							neighborPts.insert(neighborPts.end(),neighborPts_.begin(),neighborPts_.end());
					}
					if(cluster[neighborPts[j]]==0)
						cluster[neighborPts[j]]=c;
				}
			}
		}
	}
	mxDestroyArray(distMx);
	return cluster;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *xPtr, *yPtr, *eps, *minPts;
	double *outPtr;
	int noPts;

	xPtr=mxGetPr(prhs[0]);
	yPtr=mxGetPr(prhs[1]);
	eps=mxGetPr(prhs[2]);
	minPts=mxGetPr(prhs[3]);
	noPts=mxGetNumberOfElements(prhs[0]);
	plhs[0]=mxCreateDoubleMatrix(1,noPts,mxREAL);
	outPtr=mxGetPr(plhs[0]);

	std::vector<int> cluster;
	cluster=DBScan(xPtr, yPtr, *eps, (double)(*minPts), noPts);
	for (int i=0; i<noPts; i++)
		outPtr[i]=cluster[i];
}