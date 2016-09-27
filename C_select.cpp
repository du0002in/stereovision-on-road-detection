// mex file to calculate gradient
// input:   mxn gray image
// output:  gradient along x (mxn) and y (mxn)
// http://blog.teamleadnet.com/2012/07/quick-select-algorithm-find-kth-element.html

#include "mex.h"
//#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *arr, *kPtr, *out;
	int from=0, to, k;
	arr=mxGetPr(prhs[0]);
	to=mxGetNumberOfElements(prhs[0])-1;
	kPtr=mxGetPr(prhs[1]);
	k=(int)(*kPtr);

	plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
	out=mxGetPr(plhs[0]);

	while (from<to) {
		int r = from, w = to;
		int mid = arr[(r + w) / 2];
 
		// stop if the reader and writer meets
		while (r < w) {
 
			if (arr[r] >= mid) { // put the large values at the end
			int tmp = arr[w];
			arr[w] = arr[r];
			arr[r] = tmp;
			w--;
			} else { // the value is smaller than the pivot, skip
				r++;
			}
		}
 
		// if we stepped up (r++) we need to step one down
		if (arr[r] > mid)
			r--;
 
		// the r pointer is on the end of the first k elements
		if (k <= r) {
			to = r;
		} else {
			from = r + 1;
		}
	}
	*out=*(arr+k);

}





/*



    double *ImgPtr, *Grad_xPtr, *Grad_yPtr,*Grad_xPtr1,*Grad_xPtr2;
    mwSize mrows, curcols, precols, nextcols;
    mwSize ncols;
    int i, j, nt;
    clock_t start, end;
	int nthreads, tid, procs, maxt, inpar, dynamic, nested;
    ImgPtr=mxGetPr(prhs[0]);
    mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
    plhs[0]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	plhs[1]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
	plhs[2]=mxCreateDoubleMatrix(mrows, ncols, mxREAL);
    plhs[3]=mxCreateDoubleMatrix(3, 1, mxREAL);
    Grad_xPtr=mxGetPr(plhs[0]);
	Grad_xPtr1=mxGetPr(plhs[1]);
	Grad_xPtr2=mxGetPr(plhs[2]);
    Grad_yPtr=mxGetPr(plhs[3]);

	start=clock();
	#pragma omp parallel private(nthreads, tid)   
	{

// Obtain thread number    
tid = omp_get_thread_num();

// Only master thread does this    
if (tid == 0) 
{
printf("Thread %d getting environment info...\n", tid);

// Get environment information 
procs = omp_get_num_procs();
nthreads = omp_get_num_threads();
maxt = omp_get_max_threads();
inpar = omp_in_parallel();
dynamic = omp_get_dynamic();
nested = omp_get_nested();

// Print environment information 
printf("Number of processors = %d\n", procs);
printf("Number of threads = %d\n", nthreads);
printf("Max threads = %d\n", maxt);
printf("In parallel? = %d\n", inpar);
printf("Dynamic threads enabled? = %d\n", dynamic);
printf("Nested parallelism supported? = %d\n", nested);  
}
}


	/*
	#pragma omp parallel for
	{
		nt=omp_get_num_threads();
		*Grad_yPtr=nt;
		for (i=0;i<mrows*ncols;i++)
		*(Grad_xPtr+i)=sin(*(ImgPtr+i));
	}


	/*
	for (i=0;i<ncols; i++) 
	{
		for (j=0;j<mrows;j++)
		*(Grad_xPtr+i*mrows+j)=sin(*(ImgPtr+i*mrows+j)**(ImgPtr+i*mrows+j)**(ImgPtr+i*mrows+j))/(cos(*(ImgPtr+i*mrows+j))+1);
	}
	end=clock();
	*Grad_yPtr=(double)(end-start)/CLOCKS_PER_SEC;

	start=clock();
	#pragma omp parallel for private(i)
	{
		for (i=0;i<ncols; i++) 
		{
			for (j=0;j<mrows;j++)
				*(Grad_xPtr1+i*mrows+j)=sin(*(ImgPtr+i*mrows+j)**(ImgPtr+i*mrows+j)**(ImgPtr+i*mrows+j))/(cos(*(ImgPtr+i*mrows+j))+1);
		}
	}
	end=clock();
	*(Grad_yPtr+1)=(double)(end-start)/CLOCKS_PER_SEC;

	start=clock();
	#pragma omp parallel for
	{
		for (i=0;i<ncols; i++) 
		{
			for (j=0;j<mrows;j++)
				*(Grad_xPtr1+i*mrows+j)=sin(*(ImgPtr+i*mrows+j)**(ImgPtr+i*mrows+j)**(ImgPtr+i*mrows+j))/(cos(*(ImgPtr+i*mrows+j))+1);
		}
	}
	
	end=clock();
	*(Grad_yPtr+2)=(double)(end-start)/CLOCKS_PER_SEC;
}
*/

    //Grad_xPtr calculation
/*	
#pragma omp parallel for private(i)
	{
    for (i=1;i<(ncols-1);i++)
    {
        curcols=i*mrows; precols=(i-1)*mrows; nextcols=(i+1)*mrows;
        for (j=0;j<mrows;j++)
            *(Grad_xPtr+curcols+j)=0.5*(*(ImgPtr+nextcols+j)-*(ImgPtr+precols+j));
    }
	}
    curcols=(ncols-1)*mrows; precols=(ncols-2)*mrows;
	
    for (j=0;j<mrows;j++)
    {
        *(Grad_xPtr+j)=*(ImgPtr+mrows+j)-*(ImgPtr+j);
        *(Grad_xPtr+curcols+j)=*(ImgPtr+curcols+j)-*(ImgPtr+precols+j);
    }

    //Grad_yPtr calculation
	#pragma omp parallel for private(i) 
	{
    for (i=0;i<ncols;i++)
    {
        curcols=i*mrows;
        for (j=1;j<mrows-1;j++)
            *(Grad_yPtr+curcols+j)=0.5*(*(ImgPtr+curcols+j+1)-*(ImgPtr+curcols+j-1));
    }
	}
    for (i=0;i<ncols;i++)
    {
        curcols=i*mrows; nextcols=(i+1)*mrows;
        *(Grad_yPtr+curcols)=*(ImgPtr+curcols+1)-*(ImgPtr+curcols);
        *(Grad_yPtr+nextcols-1)=*(ImgPtr+nextcols-1)-*(ImgPtr+nextcols-2);
    }
}*/