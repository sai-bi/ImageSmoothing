/*				SplitBregmanROF.c (MEX version) by Tom Goldstein
 *   This code performs isotropic ROF denoising using the "Split Bregman" algorithm.
 * This version of the code has a "mex" interface, and should be compiled and called
 * through MATLAB.
 *
 *DISCLAIMER:  This code is for academic (non-commercial) use only.  Also, this code 
 *comes with absolutely NO warranty of any kind: I do my best to write reliable codes, 
 *but I take no responsibility for the reliability of the results.
 *
 *                      HOW TO COMPILE THIS CODE
 *   To compile this code, open a MATLAB terminal, and change the current directory to
 *the folder where this "c" file is contained.  Then, enter this command:
 *    >>  mex splitBregmanROF.c
 *This file has been tested under windows using lcc, and under SUSE Unix using gcc.
 *Once the file is compiled, the command "splitBregmanROF" can be used just like any
 *other MATLAB m-file.
 *
 *                      HOW TO USE THIS CODE
 * An image is denoised using the following command
 * 
 *   SplitBregmanROF(image,mu,tol);
 *
 * where:
 *   - "image" is a 2d array containing the noisy image.
 *   - "mu" is the weighting parameter for the fidelity term 
 *            (usually a value between 0.01 and 0.5 works well for images with 
 *                           pixels on the 0-255 scale).
 *   - "tol" is the stopping tolerance for the iteration.  "tol"=0.001 is reasonable for
 *            most applications.
 *  
 */

#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
typedef double num;


/*A method for isotropic TV*/
void rof_iso(num** u, num** f, num** x, num** y, num** bx, num** by ,
					 double mu, double lambda, int nGS, int nBreg, int width, int height);
	
/*A method for Anisotropic TV*/
void rof_an(num** u, num** f, num** x, num** y, num** bx, num** by ,
					 double mu, double lambda, int nGS, int nBreg, int width, int height);
		
	/*****************Minimization Methods*****************/
void gs_an(num** u, num** f, num** x, num** y, num** bx, num** by , double mu, double lambda, int width, int height, int iter);
void gs_iso(num** u, num** f, num** x, num** y, num** bx, num** by , double mu, double lambda, int width, int height, int iter);
		
	/******************Relaxation Methods*****************/
void gsU(num** u, num** f, num** x, num** y, num** bx, num** by , double mu, double lambda, int width, int height);
void gsX(num** u, num** x, num** bx , double lambda, int width, int height);
void gsY(num** u, num** y, num** by , double lambda, int width, int height);
void gsSpace(num** u, num** x, num** y, num** bx, num** by, double lambda, int width, int height);
	 
	/************************Bregman***********************/
void bregmanX(num** x,num** u, num** bx, int width, int height);
void bregmanY(num** y,num** u, num** by, int width, int height);

/**********************Memory************/
	
num** newMatrix(int rows, int cols);
void deleteMatrix(num ** a);
double** get2dArray(const mxArray *mx, int isCopy);
double copy(double** source, double** dest, int rows, int cols);
	
/***********************The MEX Interface to rof_iso()***************/
	
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
		/* get the size of the image*/
	int rows = mxGetN(prhs[0]);
	int cols = mxGetM(prhs[0]);
	
		/* get the fidelity and convergence parameters*/
	double mu =  (double)(mxGetScalar(prhs[1]));
	double lambda = 2*mu;
	double tol = (double)(mxGetScalar(prhs[2]));
	
	/* get the image, and declare memory to hold the auxillary variables*/
	double **f = get2dArray(prhs[0],0);
    double **u = newMatrix(rows,cols);
	double **x = newMatrix(rows-1,cols);
	double **y = newMatrix(rows,cols-1);
	double **bx = newMatrix(rows-1,cols);
	double **by = newMatrix(rows,cols-1);
	
    double** uOld;
    double *outArray;
    double diff;
    int count;
    int i,j;
    
     /***********Check Conditions******/  
    if (nrhs != 3) {mexErrMsgTxt("Three input arguments required.");} 
    if (nlhs > 1){mexErrMsgTxt("Too many output arguments.");}
    if (!(mxIsDouble(prhs[0]))) {mexErrMsgTxt("Input array must be of type double.");}
   
     /* Use a copy of the image as an initial guess*/
		for(i=0;i<rows;i++){
		    for(j=0;j<cols;j++){
		        u[i][j] = f[i][j];
		    }
		}
    
	/* denoise the image*/
	
	uOld = newMatrix(rows,cols);
    count=0;
	do{
		rof_iso(u,f,x,y,bx,by,mu,lambda,1,5,rows,cols);
		diff = copy(u,uOld,rows,cols);
        count++;
	}while( (diff>tol && count<1000) || count<5 );
	
	/* copy denoised image to output vector*/
	plhs[0] = mxCreateDoubleMatrix(cols, rows, mxREAL); /*mxReal is our data-type*/
	outArray = mxGetPr(plhs[0]);
	
	for(i=0;i<rows;i++){
	    for(j=0;j<cols;j++){
	        outArray[(i*cols)+j] = u[i][j];
	    }
	}
	
	/* Free memory */
    deleteMatrix(u);
	deleteMatrix(x);
	deleteMatrix(y);
	deleteMatrix(bx);
	deleteMatrix(by);
    deleteMatrix(uOld);
	
    return;
}

double** get2dArray(const mxArray *mx, int isCopy){
	double* oned = mxGetPr(mx);
	int rowLen = mxGetN(mx);
	int colLen = mxGetM(mx);
	double** rval = (double**) malloc(rowLen*sizeof(double*));
    int r;
	if(isCopy){
		double *copy = (double*)malloc(rowLen*colLen*sizeof(double));
		int i, sent = rowLen*colLen;
		for(i=0;i<sent;i++)
			copy[i]=oned[i];
		oned=copy;
	}
	
	for(r=0;r<rowLen;r++)
		rval[r] = &oned[colLen*r];
	return rval;
}







/*                IMPLEMENTATION BELOW THIS LINE                         */

/******************Isotropic TV**************/

void rof_iso(num** u, num** f, num** x, num** y, num** bx, num** by ,
								 double mu, double lambda, int nGS, int nBreg, int width, int height){
		int breg;
		for(breg=0;breg<nBreg;breg++){
				gs_iso(u,f,x,y,bx,by,mu,lambda,width, height,nGS);
				bregmanX(x,u,bx,width,height);
				bregmanY(y,u,by,width,height);
		}
}	


void gs_iso(num** u, num** f, num** x, num** y, num** bx, num** by , double mu, double lambda, int width, int height, int iter){
		int j;
		for(j=0;j<iter;j++){
			gsU(u,f,x,y,bx,by,mu,lambda,width,height);
			gsSpace(u,x,y,bx,by,lambda,width,height);
			
		}
}
	

/******************Anisotropic TV**************/
void rof_an(num** u, num** f, num** x, num** y, num** bx, num** by ,
								 double mu, double lambda, int nGS, int nBreg, int width, int height){
		int breg;
		for(breg=0;breg<nBreg;breg++){
				gs_an(u,f,x,y,bx,by,mu,lambda,width, height,nGS);
				bregmanX(x,u,bx,width,height);
				bregmanY(y,u,by,width,height);
		}
}	


void gs_an(num** u, num** f, num** x, num** y, num** bx, num** by , double mu, double lambda, int width, int height, int iter){
		int j;
		for(j=0;j<iter;j++){
			gsU(u,f,x,y,bx,by,mu,lambda,width,height);
			gsX(u,x,bx,lambda,width,height);
			gsY(u,y,by,lambda,width,height);
		}
}
	
	



/****Relaxation operators****/

void gsU(num** u, num** f, num** x, num** y, num** bx, num** by , double mu, double lambda, int width, int height){
	int w,h;
	double sum;
	int wm1,hm1;
	double normConst = 1.0/(mu+4*lambda);
	int wSent = width-1, hSent = height-1;
	for(w=1;w<wSent;w++){		/* do the central pixels*/
		wm1 = w-1;
		for(h=1;h<hSent;h++){
			hm1 = h-1;
			sum = x[wm1][h] - x[w][h]+y[w][hm1] - y[w][h]
						-bx[wm1][h] + bx[w][h]-by[w][hm1] + by[w][h];
			sum+=(u[w+1][h]+u[wm1][h]+u[w][h+1]+u[w][hm1]);
			sum*=lambda;
			sum+=mu*f[w][h];
			sum*=normConst;
			u[w][h] = sum;
		}
	}
		w=0;				/* do the left pixels*/
		for(h=1;h<hSent;h++){
			sum = - x[w][h]+y[w][h-1] - y[w][h]
						+ bx[w][h]-by[w][h-1] + by[w][h];
			sum+=(u[w+1][h]+u[w][h+1]+u[w][h-1]);
			sum*=lambda;
			sum+=mu*f[w][h];
			sum/=mu+3*lambda;
			u[w][h] = sum;
		}
		w = width-1;		/* do the right pixels*/
			for(h=1;h<hSent;h++){
				sum = x[w-1][h] +y[w][h-1] - y[w][h]
								-bx[w-1][h] -by[w][h-1] + by[w][h];
				sum+=u[w-1][h]+u[w][h+1]+u[w][h-1];
				sum*=lambda;
				sum+=mu*f[w][h];
				sum/=mu+3*lambda;
				u[w][h] = sum;
			}

		h=0;
		for(w=1;w<wSent;w++){		/* do the top pixels*/
				sum = x[w-1][h] - x[w][h] - y[w][h]
								-bx[w-1][h] + bx[w][h] + by[w][h];
				sum+=u[w+1][h]+u[w-1][h]+u[w][h+1];
				sum*=lambda;
				sum+=mu*f[w][h];
				sum/=mu+3*lambda;
				u[w][h] = sum;
			}
		h = height-1;
		for(w=1;w<wSent;w++){		/* do the bottom pixels*/
				sum = x[w-1][h] - x[w][h]+y[w][h-1] 
								-bx[w-1][h] + bx[w][h]-by[w][h-1];
				sum+=u[w+1][h]+u[w-1][h]+u[w][h-1];
				sum*=lambda;
				sum+=mu*f[w][h];
				sum/=mu+3*lambda;
				u[w][h] = sum;
			}	
			/* do the top left pixel*/		
			w=h=0;
				sum =  - x[w][h] - y[w][h]
								+ bx[w][h] + by[w][h];
				sum+=u[w+1][h]+u[w][h+1];
				sum*=lambda;
				sum+=mu*f[w][h];
				sum/=mu+2*lambda;
				u[w][h] = sum;
				/* do the top right pixel*/		
			w=width-1; h=0;
				sum = x[w-1][h]  - y[w][h]
								-bx[w-1][h] + by[w][h];
				sum+=u[w-1][h]+u[w][h+1];
				sum*=lambda;
				sum+=mu*f[w][h];
				sum/=mu+2*lambda;
				u[w][h] = sum;
				/* do the bottom left pixel*/		
			w=0; h=height-1;
				sum =  - x[w][h]+y[w][h-1] 
								+ bx[w][h]-by[w][h-1] ;
				
				sum+=u[w+1][h]+u[w][h-1];
				sum*=lambda;
				sum+=mu*f[w][h];
				sum/=mu+2*lambda;
				u[w][h] = sum;
				/* do the bottom right pixel*/		
			w=width-1; h=height-1;
				sum = x[w-1][h]+y[w][h-1] 
								-bx[w-1][h] -by[w][h-1];	
				sum+=u[w-1][h]+u[w][h-1];
				sum*=lambda;
				sum+=mu*f[w][h];
				sum/=mu+2*lambda;
				u[w][h] = sum;
}


void gsSpace(num** u, num** x, num** y, num** bx, num** by, double lambda, int width, int height){
	int w,h;
	num a,b,s;
	num flux = 1.0/lambda;
	num mflux = -1.0/lambda;
	num flux2 = flux*flux;
	num *uw,*uwp1,*bxw,*byw,*xw,*yw;
    num base;
	for(w=0;w<width-1;w++){	
		uw = u[w];uwp1=u[w+1];bxw=bx[w];byw=by[w];xw=x[w];yw=y[w];
		for(h=0;h<height-1;h++){
			a =  uwp1[h]-uw[h]+bxw[h];
			b =  uw[h+1]-uw[h]+byw[h];
			s = a*a+b*b;
			if(s<flux2){xw[h]=0;yw[h]=0;continue;}
			s = sqrt(s);
			s=(s-flux)/s;
			xw[h] = s*a;
			yw[h] = s*b;
		}
	}		

	h = height-1;
	for(w=0;w<width-1;w++){	
			base =  u[w+1][h]-u[w][h]+bx[w][h];
			if(base>flux) {x[w][h] = base-flux; continue;}
			if(base<mflux){x[w][h] = base+flux; continue;}
			x[w][h] = 0;
	}
	w = width-1;
	for(h=0;h<height-1;h++){	
		base =  u[w][h+1]-u[w][h]+by[w][h];
		if(base>flux) {y[w][h] = base-flux; continue;}
		if(base<mflux){y[w][h] = base+flux; continue;}
		y[w][h] = 0;
	}
}


void gsX(num** u, num** x, num** bx , double lambda, int width, int height){
	int w,h;
	double base;				
	const double flux = 1.0/lambda;
	const double mflux = -1.0/lambda;
	num* uwp1;
	num* uw;
	num* bxw;
	num* xw;
    width = width-1;
	for(w=0;w<width;w++){
		uwp1 = u[w+1];
		uw = u[w];
		bxw = bx[w];
		xw = x[w];
		for(h=0;h<height;h++){
			base = uwp1[h]-uw[h]+bxw[h];
			if(base>flux) {xw[h] = base-flux; continue;}
			if(base<mflux){xw[h] = base+flux; continue;}
			xw[h] = 0;
		}
	}
}

void gsY(num** u, num** y, num** by , double lambda, int width, int height){
	int w,h;
	double base;				
	const double flux = 1.0/lambda;
	const double mflux = -1.0/lambda;
	num* uw;
	num* yw;
	num* bw;
    height = height-1;
	for(w=0;w<width;w++){
		uw = u[w];
		yw = y[w];
		bw = by[w];
		for(h=0;h<height;h++){
			base = uw[h+1]-uw[h]+bw[h];
			if(base>flux) {yw[h] = base-flux; continue;}
			if(base<mflux){yw[h] = base+flux; continue;}
			yw[h] = 0;
		}
	}
}

void bregmanX(num** x,num** u, num** bx, int width, int height){
		int w,h;
		double d;
		num* uwp1,*uw,*bxw,*xw;
		for(w=0;w<width-1;w++){
			uwp1=u[w+1];uw=u[w];bxw=bx[w];xw=x[w];
			for(h=0;h<height;h++){
				d = uwp1[h]-uw[h];
				bxw[h]+= d-xw[h];		
			}
		}
	}


void bregmanY(num** y,num** u, num** by, int width, int height){
		int w,h;
		double d;
		int hSent = height-1;
		num* uw,*byw,*yw;
		for(w=0;w<width;w++){
			uw=u[w];byw=by[w];yw=y[w];
			for(h=0;h<hSent;h++){
				d = uw[h+1]-uw[h];
				byw[h]+= d-yw[h];
			}
		}
	}
	
/************************memory****************/

double copy(double** source, double** dest, int rows, int cols){
	int r,c;
	double temp,sumDiff, sum;
	for(r=0;r<rows;r++)
		for(c=0;c<cols;c++){
			temp = dest[r][c];
			sum+=temp*temp;
			temp -= source[r][c];
			sumDiff +=temp*temp;
			
			dest[r][c]=source[r][c];
		
		}
	return sqrt(sumDiff/sum);
}



num** newMatrix(int rows, int cols){
		num* a = (num*) malloc(rows*cols*sizeof(num));
		num** rval = (num**) malloc(rows*sizeof(num*));
		int j,g;
        rval[0] = a;
		for(j=1;j<rows;j++)
			rval[j] = &a[j*cols];

		for(j=0;j<rows;j++)
			for(g=0;g<cols;g++)
				rval[j][g] = 0;
		return rval;
}
void deleteMatrix(num ** a){
	free(a[0]);
	free(a);
}
