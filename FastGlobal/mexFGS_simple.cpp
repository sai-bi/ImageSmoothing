#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mex.h"

#define SQ(x) ((x)*(x))

int W, H; // image width, height
int nChannels, nChannels_guide;
double *BLFKernelI; // Kernel LUT

// Main functions
void prepareBLFKernel(double sigma);
void FGS_simple(double ***image, double ***joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation);
void solve_tridiagonal_in_place_destructive(double x[], const size_t N, const double a[], const double b[], double c[]);

// Memory management
double *** memAllocDouble3(int n,int r,int c);
double** memAllocDouble2(int r,int c);
void memFreeDouble3(double ***p);
void memFreeDouble2(double **p);


// Build LUT for bilateral kernel weight
void prepareBLFKernel(double sigma)
{
	const int MaxSizeOfFilterI = 195075;
	BLFKernelI = (double *)malloc(sizeof(double)*MaxSizeOfFilterI);

	for(int m=0;m<MaxSizeOfFilterI;m++)
		BLFKernelI[m] = exp( -sqrt((double)m)/(sigma) ); // Kernel LUT
}

// mex function call:
// x = mexFGS(input_image, guidance_image = NULL, sigma, lambda, fgs_iteration = 3, fgs_attenuation = 4);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs < 6) { mexErrMsgTxt("FGS must be called with 6 arguments. Please see FGS.m for more details."); }
  
	const mxArray *img = prhs[0], *imgGuide = prhs[1]; // input images
	mxArray *image_result; // output array
	
	// image resolution
	W = mxGetDimensions(img)[1];
	H = mxGetDimensions(img)[0];
    if(mxGetNumberOfDimensions(img) > 2)
        nChannels = mxGetDimensions(img)[2];
    else
        nChannels = 1;
    
	// FGS parameters
	double sigma = mxGetScalar(prhs[2]);
	double lambda = mxGetScalar(prhs[3]);
	int solver_iteration = (int)mxGetScalar(prhs[4]);
	double solver_attenuation = mxGetScalar(prhs[5]);

	mexPrintf("Image resolution: %d x %d x %d\n", W, H, nChannels);
	mexPrintf("Parameters:\n");
	mexPrintf("    Sigma = %f\n", sigma);
	mexPrintf("    Lambda = %f\n", lambda);
	mexPrintf("    Iteration = %d\n", solver_iteration);
	mexPrintf("    Attenuation = %f\n", solver_attenuation);
    
	// Image buffer preperation
	double*** image_filtered = memAllocDouble3(H, W, nChannels);
	double* ptr_image = (double*)mxGetPr(img);
    double* ptr_image_array = image_filtered[0][0];
	for(int y=0;y<H;y++) {
        for(int x=0;x<W;x++) {
            for(int c=0;c<nChannels;c++) {
                //image_filtered[y][x][c] = (double)ptr_image[x*H + y + c*(H*W)];
                ptr_image_array[y*W*nChannels + x*nChannels + c] = (double)ptr_image[x*H + y + c*(H*W)];
            }
        }
    }

 	double*** image_guidance = NULL;
 	if(!mxIsEmpty(imgGuide)) {        
        if(mxGetNumberOfDimensions(imgGuide) > 2)
            nChannels_guide = mxGetDimensions(imgGuide)[2];
        else
            nChannels_guide = 1;

        mexPrintf("Joint filtering mode: %d x %d x %d\n", W, H, nChannels_guide);        
 		image_guidance = memAllocDouble3(H, W, nChannels_guide);
        
 		double* ptr_guidance = (double*)mxGetPr(imgGuide);
        double* ptr_guidance_array = image_guidance[0][0];
 		for(int y=0;y<H;y++) for(int x=0;x<W;x++) for(int c=0;c<nChannels_guide;c++) //image_guidance[y][x][c] = (double)ptr_guidance[x*H + y + c*(H*W)];
            ptr_guidance_array[y*W*nChannels_guide + x*nChannels_guide + c] = (double)ptr_guidance[x*H + y + c*(H*W)];
 	} else {
        nChannels_guide = nChannels;
    }
    

	// run FGS
	sigma *= 255.0;
    
    clock_t m_begin = clock(); // time measurement;
	FGS_simple(image_filtered, image_guidance, sigma, lambda, solver_iteration, solver_attenuation);
    mexPrintf("Elapsed time is %f seconds.\n", double(clock()-m_begin)/CLOCKS_PER_SEC);
	   
	// output
	image_result = plhs[0] = mxDuplicateArray(img);
	double* ptr_output = mxGetPr(image_result);
	for(int y=0;y<H;y++) for(int x=0;x<W;x++) for(int c=0;c<nChannels;c++) ptr_output[x*H + y + c*(H*W)] = image_filtered[y][x][c];

	memFreeDouble3(image_filtered);
	if(image_guidance) memFreeDouble3(image_guidance);
}

void FGS_simple(double ***image, double ***joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation)
{
	int color_diff;

	int width = W;
	int height = H;

    if(joint_image == NULL) joint_image = image;
    
	double *a_vec = (double *)malloc(sizeof(double)*width);
	double *b_vec = (double *)malloc(sizeof(double)*width);
	double *c_vec = (double *)malloc(sizeof(double)*width);	
	double *x_vec = (double *)malloc(sizeof(double)*width);
	double *c_ori_vec = (double *)malloc(sizeof(double)*width);

	double *a2_vec = (double *)malloc(sizeof(double)*height);
	double *b2_vec = (double *)malloc(sizeof(double)*height);
	double *c2_vec = (double *)malloc(sizeof(double)*height);
	double *x2_vec = (double *)malloc(sizeof(double)*height);
	double *c2_ori_vec = (double *)malloc(sizeof(double)*height);

	prepareBLFKernel(sigma);
	
	//Variation of lambda (NEW)
	double lambda_in = 1.5*lambda*pow(4.0,solver_iteration-1)/(pow(4.0,solver_iteration)-1.0);
	for(int iter=0;iter<solver_iteration;iter++)
	{
		//for each row
		for(int i=0;i<height;i++)
		{
			memset(a_vec, 0, sizeof(double)*width);
			memset(b_vec, 0, sizeof(double)*width);
			memset(c_vec, 0, sizeof(double)*width);
			memset(c_ori_vec, 0, sizeof(double)*width);
			memset(x_vec, 0, sizeof(double)*width);
			for(int j=1;j<width;j++)
			{
                int color_diff = 0;
                // compute bilateral weight for all channels
                for(int c=0;c<nChannels_guide;c++)
                    color_diff += SQ(joint_image[i][j][c] - joint_image[i][j-1][c]);
				
				a_vec[j] = -lambda_in*BLFKernelI[color_diff];		//WLS
			}
			for(int j=0;j<width-1;j++)	c_ori_vec[j] = a_vec[j+1];
			for(int j=0;j<width;j++)	b_vec[j] = 1.f - a_vec[j] - c_ori_vec[j];		//WLS
			
            for(int c=0;c<nChannels;c++)
            {
                memcpy(c_vec, c_ori_vec, sizeof(double)*width);
                for(int j=0;j<width;j++)	x_vec[j] = image[i][j][c];			
                solve_tridiagonal_in_place_destructive(x_vec, width, a_vec, b_vec, c_vec);
                for(int j=0;j<width;j++)	image[i][j][c] = x_vec[j];                
            }
		}

		//for each column
		for(int j=0;j<width;j++)
		{
			memset(a2_vec, 0, sizeof(double)*height);
			memset(b2_vec, 0, sizeof(double)*height);
			memset(c2_vec, 0, sizeof(double)*height);
			memset(c2_ori_vec, 0, sizeof(double)*height);
			memset(x2_vec, 0, sizeof(double)*height);
			for(int i=1;i<height;i++)
			{
                int color_diff = 0;
                // compute bilateral weight for all channels
                for(int c=0;c<nChannels_guide;c++)
                    color_diff += SQ(joint_image[i][j][c] - joint_image[i-1][j][c]);

                a2_vec[i] = -lambda_in*BLFKernelI[color_diff];		//WLS
			}
			for(int i=0;i<height-1;i++)
				c2_ori_vec[i] = a2_vec[i+1];
			for(int i=0;i<height;i++)
				b2_vec[i] = 1.f - a2_vec[i] - c2_ori_vec[i];		//WLS

            for(int c=0;c<nChannels;c++)
            {
                memcpy(c2_vec, c2_ori_vec, sizeof(double)*height);
                for(int i=0;i<height;i++)	x2_vec[i] = image[i][j][c];
                solve_tridiagonal_in_place_destructive(x2_vec, height, a2_vec, b2_vec, c2_vec);
                for(int i=0;i<height;i++)	image[i][j][c] = x2_vec[i];               
            }
		}

		//Variation of lambda (NEW)
		lambda_in /= solver_attenuation;
	}	//iter	
	
	free(a_vec);
	free(b_vec);
	free(c_vec);
	free(x_vec);
	free(c_ori_vec);

	free(a2_vec);
	free(b2_vec);
	free(c2_vec);
	free(x2_vec);
	free(c2_ori_vec);
}

void solve_tridiagonal_in_place_destructive(double x[], const size_t N, const double a[], const double b[], double c[])
{
	int n;
	
	c[0] = c[0] / b[0];
	x[0] = x[0] / b[0];
	
	// loop from 1 to N - 1 inclusive 
	for (n = 1; n < N; n++) {
		double m = 1.0f / (b[n] - a[n] * c[n - 1]);
		c[n] = c[n] * m;
		x[n] = (x[n] - a[n] * x[n - 1]) * m;
	}
	
	// loop from N - 2 to 0 inclusive 
	for (n = N - 2; n >= 0; n--)
		x[n] = x[n] - c[n] * x[n + 1];
}

double *** memAllocDouble3(int n,int r,int c)
{
	int padding=10;
	double *a,**p,***pp;
	int rc=r*c;
	int i,j;
	a=(double*) malloc(sizeof(double)*(n*rc+padding));
	if(a==NULL) {mexErrMsgTxt("memAllocDouble: Memory is too huge.\n"); }
	p=(double**) malloc(sizeof(double*)*n*r);
	pp=(double***) malloc(sizeof(double**)*n);
	for(i=0;i<n;i++) 
		for(j=0;j<r;j++) 
			p[i*r+j]=&a[i*rc+j*c];
	for(i=0;i<n;i++) 
		pp[i]=&p[i*r];
	return(pp);
}

void memFreeDouble3(double ***p)
{
	if(p!=NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p=NULL;
	}
}

double** memAllocDouble2(int r,int c)
{
	int padding=10;
	double *a,**p;
	a=(double*) malloc(sizeof(double)*(r*c+padding));
	if(a==NULL) {mexErrMsgTxt("memAllocDouble: Memory is too huge.\n"); }
	p=(double**) malloc(sizeof(double*)*r);
	for(int i=0;i<r;i++) p[i]= &a[i*c];
	return(p);
}
void memFreeDouble2(double **p)
{
	if(p!=NULL)
	{
		free(p[0]);
		free(p);
		p=NULL;
	}
}
