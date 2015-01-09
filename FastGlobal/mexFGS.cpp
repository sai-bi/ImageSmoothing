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
void FGS(double*** image, double*** joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation);
void FGS_single(double*** image, double*** joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation);
void FGS_13(double*** image, double*** joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation);
void FGS_31(double*** image, double*** joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation);

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
    
    if(nChannels == 1 && nChannels_guide == 1)
        FGS_single(image_filtered, image_guidance, sigma, lambda, solver_iteration, solver_attenuation);	
    else if(nChannels == 3 && nChannels_guide == 3)
        FGS(image_filtered, image_guidance, sigma, lambda, solver_iteration, solver_attenuation);	
    else if(nChannels == 1 && nChannels_guide == 3)
        FGS_13(image_filtered, image_guidance, sigma, lambda, solver_iteration, solver_attenuation);	
    else if(nChannels == 3 && nChannels_guide == 1)
        FGS_31(image_filtered, image_guidance, sigma, lambda, solver_iteration, solver_attenuation);
    else
        mexErrMsgTxt("Unsupported image channels\n");

    mexPrintf("Elapsed time is %f seconds.\n", double(clock()-m_begin)/CLOCKS_PER_SEC);
    
	// output
	image_result = plhs[0] = mxDuplicateArray(img);
	double* ptr_output = mxGetPr(image_result);
	for(int y=0;y<H;y++) for(int x=0;x<W;x++) for(int c=0;c<nChannels;c++) ptr_output[x*H + y + c*(H*W)] = image_filtered[y][x][c];

	memFreeDouble3(image_filtered);
	if(image_guidance) memFreeDouble3(image_guidance);
}

void FGS(double*** image, double*** joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation)
{
	if(joint_image == NULL) joint_image = image;
	
	int width = W;
	int height = H;

	double *a_vec = (double *)malloc(sizeof(double)*width);
	double *b_vec = (double *)malloc(sizeof(double)*width);
	double *c_vec = (double *)malloc(sizeof(double)*width);	
	double *a2_vec = (double *)malloc(sizeof(double)*height);
	double *b2_vec = (double *)malloc(sizeof(double)*height);
	double *c2_vec = (double *)malloc(sizeof(double)*height);

	double lambda_in = 1.5*lambda*pow(4.0,solver_iteration-1)/(pow(4.0,solver_iteration)-1.0);	
	prepareBLFKernel(sigma); // prepare kernel LUT
  
	for(int iter=0;iter<solver_iteration;iter++) {
        // row
		for(int i=0;i<height;i++) {
			memset(a_vec, 0, sizeof(double)*width);
			memset(b_vec, 0, sizeof(double)*width);
			memset(c_vec, 0, sizeof(double)*width);

			double* pc = joint_image[i][0];
			double tpr = *pc++;
			double tpg = *pc++;
			double tpb = *pc++;
			double* pa_vec = &a_vec[1];
			double* pc_vec = &c_vec[0];
			for(int j=1;j<width;j++) { 
				double tcr = *pc++;
				double tcg = *pc++;
				double tcb = *pc++;
				
				// bilateral color weights
  /////////////////////////////////////////////////////////////////////////////////////
				int range = (SQ(tpr-tcr)+SQ(tpg-tcg)+SQ(tpb-tcb));
  /////////////////////////////////////////////////////////////////////////////////////
				double weight = BLFKernelI[range];
				 				
 				*pa_vec = -lambda_in*weight;
 				*pc_vec = *pa_vec; pa_vec++; pc_vec++;
 				
 				tpr = tcr; tpg = tcg; tpb = tcb;
			}

			pa_vec = a_vec;
			pc_vec = c_vec;
			double* pb_vec = b_vec;
			for(int j=0;j<width;j++) {
				*pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
			}

			// solver
			c_vec[0] = c_vec[0] / b_vec[0];
			image[i][0][0] = image[i][0][0] / b_vec[0];
			image[i][0][1] = image[i][0][1] / b_vec[0];
			image[i][0][2] = image[i][0][2] / b_vec[0];
			
			pb_vec = &b_vec[1];
			pa_vec = &a_vec[1];
			pc_vec = &c_vec[0];

			pc = image[i][0];
			tpr = *pc++;
			tpg = *pc++;
			tpb = *pc++;
			double* po = image[i][1];
			for (size_t n = 1; n < width; n++) {
				double m = 1.0 / (*pb_vec - *pa_vec * (*pc_vec++));
      			*pc_vec = *pc_vec * m; 
				
				double tcr = *pc++;
				double tcg = *pc++;
				double tcb = *pc++;
				double ncr = (tcr - *pa_vec*tpr ) * m;
				double ncg = (tcg - *pa_vec*tpg ) * m;
				double ncb = (tcb - *pa_vec*tpb ) * m;

				pb_vec++; pa_vec++;
				tpr = ncr; tpg = ncg; tpb = ncb;
				*po++ = ncr;
				*po++ = ncg;
				*po++ = ncb;
			}

			pc_vec = &c_vec[width-2];
			
			tpr = *--pc;
			tpg = *--pc;
			tpb = *--pc;
			*--po;*--po;*--po;
			
			for (size_t n = width - 1; n-- > 0; ) {
				double tcr = *--pc;
				double tcg = *--pc;
				double tcb = *--pc;

				double ncr = tcr - *pc_vec * tpr;
				double ncg = tcg - *pc_vec * tpg;
				double ncb = tcb - *pc_vec * tpb;
				
				pc_vec--;
				tpr = ncr; tpg = ncg; tpb = ncb;
				*--po = ncr;
				*--po = ncg;
				*--po = ncb;
			}
		}

		// column
		for(int j=0;j<width;j++) {
			memset(a2_vec, 0, sizeof(double)*height);
			memset(b2_vec, 0, sizeof(double)*height);
			memset(c2_vec, 0, sizeof(double)*height);

			double* pc = joint_image[0][j];
			double tpr = *pc++;
			double tpg = *pc++;
			double tpb = *pc++;
			double* pa_vec = &a2_vec[1];
			double* pc_vec = &c2_vec[0];
			for(int i=1;i<height;i++) {
				pc+=(width-1)*3;
				double tcr = *pc++;
				double tcg = *pc++;
				double tcb = *pc++;
				
  /////////////////////////////////////////////////////////////////////////////////////
				int range = (SQ(tpr-tcr)+SQ(tpg-tcg)+SQ(tpb-tcb));
  /////////////////////////////////////////////////////////////////////////////////////
				double weight = BLFKernelI[range];
				
				*pa_vec = -lambda_in*weight;				
                *pc_vec = *pa_vec; pa_vec++; pc_vec++;
                tpr = tcr; tpg = tcg; tpb = tcb;
			}

			pa_vec = a2_vec;
			pc_vec = c2_vec;
			double* pb_vec = b2_vec;
			for(int i=0;i<height;i++) {
				*pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
			}

			// solver
			c2_vec[0] = c2_vec[0] / b2_vec[0];
			image[0][j][0] = image[0][j][0] / b2_vec[0];
			image[0][j][1] = image[0][j][1] / b2_vec[0];
			image[0][j][2] = image[0][j][2] / b2_vec[0];

			pb_vec = &b2_vec[1];
			pa_vec = &a2_vec[1];
			pc_vec = &c2_vec[0];

			pc = image[0][j];
			tpr = *pc++;
			tpg = *pc++;
			tpb = *pc++;
			double* po = image[1][j];
			for (size_t n = 1; n < height; n++) {
				double m = 1.0 / (*pb_vec - *pa_vec * (*pc_vec++));
				*pc_vec = *pc_vec * m; 

				pc+=(width-1)*3;
				double tcr = *pc++;
				double tcg = *pc++;
				double tcb = *pc++;

				double ncr = (tcr - *pa_vec*tpr ) * m;
				double ncg = (tcg - *pa_vec*tpg ) * m;
				double ncb = (tcb - *pa_vec*tpb ) * m;

				pb_vec++; pa_vec++;
				tpr = ncr; tpg = ncg; tpb = ncb;
				
                *po++ = ncr;
				*po++ = ncg;
				*po++ = ncb;
                po+=(width-1)*3;			
			}

            pc_vec = &c2_vec[height-2];

			tpr = *--pc;
			tpg = *--pc;
			tpb = *--pc;
						
			po-=(width-1)*6; po-=3;
            			
			for (size_t n = height - 1; n-- > 0; ) {
				pc-=(width-1)*3;
				double tcr = *--pc;
				double tcg = *--pc;
				double tcb = *--pc;

				double ncr = tcr - *pc_vec * tpr;
				double ncg = tcg - *pc_vec * tpg;
				double ncb = tcb - *pc_vec * tpb;
 				pc_vec--;
				tpr = ncr; tpg = ncg; tpb = ncb;
				*--po = ncr;
				*--po = ncg;
				*--po = ncb;
				po-=(width-1)*3;
			}
		}

		lambda_in /= solver_attenuation; //solver_attenuation = 4 (default parameter)
	}

    free(a_vec); free(b_vec); free(c_vec);
    free(a2_vec); free(b2_vec); free(c2_vec);
}

void FGS_single(double*** image, double*** joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation)
{
	if(joint_image == NULL) joint_image = image;
	
	int width = W;
	int height = H;

	double *a_vec = (double *)malloc(sizeof(double)*width);
	double *b_vec = (double *)malloc(sizeof(double)*width);
	double *c_vec = (double *)malloc(sizeof(double)*width);	
	double *a2_vec = (double *)malloc(sizeof(double)*height);
	double *b2_vec = (double *)malloc(sizeof(double)*height);
	double *c2_vec = (double *)malloc(sizeof(double)*height);

	double lambda_in = 1.5*lambda*pow(4.0,solver_iteration-1)/(pow(4.0,solver_iteration)-1.0);	
	prepareBLFKernel(sigma); // prepare kernel LUT
  
	for(int iter=0;iter<solver_iteration;iter++) {
        // row
		for(int i=0;i<height;i++) {
			memset(a_vec, 0, sizeof(double)*width);
			memset(b_vec, 0, sizeof(double)*width);
			memset(c_vec, 0, sizeof(double)*width);

			double* pc = joint_image[i][0];
			double tpr = *pc++;
			//double tpg = *pc++;
			//double tpb = *pc++;
			double* pa_vec = &a_vec[1];
			double* pc_vec = &c_vec[0];
			for(int j=1;j<width;j++) { 
				double tcr = *pc++;
				//double tcg = *pc++;
				//double tcb = *pc++;
				
				// bilateral color weights
  /////////////////////////////////////////////////////////////////////////////////////
				int range = SQ(tpr-tcr);
  /////////////////////////////////////////////////////////////////////////////////////
				double weight = BLFKernelI[range];
				 				
 				*pa_vec = -lambda_in*weight;
 				*pc_vec = *pa_vec; pa_vec++; pc_vec++;
 				
 				tpr = tcr; //tpg = tcg; tpb = tcb;
			}

			pa_vec = a_vec;
			pc_vec = c_vec;
			double* pb_vec = b_vec;
			for(int j=0;j<width;j++) {
				*pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
			}
            
			// solver
			c_vec[0] = c_vec[0] / b_vec[0];
			image[i][0][0] = image[i][0][0] / b_vec[0];
			//image[i][0][1] = image[i][0][1] / b_vec[0];
			//image[i][0][2] = image[i][0][2] / b_vec[0];
			
			pb_vec = &b_vec[1];
			pa_vec = &a_vec[1];
			pc_vec = &c_vec[0];

			pc = image[i][0];
			tpr = *pc++;
			//tpg = *pc++;
			//tpb = *pc++;
			double* po = image[i][1];
			for (size_t n = 1; n < width; n++) {
				double m = 1.0 / (*pb_vec - *pa_vec * (*pc_vec++));
      			*pc_vec = *pc_vec * m; 
				
				double tcr = *pc++;
				//double tcg = *pc++;
				//double tcb = *pc++;
				double ncr = (tcr - *pa_vec*tpr ) * m;
				//double ncg = (tcg - *pa_vec*tpg ) * m;
				//double ncb = (tcb - *pa_vec*tpb ) * m;

				pb_vec++; pa_vec++;
				tpr = ncr; //tpg = ncg; tpb = ncb;
				*po++ = ncr;
				//*po++ = ncg;
				//*po++ = ncb;
			}

			pc_vec = &c_vec[width-2];
			
			tpr = *--pc;
			//tpg = *--pc;
			//tpb = *--pc;
			*--po;//*--po;*--po;
			
			for (size_t n = width - 1; n-- > 0; ) {
				double tcr = *--pc;
				//double tcg = *--pc;
				//double tcb = *--pc;

				double ncr = tcr - *pc_vec * tpr;
				//double ncg = tcg - *pc_vec * tpg;
				//double ncb = tcb - *pc_vec * tpb;
				
				pc_vec--;
				tpr = ncr; //tpg = ncg; tpb = ncb;
				*--po = ncr;
				//*--po = ncg;
				//*--po = ncb;
			}
		}

		// column
		for(int j=0;j<width;j++) {
			memset(a2_vec, 0, sizeof(double)*height);
			memset(b2_vec, 0, sizeof(double)*height);
			memset(c2_vec, 0, sizeof(double)*height);

			double* pc = joint_image[0][j];
			double tpr = *pc++;
			//double tpg = *pc++;
			//double tpb = *pc++;
			double* pa_vec = &a2_vec[1];
			double* pc_vec = &c2_vec[0];
			for(int i=1;i<height;i++) {
				pc+=(width-1); // pc+=(width-1)*3;
				double tcr = *pc++;
				//double tcg = *pc++;
				//double tcb = *pc++;
				
  /////////////////////////////////////////////////////////////////////////////////////
				int range = SQ(tpr-tcr);
  /////////////////////////////////////////////////////////////////////////////////////
				double weight = BLFKernelI[range];
				
				*pa_vec = -lambda_in*weight;				
                *pc_vec = *pa_vec; pa_vec++; pc_vec++;
                tpr = tcr; //tpg = tcg; tpb = tcb;
			}

			pa_vec = a2_vec;
			pc_vec = c2_vec;
			double* pb_vec = b2_vec;
			for(int i=0;i<height;i++) {
				*pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
			}

			// solver
			c2_vec[0] = c2_vec[0] / b2_vec[0];
			image[0][j][0] = image[0][j][0] / b2_vec[0];
			//image[0][j][1] = image[0][j][1] / b2_vec[0];
			//image[0][j][2] = image[0][j][2] / b2_vec[0];

			pb_vec = &b2_vec[1];
			pa_vec = &a2_vec[1];
			pc_vec = &c2_vec[0];

			pc = image[0][j];
			tpr = *pc++;
			//tpg = *pc++;
			//tpb = *pc++;
			double* po = image[1][j];
			for (size_t n = 1; n < height; n++) {
				double m = 1.0 / (*pb_vec - *pa_vec * (*pc_vec++));
				*pc_vec = *pc_vec * m; 

				pc+=(width-1); //pc+=(width-1)*3;
				double tcr = *pc++;
				//double tcg = *pc++;
				//double tcb = *pc++;

				double ncr = (tcr - *pa_vec*tpr ) * m;
				//double ncg = (tcg - *pa_vec*tpg ) * m;
				//double ncb = (tcb - *pa_vec*tpb ) * m;

				pb_vec++; pa_vec++;
				tpr = ncr; //tpg = ncg; tpb = ncb;
				
                *po++ = ncr;
				//*po++ = ncg;
				//*po++ = ncb;
                po+=(width-1); //po+=(width-1)*3;
			}

            pc_vec = &c2_vec[height-2];

			tpr = *--pc;
			//tpg = *--pc;
			//tpb = *--pc;
						
			po-=(width-1)*2; po-=1; //po-=(width-1)*6; po-=3;
            			
			for (size_t n = height - 1; n-- > 0; ) {
				pc-=(width-1); //pc-=(width-1)*3;
				double tcr = *--pc;
				//double tcg = *--pc;
				//double tcb = *--pc;

				double ncr = tcr - *pc_vec * tpr;
				//double ncg = tcg - *pc_vec * tpg;
				//double ncb = tcb - *pc_vec * tpb;
 				pc_vec--;
				tpr = ncr; //tpg = ncg; tpb = ncb;
				*--po = ncr;
				//*--po = ncg;
				//*--po = ncb;
				po-=(width-1); //po-=(width-1)*3;
			}
		}

		lambda_in /= solver_attenuation; //solver_attenuation = 4 (default parameter)
	}

    free(a_vec); free(b_vec); free(c_vec);
    free(a2_vec); free(b2_vec); free(c2_vec);
}

void FGS_13(double*** image, double*** joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation)
{
	int width = W;
	int height = H;

	double *a_vec = (double *)malloc(sizeof(double)*width);
	double *b_vec = (double *)malloc(sizeof(double)*width);
	double *c_vec = (double *)malloc(sizeof(double)*width);	
	double *a2_vec = (double *)malloc(sizeof(double)*height);
	double *b2_vec = (double *)malloc(sizeof(double)*height);
	double *c2_vec = (double *)malloc(sizeof(double)*height);

	double lambda_in = 1.5*lambda*pow(4.0,solver_iteration-1)/(pow(4.0,solver_iteration)-1.0);	
	prepareBLFKernel(sigma); // prepare kernel LUT
  
	for(int iter=0;iter<solver_iteration;iter++) {
        // row
		for(int i=0;i<height;i++) {
			memset(a_vec, 0, sizeof(double)*width);
			memset(b_vec, 0, sizeof(double)*width);
			memset(c_vec, 0, sizeof(double)*width);

			double* pc = joint_image[i][0];
			double tpr = *pc++;
			double tpg = *pc++;
			double tpb = *pc++;
			double* pa_vec = &a_vec[1];
			double* pc_vec = &c_vec[0];
			for(int j=1;j<width;j++) { 
				double tcr = *pc++;
				double tcg = *pc++;
				double tcb = *pc++;
				
				// bilateral color weights
  /////////////////////////////////////////////////////////////////////////////////////
				int range = (SQ(tpr-tcr)+SQ(tpg-tcg)+SQ(tpb-tcb));
  /////////////////////////////////////////////////////////////////////////////////////
				double weight = BLFKernelI[range];
				 				
 				*pa_vec = -lambda_in*weight;
 				*pc_vec = *pa_vec; pa_vec++; pc_vec++;
 				
 				tpr = tcr; tpg = tcg; tpb = tcb;
			}

			pa_vec = a_vec;
			pc_vec = c_vec;
			double* pb_vec = b_vec;
			for(int j=0;j<width;j++) {
				*pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
			}

			// solver
			c_vec[0] = c_vec[0] / b_vec[0];
			image[i][0][0] = image[i][0][0] / b_vec[0];
			//image[i][0][1] = image[i][0][1] / b_vec[0];
			//image[i][0][2] = image[i][0][2] / b_vec[0];
			
			pb_vec = &b_vec[1];
			pa_vec = &a_vec[1];
			pc_vec = &c_vec[0];

			pc = image[i][0];
			tpr = *pc++;
			//tpg = *pc++;
			//tpb = *pc++;
			double* po = image[i][1];
			for (size_t n = 1; n < width; n++) {
				double m = 1.0 / (*pb_vec - *pa_vec * (*pc_vec++));
      			*pc_vec = *pc_vec * m; 
				
				double tcr = *pc++;
				//double tcg = *pc++;
				//double tcb = *pc++;
				double ncr = (tcr - *pa_vec*tpr ) * m;
				//double ncg = (tcg - *pa_vec*tpg ) * m;
				//double ncb = (tcb - *pa_vec*tpb ) * m;

				pb_vec++; pa_vec++;
				tpr = ncr; //tpg = ncg; tpb = ncb;
				*po++ = ncr;
				//*po++ = ncg;
				//*po++ = ncb;
			}

			pc_vec = &c_vec[width-2];
			
			tpr = *--pc;
			//tpg = *--pc;
			//tpb = *--pc;
			*--po;//*--po;*--po;
			
			for (size_t n = width - 1; n-- > 0; ) {
				double tcr = *--pc;
				//double tcg = *--pc;
				//double tcb = *--pc;

				double ncr = tcr - *pc_vec * tpr;
				//double ncg = tcg - *pc_vec * tpg;
				//double ncb = tcb - *pc_vec * tpb;
				
				pc_vec--;
				tpr = ncr; //tpg = ncg; tpb = ncb;
				*--po = ncr;
				//*--po = ncg;
				//*--po = ncb;
			}
		}

		// column
		for(int j=0;j<width;j++) {
			memset(a2_vec, 0, sizeof(double)*height);
			memset(b2_vec, 0, sizeof(double)*height);
			memset(c2_vec, 0, sizeof(double)*height);

			double* pc = joint_image[0][j];
			double tpr = *pc++;
			double tpg = *pc++;
			double tpb = *pc++;
			double* pa_vec = &a2_vec[1];
			double* pc_vec = &c2_vec[0];
			for(int i=1;i<height;i++) {
				pc+=(width-1)*3;
				double tcr = *pc++;
				double tcg = *pc++;
				double tcb = *pc++;
				
  /////////////////////////////////////////////////////////////////////////////////////
				int range = (SQ(tpr-tcr)+SQ(tpg-tcg)+SQ(tpb-tcb));
  /////////////////////////////////////////////////////////////////////////////////////
				double weight = BLFKernelI[range];
				
				*pa_vec = -lambda_in*weight;				
                *pc_vec = *pa_vec; pa_vec++; pc_vec++;
                tpr = tcr; tpg = tcg; tpb = tcb;
			}

			pa_vec = a2_vec;
			pc_vec = c2_vec;
			double* pb_vec = b2_vec;
			for(int i=0;i<height;i++) {
				*pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
			}

			// solver
			c2_vec[0] = c2_vec[0] / b2_vec[0];
			image[0][j][0] = image[0][j][0] / b2_vec[0];
			//image[0][j][1] = image[0][j][1] / b2_vec[0];
			//image[0][j][2] = image[0][j][2] / b2_vec[0];

			pb_vec = &b2_vec[1];
			pa_vec = &a2_vec[1];
			pc_vec = &c2_vec[0];

			pc = image[0][j];
			tpr = *pc++;
			//tpg = *pc++;
			//tpb = *pc++;
			double* po = image[1][j];
			for (size_t n = 1; n < height; n++) {
				double m = 1.0 / (*pb_vec - *pa_vec * (*pc_vec++));
				*pc_vec = *pc_vec * m; 

				pc+=(width-1); //pc+=(width-1)*3;
				double tcr = *pc++;
				//double tcg = *pc++;
				//double tcb = *pc++;

				double ncr = (tcr - *pa_vec*tpr ) * m;
				//double ncg = (tcg - *pa_vec*tpg ) * m;
				//double ncb = (tcb - *pa_vec*tpb ) * m;

				pb_vec++; pa_vec++;
				tpr = ncr; //tpg = ncg; tpb = ncb;
				
                *po++ = ncr;
				//*po++ = ncg;
				//*po++ = ncb;
                po+=(width-1); //po+=(width-1)*3;
			}

            pc_vec = &c2_vec[height-2];

			tpr = *--pc;
			//tpg = *--pc;
			//tpb = *--pc;
						
			po-=(width-1)*2; po-=1; //po-=(width-1)*6; po-=3;
            			
			for (size_t n = height - 1; n-- > 0; ) {
				pc-=(width-1); //pc-=(width-1)*3;
				double tcr = *--pc;
				//double tcg = *--pc;
				//double tcb = *--pc;

				double ncr = tcr - *pc_vec * tpr;
				//double ncg = tcg - *pc_vec * tpg;
				//double ncb = tcb - *pc_vec * tpb;
 				pc_vec--;
				tpr = ncr; //tpg = ncg; tpb = ncb;
				*--po = ncr;
				//*--po = ncg;
				//*--po = ncb;
				po-=(width-1); //po-=(width-1)*3;
			}
		}

		lambda_in /= solver_attenuation; //solver_attenuation = 4 (default parameter)
	}

    free(a_vec); free(b_vec); free(c_vec);
    free(a2_vec); free(b2_vec); free(c2_vec);
}


void FGS_31(double*** image, double*** joint_image, double sigma, double lambda, int solver_iteration, double solver_attenuation)
{
	int width = W;
	int height = H;

	double *a_vec = (double *)malloc(sizeof(double)*width);
	double *b_vec = (double *)malloc(sizeof(double)*width);
	double *c_vec = (double *)malloc(sizeof(double)*width);	
	double *a2_vec = (double *)malloc(sizeof(double)*height);
	double *b2_vec = (double *)malloc(sizeof(double)*height);
	double *c2_vec = (double *)malloc(sizeof(double)*height);

	double lambda_in = 1.5*lambda*pow(4.0,solver_iteration-1)/(pow(4.0,solver_iteration)-1.0);	
	prepareBLFKernel(sigma); // prepare kernel LUT
  
	for(int iter=0;iter<solver_iteration;iter++) {
        // row
		for(int i=0;i<height;i++) {
			memset(a_vec, 0, sizeof(double)*width);
			memset(b_vec, 0, sizeof(double)*width);
			memset(c_vec, 0, sizeof(double)*width);

			double* pc = joint_image[i][0];
			double tpr = *pc++;
			double tpg;// = *pc++;
			double tpb;// = *pc++;
			double* pa_vec = &a_vec[1];
			double* pc_vec = &c_vec[0];
			for(int j=1;j<width;j++) { 
				double tcr = *pc++;
				//double tcg = *pc++;
				//double tcb = *pc++;
				
				// bilateral color weights
  /////////////////////////////////////////////////////////////////////////////////////
				int range = SQ(tpr-tcr);
  /////////////////////////////////////////////////////////////////////////////////////
				double weight = BLFKernelI[range];
				 				
 				*pa_vec = -lambda_in*weight;
 				*pc_vec = *pa_vec; pa_vec++; pc_vec++;
 				
 				tpr = tcr; //tpg = tcg; tpb = tcb;
			}

			pa_vec = a_vec;
			pc_vec = c_vec;
			double* pb_vec = b_vec;
			for(int j=0;j<width;j++) {
				*pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
			}

			// solver
			c_vec[0] = c_vec[0] / b_vec[0];
			image[i][0][0] = image[i][0][0] / b_vec[0];
			image[i][0][1] = image[i][0][1] / b_vec[0];
			image[i][0][2] = image[i][0][2] / b_vec[0];
			
			pb_vec = &b_vec[1];
			pa_vec = &a_vec[1];
			pc_vec = &c_vec[0];

			pc = image[i][0];
			tpr = *pc++;
			tpg = *pc++;
			tpb = *pc++;
			double* po = image[i][1];
			for (size_t n = 1; n < width; n++) {
				double m = 1.0 / (*pb_vec - *pa_vec * (*pc_vec++));
      			*pc_vec = *pc_vec * m; 
				
				double tcr = *pc++;
				double tcg = *pc++;
				double tcb = *pc++;
				double ncr = (tcr - *pa_vec*tpr ) * m;
				double ncg = (tcg - *pa_vec*tpg ) * m;
				double ncb = (tcb - *pa_vec*tpb ) * m;

				pb_vec++; pa_vec++;
				tpr = ncr; tpg = ncg; tpb = ncb;
				*po++ = ncr;
				*po++ = ncg;
				*po++ = ncb;
			}

			pc_vec = &c_vec[width-2];
			
			tpr = *--pc;
			tpg = *--pc;
			tpb = *--pc;
			*--po;*--po;*--po;
			
			for (size_t n = width - 1; n-- > 0; ) {
				double tcr = *--pc;
				double tcg = *--pc;
				double tcb = *--pc;

				double ncr = tcr - *pc_vec * tpr;
				double ncg = tcg - *pc_vec * tpg;
				double ncb = tcb - *pc_vec * tpb;
				
				pc_vec--;
				tpr = ncr; tpg = ncg; tpb = ncb;
				*--po = ncr;
				*--po = ncg;
				*--po = ncb;
			}
		}

		// column
		for(int j=0;j<width;j++) {
			memset(a2_vec, 0, sizeof(double)*height);
			memset(b2_vec, 0, sizeof(double)*height);
			memset(c2_vec, 0, sizeof(double)*height);

			double* pc = joint_image[0][j];
			double tpr = *pc++;
			double tpg;// = *pc++;
			double tpb;// = *pc++;
			double* pa_vec = &a2_vec[1];
			double* pc_vec = &c2_vec[0];
			for(int i=1;i<height;i++) {
				pc+=(width-1); //pc+=(width-1)*3;
				double tcr = *pc++;
				//double tcg = *pc++;
				//double tcb = *pc++;
				
  /////////////////////////////////////////////////////////////////////////////////////
				int range = SQ(tpr-tcr);
  /////////////////////////////////////////////////////////////////////////////////////
				double weight = BLFKernelI[range];
				
				*pa_vec = -lambda_in*weight;				
                *pc_vec = *pa_vec; pa_vec++; pc_vec++;
                tpr = tcr; //tpg = tcg; tpb = tcb;
			}

			pa_vec = a2_vec;
			pc_vec = c2_vec;
			double* pb_vec = b2_vec;
			for(int i=0;i<height;i++) {
				*pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
			}

			// solver
			c2_vec[0] = c2_vec[0] / b2_vec[0];
			image[0][j][0] = image[0][j][0] / b2_vec[0];
			image[0][j][1] = image[0][j][1] / b2_vec[0];
			image[0][j][2] = image[0][j][2] / b2_vec[0];

			pb_vec = &b2_vec[1];
			pa_vec = &a2_vec[1];
			pc_vec = &c2_vec[0];

			pc = image[0][j];
			tpr = *pc++;
			tpg = *pc++;
			tpb = *pc++;
			double* po = image[1][j];
			for (size_t n = 1; n < height; n++) {
				double m = 1.0 / (*pb_vec - *pa_vec * (*pc_vec++));
				*pc_vec = *pc_vec * m; 

				pc+=(width-1)*3;
				double tcr = *pc++;
				double tcg = *pc++;
				double tcb = *pc++;

				double ncr = (tcr - *pa_vec*tpr ) * m;
				double ncg = (tcg - *pa_vec*tpg ) * m;
				double ncb = (tcb - *pa_vec*tpb ) * m;

				pb_vec++; pa_vec++;
				tpr = ncr; tpg = ncg; tpb = ncb;
				
                *po++ = ncr;
				*po++ = ncg;
				*po++ = ncb;
                po+=(width-1)*3;			
			}

            pc_vec = &c2_vec[height-2];

			tpr = *--pc;
			tpg = *--pc;
			tpb = *--pc;
						
			po-=(width-1)*6; po-=3;
            			
			for (size_t n = height - 1; n-- > 0; ) {
				pc-=(width-1)*3;
				double tcr = *--pc;
				double tcg = *--pc;
				double tcb = *--pc;

				double ncr = tcr - *pc_vec * tpr;
				double ncg = tcg - *pc_vec * tpg;
				double ncb = tcb - *pc_vec * tpb;
 				pc_vec--;
				tpr = ncr; tpg = ncg; tpb = ncb;
				*--po = ncr;
				*--po = ncg;
				*--po = ncb;
				po-=(width-1)*3;
			}
		}

		lambda_in /= solver_attenuation; //solver_attenuation = 4 (default parameter)
	}

    free(a_vec); free(b_vec); free(c_vec);
    free(a2_vec); free(b2_vec); free(c2_vec);
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
