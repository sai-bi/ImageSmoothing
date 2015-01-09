/////////////////////////////////////////////////////////////////////
//
//
//  The C-mex implementation of euclidean pairwise metrics
//
//
/////////////////////////////////////////////////////////////////////

#include "mex.h"

#include <math.h>


template<typename Ty>
void euclidean(const Ty* X1, const Ty* X2, Ty* M, int d, int n1, int n2)
{
    const Ty* v2 = X2;
    for(int j = 0; j < n2; ++j)
    {        
        const Ty* v1 = X1;
        for (int i = 0; i < n1; ++i)
        {
            Ty s = 0;
            for (int k = 0; k < d; ++k) 
            {
                s += (v1[k] - v2[k])*(v1[k] - v2[k]);
            }
            *M++ = s;           
            
            v1 += d;
        }        
        v2 += d;
    }
}

template<typename Ty>
void euclidean_exp(const Ty* X1, const Ty* X2, Ty* M, int d, int n1, int n2, double sigma_r)
{
	double sigma = 1.0/sigma_r;
    const Ty* v2 = X2;
    for(int j = 0; j < n2; ++j)
    {        
        const Ty* v1 = X1;
        for (int i = 0; i < n1; ++i)
        {
            Ty s = 0;
            for (int k = 0; k < d; ++k) 
            {
                s += (v1[k] - v2[k])*(v1[k] - v2[k]);
            }
            *M++ = exp(-sigma*s);           
            
            v1 += d;
        }        
        v2 += d;
    }
}



//void euclidean(const double* X1, const double* X2, double* M, int d, int n1, int n2)
//{
//    const double* v2 = X2;
//    for(int j = 0; j < n2; ++j)
//    {        
//        const double* v1 = X1;
//        for (int i = 0; i < n1; ++i)
//        {
//		  double s = 0;
//		  for (int k = 0; k < d; ++k) 
//		  {
//			  s += (v1[k] - v2[k])*(v1[k] - v2[k]);
//		  }
//            *M++ = s;           
//            
//            v1 += d;
//        }        
//        v2 += d;
//    }
//}

/**
 * main entry
 * Input:   X1, X2, sigma(opcode)
 * Output:  D
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray* mxX1 = prhs[0];
    const mxArray* mxX2 = prhs[1];
    const mxArray* mxSigma = prhs[2];
    
    mwSize d = mxGetM(mxX1);
    mwSize n1 = mxGetN(mxX1);
    mwSize n2 = mxGetN(mxX2);
    
    double sigma_r = *((const double*)mxGetData(mxSigma));
            
    if (mxGetClassID(mxX1) == mxDOUBLE_CLASS)
    {
        mxArray* mxD = mxCreateNumericMatrix(n1, n2, mxDOUBLE_CLASS, mxREAL);
		euclidean_exp((const double*)mxGetData(mxX1), (const double*)mxGetData(mxX2), 
				  (double *)mxGetData(mxD), (int)d, (int)n1, (int)n2,sigma_r);
        plhs[0] = mxD;
    }
    else // SINGLE
    {
        mxArray* mxD = mxCreateNumericMatrix(n1, n2, mxSINGLE_CLASS, mxREAL);
        euclidean_exp((const float*)mxGetData(mxX1), (const float*)mxGetData(mxX2), 
                        (float*)mxGetData(mxD), (int)d, (int)n1, (int)n2,sigma_r);
        plhs[0] = mxD;
    }        
}




