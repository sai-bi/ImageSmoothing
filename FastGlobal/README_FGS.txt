This software provides the reference implementation of the fast global image smoother described in the paper:
 
    Fast Global Image Smoothing based on Weighted Least Squares
    D. Min, S. Choi, J. Lu, B. Ham, K. Sohn, and M. N. Do,
    IEEE Trans. Image Processing, vol. 23, no. 12, pp. 5638-5653, Dec. 2014.

Please cite the paper above if you use this software. For an up-to-date version, refer to:
https://sites.google.com/site/globalsmoothing/

Only 'Algorithm 1' of the paper was implemented in this software. Note that other functions such as sparse data interpolation, L_gamma norm smoothing with the iterative re-weighted least squares (IRLS) algorithm, and robust smoothing using aggregated data term are not provided.

It supports both a pure filtering (with 'img') and a joint filtering (with 'img' and 'guide_img').
For instance, in 'Demo.m', 
Pure filtering, 'F = FGS(img, 0.1, 30^2, [], 3, 4);'
Joint filtering, 'F = FGS(img, 0.1, 30^2, guide_img, 3, 4);'

The code was implemented using C++ with a MATLAB interface (mex file). Please compile C++ files by running 'compile.m' or using the following command: 'mex mexFGS_simple.cpp' or 'mex mexFGS.cpp'.

'mexFGS_simple.cpp': an intuitive code that works for an input image and a guidance image with an arbitrary number of color channels. The smoothing function of this code is slower than that of 'mexFGS.cpp'.

'mexFGS.cpp': an optimized code used for measuring the runtime of Table I in the paper. The smoothing function works for the following cases only.
'FGS()':	input image with 3-ch and guidance image with 3-ch
'FGS_single()':	input image with 1-ch and guidance image with 1-ch
'FGS_13()': 	input image with 1-ch and guidance image with 3-ch
'FGS_31()':	input image with 3-ch and guidance image with 1-ch
