% function OUT = TreeFilterRGB_Uint8(uint8_rgbimg,sigma,sig_s[,sig_r=0.05[,num_iter=1]])
% 
% Parameters: 
%       uint8_rgbimg:           the input image, should be 3-channel uint8 image
%       sigma:                  tree distance parameter, typical value 0.1
%       sigma_s:                spatial parameter, typical value 4
%       sigma_r (optional):     range parameter, default value 0.05
%       num_iter (optional):	number of iterations, default value 1
%
% Reference: 
%       Linchao Bao, Qingxiong Yang, Hao Yuan, 
%       "Tree Filtering: Efficient Structure-Preserving Smoothing with a
%       Minimum Spanning Tree", 
%       TIP 2014. 
% 
% Note: The demo is only for academic use. Note that the timing reported
%       in the paper is for single channel image. 
% 
% 