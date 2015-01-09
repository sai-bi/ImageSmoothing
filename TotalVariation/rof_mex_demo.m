%              rof_mex_demo.m  by Tom Goldstein
% This code tests the method defined by "splitBregmanROF.c".  This code
% must be compiled using the "mex" command before this demo will work.

 
% Step 1:  Get the test image
exact = double(imread('cameraman.tif'));
dims = size(exact);

noisy = exact+15*randn(dims);

% Step 2: Denoise the Image

clean = SplitBregmanROF(noisy,.05,0.001);

% Step 3:  Display Results

close all;
figure;
subplot(2,2,1);
imagesc(exact);
colormap(gray);
title('Original');

subplot(2,2,2);
imagesc(noisy);
colormap(gray);
title('noisy');

subplot(2,2,3);
imagesc(clean);
colormap(gray);
title('denoised');

subplot(2,2,4);
imagesc(noisy-clean);
colormap(gray);
title('difference');
