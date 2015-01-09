clear; close all;

% filter
I = imread('fishmosaic.png');
J = TreeFilterRGB_Uint8(I, 0.01, 2, 0.05, 3);
[h,w,c] = size(I);

% texture editing
T = double(imread('newtexture.png'));
T = imresize(T, [h,w]);
minval = min(T(:));
maxval = max(T(:));
TRN = double(T-minval).*140./double(maxval-minval)-70; 
M = uint8(double(J)+TRN);

% show
figure;imshow([I,J,M]);

