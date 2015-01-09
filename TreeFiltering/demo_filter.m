clear; close all;

I1 = imread('baboon.png');
tic;
J1 = TreeFilterRGB_Uint8(I1, 0.1, 4);
toc;
figure;imshow([I1,J1]);

I2 = imread('monalisamosaic.jpg');
tic;
J2 = TreeFilterRGB_Uint8(I2, 0.01, 3, 0.08, 4);
toc;
figure;imshow([I2,J2]);
