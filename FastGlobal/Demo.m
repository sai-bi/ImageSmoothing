close all;
clc;

%% Edge-preserving smoothing example

% % sigma = 0.1
% % lambda = 30^2
% % iteration = 3
% % attenuation = 4
img = imread('noisy_13.png');
F = FGS(img, 0.1, 30^2, [], 3, 4);
figure,imshow(F);


% % sigma = 0.03
% % lambda = 20^2
% % iteration = 3 (default)
% % attenuation = 4 (default)
% img = imread('lamp.jpg');
% F = FGS(img, 0.03, 20^2);
% figure,imshow(F);