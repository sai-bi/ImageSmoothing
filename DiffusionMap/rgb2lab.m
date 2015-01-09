function lab = rgb2lab( rgb )
%RGB2LAB Summary of this function goes here
%   Detailed explanation goes here

C = makecform('srgb2lab');
lab = applycform(rgb,C);
