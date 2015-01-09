function [I, r,g,b] = readLDR(filename, path, maxDim, method)


if nargin==1
   method = 'bicubic';
   path = '~/Images/testImages/';
end

if nargin==3
   method = 'bicubic';
end

fullpath = strcat(path,filename);

rc = 0.2989;
gc = 0.587;
bc = 0.114;

degamma = 2.2;

Im = imread(fullpath);
maxIm = max(Im(:));

if nargin == 1 || nargin == 2
    ;
else
    imMaxDim = max(size(Im));
    if imMaxDim > maxDim
        Im = imresize(Im, maxDim/imMaxDim, method);
    end
end
% Heuristics
m = 256.0;
if maxIm>m
    m = 65535.0;
end
Im = (double(Im)/m).^degamma;

r = 0; g = 0; b = 0;

if ndims(Im)==3
    r = Im(:,:,1); g = Im(:,:,2); b = Im(:,:,3);
    I = r*rc + g*gc + b*bc +0.0000001;
else
    I = Im;
end

    
