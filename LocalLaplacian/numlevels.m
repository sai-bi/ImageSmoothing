function nlev = numlevels(im_sz)
% number of pyramid levels for an image of size [r c]

% as many pyramid levels as possible, up to 1x1
min_d = min(im_sz(1:2));
nlev = 1;
while min_d>1
    nlev = nlev+1;
    min_d = floor((min_d+1)/2);
end

%nlev = floor(log(min(im_sz(1:2))) / log(2));
%nlev = ceil(log(min(im_sz(1:2))) / log(2));

end