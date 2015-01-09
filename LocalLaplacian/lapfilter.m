% Laplacian Filtering 
%   - public Matlab implementation for reproducibility
%   - about 30x slower than our single-thread C++ version
%
% This script implements edge-aware detail and tone manipulation as 
% described in Paris, Hasinoff, and Kautz, "Local Laplacian Filters:
% Edge-aware Image Processing with a Laplacian Pyramid", ACM 
% Transactions on Graphics (Proc. SIGGRAPH 2011), 30(4), 2011.
%
% This is a wrapper around the core algorithm (see lapfilter_core.m).
% It defines the remapping function, the color treatment, the processing 
% domain (linear or log), and it implements simple postprocessing
% for our tone mapping results.
%
% Arguments:
%   image 'I'
%   pixel-wise remapping parameters 'sigma_r', 'alpha', 'beta'
%   remapping method for color images 'colorRemapping' ['rgb' or 'lum']
%   processing domain 'domain' ['lin' or 'log']
%
% sam.hasinoff@gmail.com, April 2011
%

function R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain)

% interpret the input arguments
if size(I,3)==1
    I = repmat(I,[1 1 3]); 
end
if strcmp(domain,'log')
    sigma_r = log(sigma_r);
end

% detail remapping function
noise_level = 0.01;
function out = fd(d)
    out = d.^alpha;
    if alpha<1
        tau = smooth_step(noise_level,2*noise_level,d*sigma_r);
        out = tau.*out + (1-tau).*d;
    end
end

% edge remapping function
function out = fe(a)
    out = beta*a;
end

% define the overall pixel-wise remapping function r, using
% the threshold sigma_r for edge-detail separation
switch colorRemapping
    case 'rgb'
        % process pixels as vectors in RGB color space
        r = @(i,g0)(r_color(i,g0,sigma_r,@fd,@fe));
    case 'lum'
        % save RGB color ratios for later, process the luminance
        IY = luminance(I);
        Iratio = I ./ repmat(IY+eps,[1 1 3]);
        I = IY;
        r = @(i,g0)(r_gray(i,g0,sigma_r,@fd,@fe));
    otherwise
        error('invalid color remapping');
end

% define the processing domain
switch domain
    case 'lin',
        to_domain   = @(I) I;
        from_domain = @(R) R;
    case 'log', 
        to_domain   = @(I) log(I + eps);
        from_domain = @(R) exp(R) - eps;
    otherwise
        error('invalid domain');
end

% call the core Laplacian filtering algorithm
if alpha==1 && beta==1    
    R = I;
else
    I = to_domain(I);
    R = lapfilter_core(I,r);   
    R = from_domain(R);
end

% postprocessing
if strcmp(domain,'log') && beta<=1
    % for tone mapping, remap middle 99% of intensities to
    % fixed dynamic range using a gamma curve
    DR_desired = 100;
    prc_clip = 0.5;
    RY = luminance(R);
    Rmax_clip = prctile(RY(:),100-prc_clip);
    Rmin_clip = prctile(RY(:),prc_clip);
    DR_clip = Rmax_clip/Rmin_clip;
    exponent = log(DR_desired)/log(DR_clip);
    R = max(0,R/Rmax_clip) .^ exponent;
end
if strcmp(colorRemapping,'lum')
    % if working with luminance, reintroduce color ratios 
    R = repmat(R,[1 1 3]) .* Iratio;
end
% clip out of bounds intensities
R = max(0,R);
if beta<=1
    R = min(1,R);  
end
if strcmp(domain,'log') && beta<=1
    % for tone mapping, gamma correct linear intensities for display
    gamma_val = 2.2;
    R = R.^(1/gamma_val);
end



%% helper functions

% smooth step edge between (xmin,0) and (xmax,1)
function y = smooth_step(xmin,xmax,x)
    y = (x - xmin)/(xmax - xmin);
    y = max(0,min(1,y));
    y = y.^2.*(y-2).^2;
end

% convert RGB to grayscale intensity
function Y = luminance(I)
    switch size(I,3),
        case 1, Y = I;
        case 3, Y = (20*I(:,:,1) + 40*I(:,:,2) + I(:,:,3))/61;
    end
end

% color remapping function
function inew = r_color(i,g0,sigma_r,fd,fe)
    g0 = repmat(g0,[size(i,1) size(i,2) 1]);
    dnrm = sqrt(sum((i-g0).^2,3));
    unit = (i-g0)./repmat(eps + dnrm,[1 1 3]);
    % detail and edge processing
    rd = g0 + unit.*repmat(sigma_r*fd(dnrm/sigma_r),[1 1 3]);
    re = g0 + unit.*repmat((fe(dnrm - sigma_r) + sigma_r),[1 1 3]);
    % edge-detail separation based on sigma_r threshold
    isedge = repmat(dnrm > sigma_r,[1 1 3]);
    inew = ~isedge.*rd + isedge.*re;
end

% grayscale remapping function
function inew = r_gray(i,g0,sigma_r,fd,fe)
    dnrm = abs(i-g0);
    dsgn = sign(i-g0);
    % detail and edge processing
    rd = g0 + dsgn*sigma_r.*fd(dnrm/sigma_r);
    re = g0 + dsgn.*(fe(dnrm - sigma_r) + sigma_r);
    % edge-detail separation based on sigma_r threshold
    isedge = dnrm > sigma_r;
    inew = ~isedge.*rd + isedge.*re;
end

end
