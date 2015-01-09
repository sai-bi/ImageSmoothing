function [ out ] = ag( I, gamma )
%AG Raise I in power gamma (ApplyGamma)

if nargin==1
    gamma = 1/2.2;
end

out = real(I.^gamma);