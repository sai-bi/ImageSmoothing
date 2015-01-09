function [E, D] = decomposeFast(rgb, sigma_r, nvec, nsamples)


% Some people advocate LAB...
X = patchToPoints(rgb2lab(rgb)/100);
% X = patchToPoints(rgb);

if(~exist('nsamples', 'var')),
    n = 50;
else
    n = nsamples;
end

subs = round(size(X,1)/n);
find = 1:size(X,1);

ind1 = 1:subs:size(X,1);
ind2 = setdiff(find, ind1);

p = randperm(size(X,1));
ind1 = p(ind1);
ind2 = p(ind2);

Ximp = X(ind1,:)';
Xoth = X(ind2,:)';

A = dist(Ximp,Ximp,sigma_r); 
B = dist(Ximp,Xoth,sigma_r);

% Nystrom
n = size(B,1);
m = size(B,2);
Up = [A;B'];

d1 = sum(Up, 1);
d2 = sum(B,1) + sum(B,2)'*pinv(A)*B;
dhat = sqrt(1./[d1 d2])';

A = A.*(dhat(1:n)*dhat(1:n)');
B = B.*(dhat(1:n)*dhat(n+(1:m))');

Asi = sqrtm(pinv(A));
Q = A+Asi*(B*B')*Asi;

% Apparently this is faster
if (true)
%     [U L] = eig(Q);
    [U,L,T] = svd(Q);
else
    opts.issym = 1;
    [U L]  = eigs(Q,nvec+1,'lm',opts) ;
end

V = Up*(Asi*U*diag(1./sqrt(abs(diag(L)))));

E=[]; pres=[]; D=[];
for i = 2:nvec+1
    res = V(:,i)./V(:,1);
    pres(ind1) = res(1:n);
    pres(ind2) = res(n+1:end);
    E(:,i-1) = real(pres);
end
D = diag(L);
D = D(2:nvec+1);

