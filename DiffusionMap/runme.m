%% Compute diffusion maps

% Compute diffusion maps and use them as per-pixel affinity.

% Make sure you compiled dist.cpp for your platform:
% mex dist.cpp

clear 

rgb = im2double(imread('poppy.jpg'));
rgb = imresize(rgb, 0.3);
% Optional: Remove "gamma".
rgb = ag(rgb, 2.2); 
[h,w,~] = size(rgb);

% Try other sigma values (0.05)
GAUSS_SIGMA = 0.1;  
NUM_VECTORS = 7;
NUM_SAMPLES = 40;
[E, D] = decomposeFast(rgb, GAUSS_SIGMA, NUM_VECTORS, NUM_SAMPLES);

% Plot eigenvectors
figure(1)
subplot(2,4,1); imshow(ag(rgb)); title('Original Image');
for i=2:min(8,NUM_VECTORS+1)
    subplot(2,4,i), imagesc(real(reshape(E(:,i-1),h,w))), title(D(i-1))
    title(['EValue' num2str(i) ': ' num2str(D(i-1),4)]);
    axis off
end

%% Fading gradients with different t values

% Let's stack first three eigenvectors as an image and visaulize the
% magnitude of its gradients. This is how we got the WLS with Diffusion
% Maps set of results (Section 4).

orig_map=reshape(E,h,w,NUM_VECTORS);
map = zeros([h w 3]);
D = D./D(1);

t = [1, 2, 4, 16];

alpha = 1.0;
for tt = 1:length(t)
    DD = ag(D,t(tt));
    for i=1:NUM_VECTORS
        map(:,:,i) = orig_map(:,:,i)*ag(DD(i),t(tt));
    end 
    dy = diff(map, 1, 1);
    dy = sum(abs(dy),3).^alpha;
    dy = padarray(dy, [1 0], 'post');
    dx = diff(map, 1, 2); 
    dx = sum(abs(dx),3).^alpha;
    dx = padarray(dx, [0 1], 'post');
    grad = dx+dy;
    geoGrad(:,:,tt) = grad./max2(grad);
end

figure(2);
subplot(2,2,1), imshow(geoGrad(:,:,1)), title(['DM Gradient, t=',num2str(t(1))])
subplot(2,2,2), imshow(geoGrad(:,:,2)), title(['DM Gradient, t=',num2str(t(2))])
subplot(2,2,3), imshow(geoGrad(:,:,3)), title(['DM Gradient, t=',num2str(t(3))])
subplot(2,2,4), imshow(geoGrad(:,:,4)), title(['DM Gradient, t=',num2str(t(4))])


