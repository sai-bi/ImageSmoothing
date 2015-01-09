%% example: detail manipulation

fn_in = 'input_png/flower.png';
I = double( imread(fn_in) )/255;
I = min(1,max(0, imresize(I,1/4) ));  % downscale, Matlab version is slow

figure; imshow(I); set(gcf,'name',fn_in)

sigma_r = 0.4;
alpha = 0.25;
beta = 1;
colorRemapping = 'rgb';
domain = 'lin';
R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain);
figure; clf; imshow(R);

%% example: tone manipulation

fn_in = 'input_hdr/doll.hdr';
I = double( hdrread(fn_in) );
I = I(245:500,50:272,:);  % crop, Matlab version is slow

Igamma = lapfilter(I,0,1,1,'lum','log');  % simple gamma tonemapping
figure; imshow(Igamma); set(gcf,'name',fn_in)

sigma_r = 2.5;
alpha = 1;
beta = 0;
colorRemapping = 'lum';
domain = 'log';
R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain);
figure; clf; imshow(R);

%% generate supplementary material: detail manipulation

beta = 1;
colorRemapping = 'rgb';
domain = 'lin';

srcdir = 'input_png/';
fnlist = dir([srcdir,'*.png']);
mkdir('results_detail')
for F = 1:length(fnlist)    
    fn_in = [srcdir,fnlist(F).name]
    I = double( imread(fn_in) )/255;
    %figure; imshow(I); set(gcf,'name',fn_in)
    
    for sigma_r = [0.1 0.2 0.4]
        for alpha = [0.25 0.5 2 4]
            R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain);
            fn_out = sprintf('results_detail/%s_%s_%s_s%g_a%g_b%g.png',fnlist(F).name(1:end-4),colorRemapping,domain,sigma_r,alpha,beta)
            imwrite(R,fn_out)
            %figure; clf; imshow(R); set(gcf,'name',fn_out)
        end
    end    
end

%% generate supplementary material: tone manipulation

sigma_r = 2.5;
colorRemapping = 'lum';
domain = 'log';

srcdir = 'input_hdr/';
fnlist = dir([srcdir,'*.hdr']);
mkdir('results_hdr')
for F = 1:length(fnlist)
    fn_in = [srcdir,fnlist(F).name]
    I = double( hdrread(fn_in) );
    %figure; imshow(I); set(gcf,'name',fn_in)
    
    for alpha = [0.25 0.5 0.75 1]
        for beta = [0 0.3 0.6]
            R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain);
            fn_out = sprintf('results_hdr/%s_%s_%s_s%g_a%g_b%g.png',fnlist(F).name(1:end-4),colorRemapping,domain,sigma_r,alpha,beta)
            imwrite(R,fn_out)
            %figure; clf; imshow(R); set(gcf,'name',fn_out)
        end
    end

    alpha = 1;
    beta = 1;
    R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain);
    fn_out = sprintf('results_hdr/%s.png',fnlist(F).name(1:end-4))
    imwrite(R,fn_out)
end

srcdir = 'input_hdr_large/';
fnlist = dir([srcdir,'*.hdr']);
mkdir('results_hdr_large')
for F = 1:length(fnlist)
    fn_in = [srcdir,fnlist(F).name]
    I = double( hdrread(fn_in) );
    %figure; imshow(I); set(gcf,'name',fn_in)
    
    for alpha = [0.25 1]
        for beta = [0]
            R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain);
            fn_out = sprintf('results_hdr_large/%s_%s_%s_s%g_a%g_b%g.png',fnlist(F).name(1:end-4),colorRemapping,domain,sigma_r,alpha,beta)
            imwrite(R,fn_out)
            %figure; clf; imshow(R); set(gcf,'name',fn_out)
        end
    end

    alpha = 1;
    beta = 1;
    R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain);
    fn_out = sprintf('results_hdr_large/%s.png',fnlist(F).name(1:end-4))
    imwrite(R,fn_out)
end

%% generate supplementary material: inverse tone mapping

sigma_r = 2.5;
colorRemapping = 'lum';
domain = 'log';
alpha = 1;
beta = 2.5;

srcdir = 'input_png/';
fnlist = dir([srcdir,'*.png']);
mkdir('results_ihdr')
for F = 1:length(fnlist)
    fn_in = [srcdir,fnlist(F).name]
    I = (double( imread(fn_in) )/255) .^ 2.2;  % apply gamma to linearize
    %figure; imshow(I); set(gcf,'name',fn_in)
    
	R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain);
    fn_out = sprintf('results_ihdr/%s_%s_%s_s%g_a%g_b%g.hdr',fnlist(F).name(1:end-4),colorRemapping,domain,sigma_r,alpha,beta)
    hdrwrite(R,fn_out)
end
