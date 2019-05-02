function [n,yl2,yls,cs]=crvlt_errorplot(shape, smooth, wavelet)

% [n,yl2,yls,cs] = crvlt_errorplot(shape, smooth, wavelet)
%
% Computes subsequent reconstructions using different numbers of 
% coefficients, using either curvelets or wavelets, and returns
% n     vector with the number of coefficients
% yl2   the l^2-errors of the reconstructions
% yls   the relative l^2-error on the 0.5-level set of the image
% cs    the coefficient sizes at positions in n
%
% Input arguments:
% shape     string that determines the image to use. May take the values
%           'drop', 'smcirc', 'lgcirc', 'ring', 'front', 'Lena', 'inpaper'
% smooth    if true, then the image is slightly smoothed by a Gaussian
%           filter
% wavelet   if true, use 'db3' wavelets, otherwise use curvelets
%
% If no output arguments are given, the results are plotted.
%


if nargin==0,
    shape = 'drop';
    smooth = 1;
    wavelet=0;
end

N = 1024;
x=0:1/N:1-10*eps;
[X,Y]=meshgrid(x,x);

switch shape,
    case 'drop' 
        im = eggshape(2, N, [0.5; 0.35], [0.55; 0.7], 0.25, 0.01);
    case 'smcirc'
        Rsm = 0.1;
        im = double((X-0.5).^2 + (Y-0.5).^2 < Rsm^2);
    case 'lgcirc'
        Rlg = 0.4;
        im = double((X-0.5).^2 + (Y-0.5).^2 < Rlg^2);
    case 'ring'
        Rsm = 0.1;
        im0 = double((X-0.5).^2 + (Y-0.5).^2 < Rsm^2);
        im = ([diff(im0,1,1); zeros(1,N)] + [diff(im0,1,2)  zeros(N,1)]);
        im = double(abs(im) > 0);
    case 'front'
        im = double(X < 0.5);
    case 'Lena'
        im = double(imread('Lena.jpg'));
        N = size(im,1);
        x=0:1/N:1-10*eps;
        [X,Y]=meshgrid(x,x);
    case 'inpaper'
        im = exp(-20*((X-0.5).^2 + (Y-0.5).^2));
        im(X.^2 + (Y-1).^2 < 0.66^2) = 0.1 * im(X.^2 + (Y-1).^2 < 0.66^2);
    otherwise
        error('Unknown shape type')
end

% smooth the function
if smooth,
    sig2 = (2/N)^2;
    gauss=exp(-((X-0.5).^2 + (Y-0.5).^2)/sig2);
    ims=ifftshift(ifft2(fft2(gauss).*fft2(im))/sum(gauss(:)));
else
    ims = im;
end
    
wtype = 'db3';
Nwav = 3;
dwtmode('per');  % set to peridodic extension
[D, S] = wavedec2(ims, Nwav, wtype);
if wavelet,
    ds = sort(abs(D),2,'descend');
else
    C = fdct_wrapping(ims, 0);
    % Choose the same nr of coeffs as in wavelet case
    % and rescale to 'compensate for redundancy'
    Cth = crvlt_keeplargest(C, length(D));
    oldnrm = norm(ims(:));
    newnrm = crvlt_norm(Cth);
    C = crvlt_rescale(C, oldnrm/newnrm);
end

%thresholds = [0.5 0.2 0.1 0.05 0.01 0.001 0.0001];
nterms = [10 50 100 500 1000 5000 20000 100000 1000000];
imnorm = norm(ims(:));
immaxmin = max(ims(:)) - min(ims(:));

disp('.')

n = [];
yl2 = [];
yls = [];
psnr = [];
cs = [];
%for thrsh=thresholds,
for nt = nterms,
    if wavelet,
        thrsh = ds(nt+1);
        %Dth = wthcoef2('t', D, S, 1:Nwav, thrsh*ones(1,Nwav), 'h');
        Dth = wthresh(D,'h',thrsh);
        n = [n nnz(Dth)];
        cs = ds(n);
        imth = waverec2(Dth,S,wtype);
    else
        %Cth = crvlt_thresh(C, thrsh);  % Threshold coefficients
        [Cth, thrsh] = crvlt_keeplargest(C, nt);
        n = [n crvlt_countnnz(Cth)];
        cs = [cs thrsh];
        imth = real(ifdct_wrapping(Cth, 0));
    end
    
    % Save L2-error
    yl2 = [yl2 norm(imth(:) - ims(:))/imnorm];
    psnr = [psnr 20 * log10(immaxmin/(norm(imth(:) - ims(:))/sqrt(numel(ims))))];
    
    % Save level-curve position error
    % find 0.5-level curve
    cont = contourc(x,x,imth, [0.5 0.5]);
    if size(cont, 2) > 1,
        cont = cont(:, 2:cont(2,1));  % ignore level and size info, only get first curve
        % evaluate original function on level-curve points
        cvals = interp2(X, Y, ims, cont(1,:), cont(2,:), 'linear') - 0.5;
        yls = [yls norm(cvals(:))/(0.5*sqrt(length(cont)))];
    else
       yls = [yls NaN]; 
    end
end
   
% plot results
if nargout==0,
    if wavelet,
        linestyle = 'ko-';
    else
        linestyle = 'r*-';
    end
    figure(1)
    loglog(n, psnr, linestyle)
    title('PSNR value depending on number of coefficients')
    xlabel('Number of coefficients')
    ylabel('PSNR')
    hold on
    
    figure(2)
    loglog(n, yls, linestyle)
    title('relative L2 error on level set')
    xlabel('Number of coefficients')
    ylabel('relative L2-error on level set')
    hold on
    plot([N^2 N^2],[1 1e-5],'b--')
    loglog(n, 10^2./n, 'b--')
    loglog(n, 10^4./n.^2, 'b--')
    
    figure(3), clf
    subplot(1,3,1)
    imagesc(x,x,ims), set(gca,'YDir','normal'), axis equal tight, colormap gray
    subplot(1,3,2)
    plot3(cont(1,:), cont(2,:), real(cvals), 'r'), grid on
    subplot(1,3,3)
    mesh(X,Y,real(imth-ims)), shading flat
    
    figure(4)
    loglog(n,cs,linestyle)
    title('Coefficient magnitudes')
    xlabel('number of coefficients')
    ylabel('coefficient magnitude')
    hold on
end

    
