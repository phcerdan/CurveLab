
% Shows the directional field extracted by crvlt_extractdirs
% on top of the magnitude of the field

im = imread('Lena.jpg');

disp('Applying curvelet transform...')
[m,n]=size(im);
tic
C = fdct_wrapping_mex(m, n, 5, 16, 0, double(im));
toc

disp('Extracting direction field...')
lev = 3;  % specify level(s) to extract
tic
fld = crvlt_extractdirs(C, lev, 1);
toc

mag = sqrt(fld{1}{1}.^2 + fld{1}{2}.^2);

figure(1), clf
subplot(1,2,1)
imagesc(im), axis image, colormap gray
subplot(1,2,2)
imagesc(mag), axis equal tight, colormap gray, colorbar, hold on
quiver(fld{1}{1}, -fld{1}{2})  % minus on y-component since imagesc() reverses y-axis


