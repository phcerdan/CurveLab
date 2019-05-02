
% crvlt_levelseterror.m
%
% Thresholds the curvelet transform of a simple image and sets evaluates
% the reconstruction at the 0.5-level set of the original image, plotting
% the results
%

% define grid
N=512;
x=0:1/N:1-10*eps;
[X,Y]=meshgrid(x,x);

% different shapes
front = double(X < 0.5);

Rsm=0.1;
Rlg=0.35;
circsm = double((X-0.5).^2 + (Y-0.5).^2 < Rsm^2);
circlg = double((X-0.5).^2 + (Y-0.5).^2 < Rlg^2);

drop = eggshape(2, N, [0.5; 0.35], [0.55; 0.7], 0.25, 0.01);

inpaper = exp(-20*((X-0.5).^2 + (Y-0.5).^2));
inpaper(X.^2 + (Y-1).^2 < 0.66^2) = 0.1 * inpaper(X.^2 + (Y-1).^2 < 0.66^2);

% choose a shape
im = drop;

% smooth the function
sig2 = (1/N)^2;
gauss=exp(-((X-0.5).^2 + (Y-0.5).^2)/sig2);
ims=ifftshift(ifft2(fft2(gauss).*fft2(im))/sum(gauss(:)));

% compute curvelets
C = fdct_wrapping(ims, 0);

% select levels, and threshold
levs = 1:2;
Ccbc2 = crvlt_getlevels(C, levs);
Cth = crvlt_thresh(C, 0.1);
imlev = real(ifdct_wrapping(Ccbc2, 0));
imth = real(ifdct_wrapping(Cth, 0));

% find 0.5-level curve
cont = contourc(x,x,imlev, [0.5 0.5]);
cont = cont(:, 2:end);  % ignore level and size info
contth = contourc(x,x,imth, [0.5 0.5]);
contth = contth(:, 2:end);  % ignore level and size info

% evaluate original function on level-curve points
cvals = interp2(X, Y, ims, cont(1,:), cont(2,:), 'linear') - 0.5;
cvalsth = interp2(X, Y, ims, contth(1,:), contth(2,:), 'linear') - 0.5;
figure(1), clf
plot3(cont(1,:), cont(2,:), cvals, 'r');
axis([0 1 0 1 min(cvals) max(cvals)+eps]), grid on
hold on
plot3(contth(1,:), contth(2,:), cvalsth, 'b');

figure(2), clf
subplot(1,3,1)
imagesc(x,x,ims), set(gca,'YDir','normal'), colormap gray, axis equal tight, title('original')
subplot(1,3,2)
imagesc(x,x,imlev), set(gca,'YDir','normal'), colormap gray, axis equal tight, title(['reconstruction, levels ', num2str(levs)])
subplot(1,3,3)
imagesc(x,x,imth), set(gca,'YDir','normal'), colormap gray, axis equal tight, title(sprintf('reconstruction, %d largest coeffs', crvlt_countnnz(Cth)))
plotcurveletpos(Cth)


disp(sprintf('Curv: Levelset-err: %f, L2-err: %f, Ncoeffs: %d', max(abs(cvals)), norm(imlev(:)-ims(:))/sqrt(prod(size(ims))), crvlt_countnnz(Ccbc2)))
disp(sprintf('Thrsh: Levelset-err: %f, L2-err: %f, Ncoeffs: %d', max(abs(cvalsth)), norm(imth(:)-ims(:))/sqrt(prod(size(ims))), crvlt_countnnz(Cth)))
