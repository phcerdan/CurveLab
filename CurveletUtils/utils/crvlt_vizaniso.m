

function img = crvlt_vizaniso(A)

% IMG = CRVLT_VIZANISO(A)
%
%   visualizes the anisotropy measure A computed by crvlt_anisomeas
%   and returns the result in the image IMG. The visualization is based on
%   a circular decomposition similar to the decomposition of Fourier space
%   in the curvelet transform.
%
% EXAMPLES
%   To compute the anisotropy measure and display it
%       img = imread('Lena.jpg');
%       C = fdct_wrapping(img, 0);
%       A = crvlt_anisomeas(C);
%       aimg = vizaniso(A);
%       imagesc(aimg); axis equal tight; axis off; colormap bone
%


dr = 1;
xmax = dr * length(A);
x = linspace(-xmax,xmax,20*length(A)+1);
[X,Y] = meshgrid(x,x);
R2 = X.^2 + Y.^2;
TH = atan2(Y,X);
img = zeros(size(R2));

for k=1:length(A),
    ndirs = length(A{k});
    if ndirs == 1,
        hiang = pi;
        lowang = -pi;
    else
        hiang = 3*pi/4 - (0:ndirs-1) * 2*pi/ndirs;
        lowang = hiang - 2*pi/ndirs;
        hiang = forceinpipi(hiang);
        lowang = forceinpipi(lowang);
        hiang(lowang > hiang) = hiang(lowang > hiang) + 2*pi;
    end
    
    for l=1:ndirs,
        img(R2 >= ((k-1)*dr)^2 & R2 < (k*dr)^2 & TH >= lowang(l)-10*eps & TH < hiang(l)+10*eps) = A{k}(l);
    end
end

img = flipud(img);
    
    
   
%% Force angles into [-pi, pi) range
function ang = forceinpipi(ang)

ang = mod(ang + 3*pi, 2*pi) - pi;
