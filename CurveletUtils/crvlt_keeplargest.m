% Help file for CRVLT_KEEPLARGEST mex-file
%
%   [CMOD, THVAL] = CRVLT_KEEPLARGEST(C, NCOEFFS)
%
%   If NCOEFFS >= 1:
%   Retains the NCOEFFS largest (in magnitude) curvelet coefficients in the
%   curvelet data structure C, setting the other coefficients to zero.
%   CMOD has the same dimensions as C, and NCOEFFS non-zero entries.
%   THVAL contains the threshold value, i.e. the magnitude of the smallest 
%   non-zero coefficient in CMOD.
%  
%   If 0 <= NCOEFFS < 1:
%   The value NCOEFFS is interpreted as a fraction of the non-zero entries in
%   NCOEFFS*100 percent of the non-zero entries are
%   zero.
%
%
% EXAMPLES:
%
%   To compute the reconstruction of the image 'Lena.jpg' using the
%   1000 largest coefficients.
%     im = imread('Lena.jpg');
%     C = fdct_wrapping(im, 0);
%     Cmod = crvlt_keeplargest(C, 1000);
%     immod = idfct_wrapping(Cmod, 0);
%     imagesc(immod), axis image, colormap gray 
%	
%


