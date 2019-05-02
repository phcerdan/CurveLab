% Help for the CRVLT_EXTRACTDIRS mex-file.
%
% FIELDS = CRVLT_EXTRACTDIRS(C, LEVELS)
%   extracts two directional fields (pointing along the ridges of the curvelets, 
%   i.e. along the edges in the image) from the curvelet data in C, using only
%   data from the levels specified in the vector levels (one-based indices).
%  
%   The absolute values of the curvelet coefficients at each grid point (on the
%   grid defined by the finest level in C) are compared and the direction
%   corresponding to the largest coefficient is selected and returned as a 
%   vector with x- and y-components, with magnitude equal to the magnitude of
%   the selected coefficient.
%   Then the second largest coefficient is selected and its direction vector
%   saved in the second field, in case it is more than three direction steps from 
%   the first selected direction and is a local maximum among the directions. If 
%   no such direction is found, a zero-vector is returned.
%   The coefficients from the selected levels are added togther (mapped to the 
%   finest grid), and the coefficients from neighboring directions and grid
%   points are also added to provide some smoothing.
%
%   The 'FIELDS' data structure is a cell array, where each item contains one 
%   directional field. Each directional field is a cell array with two 
%   components, the x- and y-components of the field, each being a
%   two-dimensional
%   array of the same size as the finest grid size of the specified levels.
%   Thus, 'FIELDS{2}{2}(:,:)' gives the y-component of the second field. 
% 
%
% FIELDS = CRVLT_EXTRACTDIRS(C, LEVELS, NFIELDS, CSUMSIZE, CRVLTSIZE, DIRGAP) 
%   sets additional parameters as follows:
%   NFIELDS 	specifies the number of fields to extract (default 2)
%   CSUMSIZE	specifies the number of directions left and right to sum over
%  		before selecting maximum (default 1)
%   CRVLTSIZE	sets the number of grid points away that a curvelet influences
%  		(on its own grid), meaning that this many coefficient magnitudes
%   		in all directions are added to the current coefficient before 
%  		selecting maximum (default 1).
%   DIRGAP	sets the required spacing between subsequent maximum direction
%  		selections (default 3)
%   
%   Empty matrices may be passed as argument where default parameter values are
%   desired.
%
%
% EXAMPLES:
%
%   To compute two directional fields using levels 3 and 4 from the curvelet
%   coefficients computed from the image 'Lena.jpg', and plot the 
%   magnitude of the field and then the directional field as arrows on top.
%     im = imread('Lena.jpg);
%     C = fdct_wrapping(im, 0);
%     fields = crvlt_extractdirs(C, 3:4)
%     mag = sqrt(fld{1}{1}.^2 + fld{1}{2}.^2);
%     imagesc(mag), axis image, colormap gray, hold on
%     quiver(fields{1}{1}, -fields{1}{2})
%
%   Note the minus sign for the second component, due to the reversed 
%   y-axis produced by imagesc().
%
%
%

