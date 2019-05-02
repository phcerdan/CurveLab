
% Example of the use of crvlt_extractdir and curvecanny_multi
%
% Tries to extract edges in the loaded image.
%

%im = imread('../Curvelets/Images/dmso b 1.jpg');
%im = double(im(:,:,1));

disp('Applying curvelet transform...')
tic
C = fdct_wrapping(im, 0);
toc

% Parameters for direction-field extraction
cavg = 1;   % how many steps away to average directions
crvltsize = 1; % how many steps away to sum in space
nrdirfld = 1; % nr of directional fields to use (largest magnitude, plus next to largest, plus...)
dirspace = 3; % minimum distance (in nr of directions) between selected dirs (if nrdirfld > 1)

levs = length(C)-2:length(C)-1;  % which levels to use

% extract the directional field
disp('Extracting directional field...')
tic
fld = crvlt_extractdirs(C, levs, nrdirfld, cavg, crvltsize, dirspace);
mag = sqrt(fld{1}{1}.^2 + fld{1}{2}.^2);
toc
    
% Parameters for canny algorithm and postprocessing
thlow = 0.25;  % threshold for weak pixels
thhigh = 0.35; % threshold for strong pixels
extlen = 5;    % number of pixels to extend along direction
extth = [0.3 0.7]; % thresholds for valid pixels for edge extension [fraction of max, fraction of max at pixel]
thinedge = 0;  % set to 1 to thin edges
dilatethin = 0; % set to 1 to dilate and then thin edges (to connect close-by components)
despur = 0;     % set to 1 to remove spurs (=tangling edges)

% get the edges
disp('Extracting edges...') 
tic
e = curvecanny_multi(fld, [thlow thhigh], extlen, extth, thinedge, dilatethin, despur);
%e = e>=2;
toc

% interpolate to the image grid
E = interp2(linspace(0,1,size(e,2))', linspace(0,1,size(e,1)), double(e), ...
    linspace(0,1,size(im,2))', linspace(0,1,size(im,1)), 'nearest');
E = E > 0.5;
% E = bwmorph(E, 'thin', 3);  % uncomment to thin the edges on image grid


% plot the edges
figure(1), clf
subplot(1,2,1)
imagesc(im), colormap gray; axis equal tight;
hold on
[Ex,Ey] = find(E);
plot(Ey,Ex,'g.')
title('Image and detected edges')

subplot(1,2,2)
imagesc(mag), axis equal tight;
hold on
quiver(fld{1}{1}, -fld{1}{2}, 'r')   % plot arrows (minus on y-component because imagesc reverses y-axis)
title('Magnitude and direction of maximal curvelet coefficient')
