function [eout,thresh,mag] = curvecanny_multi(varargin)

% [EOUT,THRESH,MAG] = CURVECANNY_MULTI(FLDS, THRSH, EXTLEN, EXTTH, THIN, DILTHIN, DESPUR)
% 
% Applies a Canny edge detector type algorithm to the field data extracted
% by CRVLT_EXTRACTDIRS, tracing along the ridges of the curvelet
% magnitude image. Makes use of multiple fields if possible.
% 
% Postprocessing tries to extend edges along field directions, and join
% close edges.
%
% Returns:
% EOUT      an image with value 1 on the found edges, and 0 elsewhere
%
% THRESH    the actual thresholds used (useful if default values are used)
%
% MAG       the magnitude of the field
%
%
% Optional parameters:
% THRSH     thresholds for the Canny algorithm (fractions of maximum image
%           intensity)
%           THRSH(1) is the 'weak' threshold and THRSH(2) the 'strong'
%           if length(THRSH)==1, this specifies the 'strong' threshold, and
%           the 'weak' is computed. If no threshold is given, the strong
%           threshold is computed, making a fraction of the pixels edges.
% 
% EXTLEN    the maximum distance to extend the edges
%
% EXTTH     thresholds for edge extension. 
%           EXTTH(1) is the threshold for valid pixels, as fraction of max
%           directional magnitude. 
%           EXTTH(2) is the fraction of the first direction magnitude at
%           each pixel, that subsequent magnitudes must be larger than.
%
% THIN      if nonzero, make the edges thinner (by THIN steps) after extending them
%
% DILTHIN   if nonzero, dilate and thin (by DILTHIN steps), in order to join closeby edges
%
% DESPUR    if nonzero, remove tangling edges (by DESPUR steps)
%


% parse input arguments

if nargin < 1,
  error('Too few input arguments')
end

a = varargin{1};
% default values
thresh = [];
extlen = 5;
extth = [0.3 0.7];
thinedge = 1;
dilatethin = 0;
despur = 0;
% set values if sent as arguments
switch nargin,
    case 2,
        if ~isempty(varargin{2}), thresh = varargin{2}; end
    case 3,
        if ~isempty(varargin{2}), thresh = varargin{2}; end
        if ~isempty(varargin{3}), extlen = varargin{3}; end
    case 4,
        if ~isempty(varargin{2}), thresh = varargin{2}; end
        if ~isempty(varargin{3}), extlen = varargin{3}; end
        if ~isempty(varargin{4}), extth = varargin{4}; end
    case 5,
        if ~isempty(varargin{2}), thresh = varargin{2}; end
        if ~isempty(varargin{3}), extlen = varargin{3}; end
        if ~isempty(varargin{4}), extth = varargin{4}; end
        if ~isempty(varargin{5}), thinedge = varargin{4}; end
    case 6,
        if ~isempty(varargin{2}), thresh = varargin{2}; end
        if ~isempty(varargin{3}), extlen = varargin{3}; end
        if ~isempty(varargin{4}), extth = varargin{4}; end
        if ~isempty(varargin{5}), thinedge = varargin{5}; end
        if ~isempty(varargin{6}), dilatethin = varargin{6}; end
    case 7,
        if ~isempty(varargin{2}), thresh = varargin{2}; end
        if ~isempty(varargin{3}), extlen = varargin{3}; end
        if ~isempty(varargin{4}), extth = varargin{4}; end
        if ~isempty(varargin{5}), thinedge = varargin{5}; end
        if ~isempty(varargin{6}), dilatethin = varargin{6}; end
        if ~isempty(varargin{7}), despur = varargin{7}; end
end



% Get sizes and compute basic values

[m,n]=size(a{1}{1});

% The output edge map:
e = false(m,n);

% Magic numbers
PercentOfPixelsNotEdges = .9; % Used for selecting thresholds
ThresholdRatio = .4;          % Low thresh is this fraction of the high.

mag = cell(1,length(a));
for k=1:length(a),
    mag{k} = sqrt(a{k}{1}.^2 + a{k}{2}.^2);
end
magmax = max(mag{1}(:));
if magmax>0
    for k=1:length(mag),
      mag{k} = mag{k} / magmax;   % normalize by first magnitudes
    end
end



% Select the thresholds for Canny step
if isempty(thresh)
  counts=imhist(mag{1}, 64);
  highThresh = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
    1,'first') / 64;
  lowThresh = ThresholdRatio*highThresh;
  thresh = [lowThresh highThresh];
elseif length(thresh)==1
  highThresh = thresh;
  if thresh>=1
    eid = sprintf('Images:%s:thresholdMustBeLessThanOne', mfilename);
    msg = 'The threshold must be less than 1.';
    error(eid,'%s',msg);
  end
  lowThresh = ThresholdRatio*thresh;
  thresh = [lowThresh highThresh];
elseif length(thresh)==2
  lowThresh = thresh(1);
  highThresh = thresh(2);
  if (lowThresh >= highThresh) || (highThresh >= 1)
    eid = sprintf('Images:%s:thresholdOutOfRange', mfilename);
    msg = 'Thresh must be [low high], where low < high < 1.';
    error(eid,'%s',msg);
  end
end

% Do the Canny step
% The next step is to do the non-maximum supression.
% We will accrue indices which specify ON pixels in strong edgemap
% The array e will become the weak edge map.
idxStrong = [];
for k=1,%k=1:length(mag),
    for dir = 1:4
        idxLocalMax = cannyFindLocalMaxima(dir,a{k}{1},a{k}{2},mag{k});
        idxWeak = idxLocalMax(mag{k}(idxLocalMax) > lowThresh);
        e(idxWeak)=1;
        idxStrong = [idxStrong; idxWeak(mag{k}(idxWeak) > highThresh)];
    end
end

es = [];
if ~isempty(idxStrong) % result is all zeros if idxStrong is empty
  rstrong = rem(idxStrong-1, m)+1;
  cstrong = floor((idxStrong-1)/m)+1;
  es = bwselect(e, cstrong, rstrong, 8);
  %es = bwmorph(es, 'thin', 1);  % Thin double (or triple) pixel wide contours
end
e = double(e);  % set weak pixels to 1
e(es>0) = 3;      % set (extended) strong pixels to 3


% extend edges along major directions
if extlen > 0,
    e = extendalongdir(e,a,mag,extlen, extth);
end
e = (e>1);

% thin and/or dilate and thin
if thinedge,
    e = bwmorph(e, 'thin', thinedge);
end
if dilatethin,
    e = bwmorph(e, 'dilate', 1);
    e = bwmorph(e, 'thin', dilatethin);
end
if despur,
    e = bwmorph(e, 'spur', despur);
    e = bwmorph(e, 'clean', 1);
end


if nargout==0,
  imshow(e);
else
  eout = e;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : cannyFindLocalMaxima
%
function idxLocalMax = cannyFindLocalMaxima(direction,ix,iy,mag)
%
% This sub-function helps with the non-maximum supression in the Canny
% edge detector.  The input parameters are:
% 
%   direction - the index of which direction the gradient is pointing, 
%               read from the diagram below. direction is 1, 2, 3, or 4.
%   ix        - input image filtered by derivative of gaussian along x 
%   iy        - input image filtered by derivative of gaussian along y
%   mag       - the gradient magnitude image
%
%    there are 4 cases:
%
%                         The X marks the pixel in question, and each
%         3     2         of the quadrants for the gradient vector
%       O----0----0       fall into two cases, divided by the 45 
%     4 |         | 1     degree line.  In one case the gradient
%       |         |       vector is more horizontal, and in the other
%       O    X    O       it is more vertical.  There are eight 
%       |         |       divisions, but for the non-maximum supression  
%    (1)|         |(4)    we are only worried about 4 of them since we 
%       O----O----O       use symmetric points about the center pixel.
%        (2)   (3)        


[m,n] = size(mag);

% Find the indices of all points whose gradient (specified by the 
% vector (ix,iy)) is going in the direction we're looking at.  

switch direction
 case 1
  idx = find((iy<=0 & ix>-iy)  | (iy>=0 & ix<-iy));
 case 2
  idx = find((ix>0 & -iy>=ix)  | (ix<0 & -iy<=ix));
 case 3
  idx = find((ix<=0 & ix>iy) | (ix>=0 & ix<iy));
 case 4
  idx = find((iy<0 & ix<=iy) | (iy>0 & ix>=iy));
end

% Exclude the exterior pixels
if ~isempty(idx)
  v = mod(idx,m);
  extIdx = find(v==1 | v==0 | v==2 | v==m-1 | idx<=2*m | (idx>(n-2)*m));
  idx(extIdx) = [];
end

ixv = ix(idx);  
iyv = iy(idx);   
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels
switch direction
 case 1
  d = abs(iyv./ixv);
  gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d; 
  gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d; 
 case 2
  d = abs(ixv./iyv);
  gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d; 
  gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d; 
 case 3
  d = abs(ixv./iyv);
  gradmag1 = mag(idx+m).*(1-d) + mag(idx+m+1).*d; 
  gradmag2 = mag(idx-m).*(1-d) + mag(idx-m-1).*d; 
 case 4
  d = abs(iyv./ixv);
  gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d; 
  gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d; 
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2); 




function e = extendalongdir(e, flds, mags, extlen, extth)
% extend the pixels marked in e along directions in flds

angle = cell(size(flds));
valid = cell(size(flds));
th = extth(1) * max(mags{1}(:));
for m=1:length(flds),
    angle{m} = atan2(flds{m}{2},flds{m}{1});  % angles in [-pi,pi]
    valid{m} = (mags{m} > th) & (mags{m} > extth(2) * mags{1});
end

wstep = size(e,1);
niter = extlen;
dirmap = zeros(numel(e),8);

[elbl, Nc] = bwlabel(e,8);
cmpmap = zeros(numel(e), Nc);
flbl = find(elbl(:));
cmpmap(flbl + numel(e)*(elbl(flbl)-1)) = 1;

vhhalfangle = pi/6; % pi/8;
diaghalfangle = pi/6;

eskel = bwmorph(e==3, 'thin', Inf);
ef = imfilter(int32(eskel), ones(3,3));
dsum = (ef==2 & eskel);  % select only end of line segments
% extend in directions
for k=1:niter,
    for m=1:length(flds),
        % find directions
        ns = find(((angle{m} > pi/2-vhhalfangle & angle{m} <= pi/2 + vhhalfangle) | (angle{m} > -pi/2-vhhalfangle & angle{m} <= -pi/2+vhhalfangle)) & dsum & valid{m});
        ew = find(((angle{m} > -vhhalfangle & angle{m} <= vhhalfangle) | (angle{m} > pi-vhhalfangle | angle{m} <= -pi+vhhalfangle)) & dsum & valid{m});
        nwse = find(((angle{m} > 3*pi/4-diaghalfangle & angle{m} <= 3*pi/4+diaghalfangle) | (angle{m} > -pi/4-diaghalfangle & angle{m} <= -pi/4+diaghalfangle)) & dsum & valid{m});
        nesw = find(((angle{m} > pi/4-diaghalfangle & angle{m} <= pi/4+diaghalfangle) | (angle{m} > -3*pi/4-diaghalfangle & angle{m} <= -3*pi/4+diaghalfangle)) & dsum & valid{m});

        % assert that we stay in the image
        ns = ns(mod(ns,size(e,1))>1 & mod(ns,size(e,1)) < size(e,1));
        ew = ew(ew>size(e,1) & (ew < numel(e) - size(e,1)));
        nwse = nwse(nwse>size(e,1) & (nwse < numel(e) - size(e,1)) & mod(nwse,size(e,1))>1 & mod(nwse,size(e,1)) < size(e,1));
        nesw = nesw(nesw>size(e,1) & (nesw < numel(e) - size(e,1)) & mod(nesw,size(e,1))>1 & mod(nesw,size(e,1)) < size(e,1));

        % set directions that we went to get to new points, allowing only
        % "forward" movement
        curval = niter+1-k;
        for j = 1:length(ns),
            ldir = find(dirmap(ns(j),:) == curval+1);
            if any(ismember(ldir, [1 2 3 7 8])) || k==1,
                dirmap(ns(j)-1,1) = max(dirmap(ns(j)-1,1),curval);
                cmpmap(ns(j)-1, :) = max(cmpmap(ns(j)-1,:), cmpmap(ns(j),:));
            end
            if any(ismember(ldir, [3 4 5 6 7])) || k==1,
                dirmap(ns(j)+1,5) = max(dirmap(ns(j)+1,5),curval);
                cmpmap(ns(j)+1, :) = max(cmpmap(ns(j)+1,:), cmpmap(ns(j),:));
            end
        end
        for j = 1:length(ew),
            ldir = find(dirmap(ew(j),:) == curval+1);
            if any(ismember(ldir, [1 2 3 4 5])) || k==1,
                dirmap(ew(j)+wstep,3) = max(dirmap(ew(j)+wstep,3),curval);
                cmpmap(ew(j)+wstep, :) = max(cmpmap(ew(j)+wstep,:), cmpmap(ew(j),:));
            end
            if any(ismember(ldir, [5 6 7 8 1])) || k==1,
                dirmap(ew(j)-wstep,7) = max(dirmap(ew(j)-wstep,7),curval);
                cmpmap(ew(j)-wstep, :) = max(cmpmap(ew(j)-wstep,:), cmpmap(ew(j),:));
            end
        end
        for j=1:length(nwse),
            ldir = find(dirmap(nwse(j),:) == curval+1);
            if any(ismember(ldir, [6 7 8 1 2])) || k==1,
                dirmap(nwse(j)-1-wstep,8) = max(dirmap(nwse(j)-1-wstep,8), curval);
                cmpmap(nwse(j)-1-wstep, :) = max(cmpmap(nwse(j)-1-wstep,:), cmpmap(nwse(j),:));
            end
            if any(ismember(ldir, [2 3 4 5 6])) || k==1,
                dirmap(nwse(j)+1+wstep,4) = max(dirmap(nwse(j)+1+wstep,4), curval);
                cmpmap(nwse(j)+1+wstep, :) = max(cmpmap(nwse(j)+1+wstep,:), cmpmap(nwse(j),:));
            end
        end
        for j=1:length(nesw),
            ldir = find(dirmap(nesw(j),:) == curval+1);
            if any(ismember(ldir, [8 1 2 3 4])) || k==1,
                dirmap(nesw(j)-1+wstep,2) = max(dirmap(nesw(j)-1+wstep,2), curval);
                cmpmap(nesw(j)-1+wstep, :) = max(cmpmap(nesw(j)-1+wstep,:), cmpmap(nesw(j),:));
            end
            if any(ismember(ldir, [4 5 6 7 8])) || k==1,
                dirmap(nesw(j)+1-wstep,6) = max(dirmap(nesw(j)+1-wstep,6), curval);
                cmpmap(nesw(j)+1-wstep, :) = max(cmpmap(nesw(j)+1-wstep,:), cmpmap(nesw(j),:));
            end
        end
        
        dsum = reshape(any(dirmap==curval,2),size(e));% .* (e==0); %dsum + int32(reshape(sum(dirmap,2),size(dsum)) > 0);
    end
end

exval = 2;
% trace back where we hit something
dd = reshape(max(dirmap,[],2) .* (sum(cmpmap, 2) >= 2),size(e));
%e(dd > 0 & e>=2) = exval;
e(dd > 0) = exval;
for k=1:niter,
    pixchoice = (e(:) == exval);
    s = find(dirmap(:,1)==k & pixchoice);
    n = find(dirmap(:,5)==k & pixchoice);
    ea = find(dirmap(:,7)==k & pixchoice);
    w = find(dirmap(:,3)==k & pixchoice);
    se = find(dirmap(:,8)==k & pixchoice);
    nw = find(dirmap(:,4)==k & pixchoice);
    sw = find(dirmap(:,2)==k & pixchoice);
    ne = find(dirmap(:,6)==k & pixchoice);
   
    e(s+1) = max(e(s+1),exval);
    e(n-1) = max(e(n-1),exval);
    e(ea+wstep) = max(e(ea+wstep),exval);
    e(w-wstep) = max(e(w-wstep),exval);
    e(se+1+wstep) = max(e(se+1+wstep),exval);
    e(nw-1-wstep) = max(e(nw-1-wstep),exval);
    e(sw+1-wstep) = max(e(sw+1-wstep),exval);
    e(ne-1+wstep) = max(e(ne-1+wstep),exval);
end

return


