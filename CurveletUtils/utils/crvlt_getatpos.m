function [Ct, cfs] = crvlt_getatpos(C, x, y, lev)

% [Ct, cfs] = crvlt_getatpos(C, x, y, lev)
%
% Returns curvelet coefficients at positions in 
% vectors x,y (assuming domain is [0,1) x [0,1))
% at levels in vector lev
% Ct has same dimensions as C, with non-zeros at
% specified places
% cfs gives the coefficients as a vector of size(x)
% for each level and direction 


if (nargin < 4)
  lev = 1:length(C);
end

Ct = zerodct(C);
cfs = cell(1,length(x));
for p=1:length(x)
  cfs{p} = cell(1,length(C));
end
for j=lev,
    for p=1:length(x)
      cfs{p}{j} = zeros(1,length(C{j}));
    end
    for k=1:length(C{j});
        xmap = round(x * size(C{j}{k},2));
        ymap = round(y * size(C{j}{k},1));
        Ct{j}{k}(ymap,xmap) = C{j}{k}(ymap,xmap);
        for p=1:length(x)
            cfs{p}{j}(k) = Ct{j}{k}(ymap(p),xmap(p));
        end
    end
end
