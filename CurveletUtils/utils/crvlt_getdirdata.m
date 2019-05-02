function [dirs,vals] = crvlt_getdirdata(C, pos, levels)

% [dirs, vals] = crvlt_getdirdata(C, pos, levels)
%
% Gets coefficients for all directions at pixel pos = [y x]
% and on specified levels of curvelet data C.
%
% dirs is a cell array of the same length as levels, where
% 	dirs{k} is an array from 1 to the number of directions
% 	on level levels{k}
%
% vals is a cell array which for each specified level contains
%	the coefficients corresponding to dirs

dirs = cell(1, length(levels));
vals = dirs;
for k=1:length(levels)
    dirs{k} = 1:length(C{levels(k)});
    v = zeros(1,length(C{levels(k)}));
    for dir = dirs{k},
        idx = maptolev(pos, size(C{end}{1}), C, levels(k), dir);
        v(dir) = C{levels(k)}{dir}(idx(1) + 1,idx(2) + 1);
    end
    vals{k} = v;
end




function idx = maptolev(pos, dims, C, lev, dir)

outdim = size(C{lev}{dir});
idx = round((pos - [0.5 0.5])./dims.*outdim);

