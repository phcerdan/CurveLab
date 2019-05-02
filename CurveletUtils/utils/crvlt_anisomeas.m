function A = crvlt_anisomeas(C)

%
% A = CRVLT_ANISOMEAS(C)
%
%   computes the space-scale anisotrophy measure from the curvelet data C.
%   This measure is defined by 
%   A(j,l) = sum_k |c_{j,l,k}|^2 / < sum_k |c_{j,l,k}|^2 >_l
%   where < >_l denotes the mean over angle indices l.
%   k is the spacial index and j the level index.
%
%   A is a cell array with one item per level in C, each item being a row
%   array with one entry per direction at that level.
%
%

A = cell(size(C));
for j=1:length(C),
    A{j} = zeros(1,length(C{j}));
    for l=1:length(C{j}),
        A{j}(l) = norm(C{j}{l}(:))^2;
    end
    A{j} = A{j} / sum(A{j}) * length(A{j});
end
