function [nrnz,nel]=countnnz(C, lev)

% Counts the number of non-zero coefficients, and the total number of
% coefficients in the curvelet-coefficient cell array C

if nargin < 2,
  lev = 1:length(C);
end

nrnz=0;
nel=0;
for k=lev,
    for j=1:length(C{k}),
        nrnz = nrnz + nnz(C{k}{j});
        nel = nel + prod(size(C{k}{j}));
    end
end
