function nz = crvlt_nnzperlev(C)

% nz = crvlt_nnzperlev(C)
%
% returns the number of non-zeros on each level
% of the curvelet data in C
% nz is a vector of the same length as the number of levels in C
%

nz = zeros(length(C),1);
for j=1:length(C),
    for l=1:length(C{j}),
        nz(j) = nz(j) + nnz(C{j}{l});
    end
end

