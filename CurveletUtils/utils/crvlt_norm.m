function nrm = crvlt_norm(C)

% nrm = crvlt_norm(C)
%
% Returns the l^2-norm of the coefficients in C, i.e.
% sqrt(sum_jlk |c_jlk|^2)
%

nrm = 0;
for j=1:length(C),
    for l=1:length(C{j}),
        nrm = nrm + sum(abs(C{j}{l}(:)).^2);
    end
end
nrm = sqrt(nrm);

