function csort=crvlt_sortcoeffs(C)

% csort = crvlt_sortcoeffs(C)
%
% Sorts all coefficients in curvelet data structure 'C' in descending order
% by their magnitude and returns the magnitudes in a column vector.
%


csort = [];
for j=1:length(C),
    for l=1:length(C{j}),
        csort = [csort; abs(C{j}{l}(:))];
    end
end

csort = sort(csort,1,'descend');

