function Csc = crvlt_rescale(C, sc)

% Csc = crvlt_rescale(C, sc)
%
% multiplies all curvelet coefficients in C by the factor sc,
% returning the result in Csc
%

Csc = C;
for j=1:length(C),
    for l=1:length(C{j}),
        Csc{j}{l} = sc * Csc{j}{l};
    end
end


   
