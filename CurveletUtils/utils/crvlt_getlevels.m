function Co = crvlt_getlevels(C, lev)

%  Co = crvlt_getlevels(C, lev)
%  
%  extracts the levels specified in vector lev from
%  curvelet data C, by setting all other coefficients to zero.
%  Co has same dimensions as C, but has only zeros on levels not
%  in lev.
%

Co = C;
for ll=setdiff(1:length(C), lev),
    for k=1:length(C{ll}),
        Co{ll}{k} = zeros(size(Co{ll}{k}));
    end
end

