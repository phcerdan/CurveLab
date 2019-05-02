function [y,n]=crvlt_getnorms(C)

% [y,n] = crvlt_getnorms(C)
%
% gets the l^2-norm of the curvelet coefficients at each
% level
% y is a vector with l^2-norms
% n is a vector with the level numbers

y=zeros(1,length(C));
n=1:length(C);
for j=1:length(C),
    for l=1:length(C{j}),
        y(j) = y(j) + sum(abs(C{j}{l}(:)).^2);
    end
    y(j) = sqrt(y(j));
end

