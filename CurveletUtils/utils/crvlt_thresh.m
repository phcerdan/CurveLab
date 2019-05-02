function Cth = crvlt_thresh(C, th)

% Cth = crvlt_thresh(C, th)
%
% Performs hard thresholding on the curvelet coefficients 'C', using the
% threshold 'th'.
% 

Cth=C;
for j=1:length(Cth)
    for l=1:length(Cth{j}),
        Cth{j}{l} = C{j}{l} .* (abs(C{j}{l}) > th);
    end
end
