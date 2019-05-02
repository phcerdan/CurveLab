function Cz=crvlt_zeros(C)

% Cz = crvlt_zeros(C)
%
% Creates an all-zero curvelet coefficient cell array with the same sizes
% as the given cell array
%

Cz=cell(size(C));
for k=1:length(Cz),
    Cz{k}=cell(size(C{k}));
    for j=1:length(Cz{k}),
        Cz{k}{j} = zeros(size(C{k}{j}));
    end
end

        
