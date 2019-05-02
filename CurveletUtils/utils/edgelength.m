function len=edgelength(edges)

% Measures the length of a connected component, specified by 1's in the
% black/white image edges.

len = 0;
sq2 = sqrt(2);
iedges = int32(edges);
iedges([1 end]) = 0;  % discard first and last element, because of Csel selection below

I = find(iedges);
nrows = size(iedges, 1);
conn4 = [-1 1 -nrows nrows];
conn8 = [-nrows-1 -nrows+1 nrows-1 nrows+1];

newI = I(1);
while ~isempty(newI),
    iedges(newI) = 2;   % mark as already visited
 
    % find axially and diagonally connected pixels
    CONN4 = repmat(conn4, [length(newI) 1]);
    CONN8 = repmat(conn8, [length(newI) 1]);
    NEWI = repmat(newI', [1 4]);

    Csel = iedges(min(max(NEWI + CONN4,1), numel(iedges))) == 1;
    I4 = NEWI(Csel) + CONN4(Csel);
    I4 = I4(:)';

    Csel = iedges(min(max(NEWI + CONN8,1), numel(iedges))) == 1;
    I8 = NEWI(Csel) + CONN8(Csel);
    I8 = I8(:)';
    I8 = I8(~ismember(I8,I4));  % exclude doublets
  
    newI = unique([I4 I8]);
    
    len = len + 1*length(I4);
    len = len + sq2*length(I8);
end

