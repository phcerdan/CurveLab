function plotwaveletpos(D, S, axh, marksize)

% plotwaveletpos(D, S, axh, marksize)
%
% Attempts to plot the wavelet positions given by the wavelet data [D,S]
% in the axes 'axh', using dots of size 'marksize'
% 
% This function has NOT been tested for all wavelet types and the results
% are a bit dubious. Assumes that the wavelet coefficients were computed
% with dwtmode('ppd'), i.e. periodic images.
%

if nargin < 2,
    error('Too few input arguments!')
elseif nargin==2,
    axh = gca;
    marksize = 6;
elseif nargin == 3,
    marksize = 1;
end

nlevs = size(S,1)-1;
Dlev = cell(1, nlevs);
levlength = prod(S,2);
ofs = [1 3*ones(1,length(levlength)-1)].*levlength';  % account for horiz, vert, diag on detail levels
ofs = [0 cumsum(ofs)];

color={'b.','r.','c.','g.'};
for j=5:nlevs, color{j}='m.'; end

axes(axh)
hold on

for k=1:nlevs,
    Dlev{k} = abs(reshape(D(ofs(k)+(1:levlength(k))), S(k,:))); 
    if k > 1,
        for j=1:2,
            Dlev{k} = Dlev{k} + abs(reshape(D(ofs(k)+j*levlength(k) + (1:levlength(k))), S(k,:)));
        end
    end

    scale = (S(end,:)-1) ./ (S(k,:) - 1);
    xv = (0:S(k,2)-1) * scale(2) + 1;
    yv = (0:S(k,1)-1) * scale(1) + 1;
    [X,Y] = meshgrid(xv,yv);
    I = Dlev{k} > 0;
    plot(axh,X(I),Y(I),color{k},'MarkerSize',marksize)
end
