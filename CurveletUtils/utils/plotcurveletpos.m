function plotcurveletpos(C, axh, arrows, marksize)

% if arrows is nonzero, then marksize means scale arrows or not

if nargin<2
  axh=gca;
  arrows = 1;
  marksize = 1;
elseif nargin==2,
  arrows = 1;
  marksize = 1;
elseif nargin==3,
  marksize = 6;
end

if arrows
    scalearrows = (marksize > 0);
    marksize = 6;
    autoscaling = 1;
end

axes(axh);
xlim = get(axh, 'XLim');
ylim = get(axh, 'YLim');
hold on
if arrows,
    color={'b.','rx','cx','gx'};
    for j=5:length(C), color{j}='mx'; end
else
    color={'b.','r.','c.','g.'};
    for j=5:length(C), color{j}='m.'; end
end

if (ndims(C{1}{1}) == 2)
    x = ((0:size(C{1}{1},2)-1))/size(C{1}{1},2) * (xlim(2) - xlim(1)) + xlim(1);
    y = ((0:size(C{1}{1},1)-1))/size(C{1}{1},1) * (ylim(2) - ylim(1)) + ylim(1);
    [X,Y]=meshgrid(x,y);
    xx = X(abs(C{1}{1}) > 0);
    yy = Y(abs(C{1}{1}) > 0);
    plot(xx,yy,color{1},'MarkerSize',marksize);
    for j=length(C)-1:-1:2,
        XX=[]; YY=[]; V=[]; Xsl=[]; Ysl=[];
        for l=1:length(C{j}),
            x = ((0:size(C{j}{l},2)-1))/size(C{j}{l},2) * (xlim(2) - xlim(1)) + xlim(1);
            y = ((0:size(C{j}{l},1)-1))/size(C{j}{l},1) * (ylim(2) - ylim(1)) + ylim(1);
            [X,Y]=meshgrid(x,y);
            xx = X(abs(C{j}{l}) > 0);
            yy = Y(abs(C{j}{l}) > 0);
            v = abs(C{j}{l}(abs(C{j}{l}) > 0));
            XX = [XX; xx(:)];
            YY = [YY; yy(:)];
            V = [V; v(:)];
            [xsl,ysl] = cslope(l,length(C{j}));
            Xsl = [Xsl; xsl*ones(prod(size(xx)),1)];
            Ysl = [Ysl; ysl*ones(prod(size(xx)),1)];
        end
        if ~isempty(XX),
            if arrows,
                if scalearrows,
                    quiver(XX, YY, Xsl.*V, Ysl.*V, autoscaling, color{j});
                else
                    quiver(XX, YY, 10*Xsl.*ones(size(V)), 10*Ysl.*ones(size(V)), 0*autoscaling, color{j});
                end
            else
                plot(XX(:),YY(:),color{j},'MarkerSize',marksize)
            end
        end
    end
    x = ((0:size(C{end}{1},2)-1))/size(C{end}{1},2) * (xlim(2) - xlim(1)) + xlim(1);
    y = ((0:size(C{end}{1},1)-1))/size(C{end}{1},1) * (ylim(2) - ylim(1)) + ylim(1);
    [X,Y]=meshgrid(x,y);
    xx = X(abs(C{end}{1}) > 0);
    yy = Y(abs(C{end}{1}) > 0);
    plot(xx,yy,color{end},'MarkerSize',marksize);
else
    for j=1:length(C),
        for l=1:length(C{j}),
            x = ((0:size(C{j}{l},1)-1) + 0.5)/size(C{j}{l},1);
            y = ((0:size(C{j}{l},2)-1) + 0.5)/size(C{j}{l},2);
            z = ((0:size(C{j}{l},3)-1) + 0.5)/size(C{j}{l},3);
            [X,Y,Z]=ndgrid(x,y,z);
            xx = X(abs(C{j}{l}) > 0);
            yy = Y(abs(C{j}{l}) > 0);
            zz = Z(abs(C{j}{l}) > 0);
            plot3(xx,yy,zz,color{j});
        end
    end
    view(3)
end


% This function gives slopes when image is displayed with y-axis upwards
% as with "imagesc(1:N,1:N,im), set(gca,'YDir','normal')", or with quiver
function [xsl,ysl] = cslope(idx,totnr)

idx = mod(idx-1,totnr/2);
if idx >= totnr/4,
  idx = idx - totnr/4;
  xsl = 1;
  ysl = (idx - (totnr/8) + 0.5)/(totnr/8);  
else
  ysl = 1;
  xsl = ((totnr/8) - idx - 0.5)/(totnr/8);
end
nrm=sqrt(xsl^2+ysl^2);
xsl=xsl/nrm;
ysl=ysl/nrm;

