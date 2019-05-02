function F=eggshape(dim, N, x0, x1, R, r)

% F = eggshape(dim, N, x0, x1, R, r)
%
% Produces the characteristic function of an 'egg'
% in dim (=2 or 3) dimensions with N^dim grid points
% in the domain [0,1)^dim
% the 'egg' is the convex hull of the union of two discs
% with radii R and r and centered at x0 and x1 respectively
%

if r > R,
    tmp = r; r = R; R = tmp;
    tmp = x0; x0 = x1; x1 = tmp;
end

x=0:1./N:1-0.5./N;
if dim==2,
  [X,Y] = meshgrid(x,x);
 
  % circular parts
  F = ((X - x0(1)).^2 + (Y - x0(2)).^2 <= R^2);
  F = F | ((X - x1(1)).^2 + (Y - x1(2)).^2 <= r^2);
 
  % straight parts
  dx = (x1-x0)/norm(x1-x0);
  dxo = [-dx(2); dx(1)];
  edx = (R-r)/norm(x1-x0);
  edxo = sqrt(1-edx^2);
  e{1} = [dx dxo] * [edx; edxo];
  e{2} = [dx -dxo] * [edx; edxo];  
  eo{1} = [-e{1}(2); e{1}(1)];
  eo{2} = [e{2}(2); -e{2}(1)];
  G = ones(size(F));
  for k=1:2,   
      pe0 = ((X - x0(1)) * e{k}(1) + (Y - x0(2)) * e{k}(2));
      pe0o = ((X - x0(1)) * eo{k}(1) + (Y - x0(2)) * eo{k}(2));
      pe1 = ((X - x1(1)) * e{k}(1) + (Y - x1(2)) * e{k}(2));
      pe1o = ((X - x1(1)) * eo{k}(1) + (Y - x1(2)) * eo{k}(2));
      F = F | (pe0 <= R & pe1 >= -r & pe1o >= 0 & pe0o <=0);
      G = G & (pe1 <= r & pe0 >= 0 & pe1o >= 0 & pe0o <= 0);
  end
  F = F | G;
else
    
    
end

F = double(F);

