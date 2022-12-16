% it loads the grid for Chris' model (assuming 3D)
function [x, y, z, x1d, z1d, sz, NX, NZ, dx, Lz, Lx, wtz, wt] = spinsgrid3d()
x = xgrid_reader();
y = zgrid_reader();
z = zgrid_reader();
slice_grids;
% z1d=squeeze(z(1,:));
% x1d=squeeze(x(:,1));
% sz=size(x);
% NX=sz(1);
% NZ=sz(2);
% 
% sz=size(x);
% dx=x(2,2)-x(1,1);
% Lz=max(abs(z(:)));
% Lx=sz(1)*dx;
% sz=size(x);
% [wtz,wt]=clencurt(sz(2)-1); wt=wt*Lz/2;
