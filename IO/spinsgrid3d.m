function [x, y, z, xv, zv, xh, yh, yp, zp, sz, NX, NY, NZ, dx, Lz, Lx] = spinsgrid3d()
x = xgrid_reader();
y = ygrid_reader();
z = zgrid_reader();

if nargout > 3
    sz = size(x);
    NX = sz(1);
    NY = sz(2);
    NZ = sz(3);
    
    xv=xzslice(x,1);
    zv=xzslice(z,1);
    xh=xyslice(x,NZ);
    yh=xyslice(y,NZ);
    if nargout > 7
        yp=yzslice(y,1);
        zp=yzslice(z,1);
    end
end
% z1d=squeeze(z(1,:));
% x1d=squeeze(x(:,1));
%
% sz=size(x);
% dx=x(2,2)-x(1,1);
% Lz=max(abs(z(:)));
% Lx=sz(1)*dx;
% sz=size(x);
% [wtz,wt]=clencurt(sz(2)-1); wt=wt*Lz/2;