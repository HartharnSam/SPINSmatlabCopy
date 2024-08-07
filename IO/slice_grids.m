function [xv, zv, xh, yh, yp, zp] = slice_grids(x, y, z)
% This M-file creates the various 2D grids for 2d slices in SPINS
xv=xzslice(x,1);
zv=xzslice(z,1);
if nargin > 2
    xh=xyslice(x,1);
    yh=xyslice(y,1);
    yp=yzslice(y,1);
    zp=yzslice(z,1);
end
