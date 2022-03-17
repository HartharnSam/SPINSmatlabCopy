function [h xout yout] = spins_contour(x,y,z,n);
% my_nice_contour -- creates a contour plot using line objects
%
% syntax: [handle x_out y_out] = my_nice_contour(x,y,z,n)
% x,y -- input two-dimensional grid, of the syntax used by 
%        the MATLAB function contours
% z   -- input function to be contoured (nominally the z-values)
% n   -- number (if scalar) or levels (if vector) to contour
% outputs:
% h   -- LINE HANDLES of output lines (drawn on current figures)
%        By returning handles, linestyles and thicknesses can be
%        changed after the fact with get/set(h,'property','value')
% x_out, y_out
%     -- Cell array (of the same length as h) of the values plotted.
%        Use this if further processing of the contour levels is
%        necessary.


xx = x;
yy = y;

% The matlab function "contours" ends up doing the number-crunching
% work for us to find the contours.  Its return value is a matrix
% consisting of countour levels and the (x,y) points along those
% levels, which we can then feed to "line" for plotting.
cmat = contours(xx,yy,z,n);

ind = 1;
h = [];
xout = {};
yout = {};
while (ind < size(cmat,2))
   len = cmat(2,ind);
   xs = cmat(1,ind+(1:len));
   xout{ind} = xs;
   ys = cmat(2,ind+(1:len));
   yout{ind} = ys;
   h(end+1) = line(xs,ys);
   ind = ind + len+1;
end
