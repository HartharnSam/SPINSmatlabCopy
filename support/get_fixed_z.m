function [xgrid, ygrid, fixedzfield] = get_fixed_z(x, y, z, field, depth)
%  GET_FIXED_Z  read in 3D field and grid, return the 2D field at a constant depth.
%
%  Usage:
%	[xgrid, ygrid, fixedzfield] = get_fixed_z(x, y, z, field, depth)
%
%  Inputs:
%    x, y, z, and 'field'	- three dimensional fields
%    'depth'			- a real number
%
%  Outputs:
%    xgrid, ygrid	- horizontal cross-sections of the grid
%    fixedzfield	- interpolated value of 'field' at 'depth'
%
%  David Deepwell, 2015
%  adapted from code by Marek Stastna 

    % create empty arrays
    sz = size(x);
    fixedzfield = zeros(sz(1),sz(2));
    xgrid = fixedzfield;
    ygrid = xgrid;

    % interpolate between neighbours at selected depth
    for xi = 1:sz(1)
        xgrid(xi,:) = squeeze(x(xi,:,1));
        ygrid(xi,:) = squeeze(y(xi,:,1));
        fixedzfield(xi,:) = interp1(squeeze(z(xi,1,:)), squeeze(field(xi,:,:))',...
                                   depth,'spline',NaN);
    end
end

