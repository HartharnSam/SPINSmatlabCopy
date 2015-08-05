function [xgrid, ygrid, fixedzdata] = get_fixed_z(x,y,z,data,depth)
% GET_FIXED_Z  read in 3D data and grid, return the 2D data at a constant depth
%
%   used in spins_plot2dmapped
%
%   David Deepwell, 2015

    % create empty arrays
    sz = size(x);
    fixedzdata = zeros(sz(1),sz(2));
    xgrid = fixedzdata;
    ygrid = xgrid;

    % interpolate between neighbours at selected depth
    for xi = 1:sz(1)
        xgrid(xi,:) = squeeze(x(xi,:,1));
        ygrid(xi,:) = squeeze(y(xi,:,1));
        fixedzdata(xi,:) = interp1(squeeze(z(xi,1,:)), squeeze(data(xi,:,:))',...
                                   depth,'spline',NaN);
    end
end

