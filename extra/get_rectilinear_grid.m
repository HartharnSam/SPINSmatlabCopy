function gd_rect = get_rectilinear_grid(gd)
%  GET_RECTILINEAR_GRID  Take a mapped grid and interpolate onto
%                        a rectilinear one
%
%  Usage:
%    gd_rect = get_rectilinear_grid(gd);
%
%  Inputs:
%    gd      - Grid structure (from spins_grid or spins_gridparams)
%
%  Outputs:
%    gd_rect - Structure containing the rectilinear grid
%
%  David Deepwell, 2019

% grid resolution
Ndims = length(size(gd.x));
if Ndims == 2
    [Nx Nz] = size(gd.x);
elseif Ndims == 3
    [Nx Ny Nz] = size(gd.x);
end

% grid limits
min_x = min(gd.x(:));
max_x = max(gd.x(:));
min_z = min(gd.z(:));
max_z = max(gd.z(:));

% create rectilinear grid
x1d = linspace(min_x, max_x, Nx);
z1d = linspace(min_z, max_z, Nz);
[z, x] = meshgrid(z1d, x1d);

% place new grids into structure
gd_rect.x = x;
gd_rect.z = z;

% Adjust for 3D
if Ndims == 3
    % extend grids into the span-wise dimension
    x3d = repmat(reshape(x, [Nx 1 Nz]), [1 Ny 1]);
    z3d = repmat(reshape(z, [Nx 1 Nz]), [1 Ny 1]);

    % place new grids into structure
    gd_rect.y = gd.y;
    gd_rect.x = x3d;
    gd_rect.z = z3d;
end
