function max_dxyz = max_grid_spacing(gdpar)
%  MAX_GRID_SPACING Return the maximum grid spacing in each dimension
%
%  Usage:
%    max_dxyz = max_grid_spacing()
%
%  Inputs:
%    n/a
%
%  Outputs:
%    'max_dxyz'	- a vector of the max grid spacing in each dimension
%
%  David Deepwell, 2018

% separate grid and parameters
split_gdpar
% vectorize grid
gdvec = get_vector_grid(gd);

% max x-grid spacing
dxs = diff(gdvec.x);
dx_max = max(dxs);

% max y-grid spacing
if params.ndims == 3
    dys = diff(gdvec.y);
    dy_max = max(dys);
else
    dy_max = 0;
end

% max z-grid spacing
if strcmp(params.mapped_grid,'false')
    % if mapped
    if params.ndims == 3
        dzs = diff(gd.z, 1, 3);
    else
        dzs = diff(gd.z, 1, 2);
    end
    dz_max = max(dzs(:));
else
    % if unmapped
    dzs = diff(gdvec.z);
    dz_max = max(dzs);
end

max_dxyz = [dx_max dy_max dz_max];
