function gdpl = get_planar_grid(gd)
%  GET_PLANAR_GRID  Returns three cross-sections of a full 3D vector field
%
%  Usage:
%    gdpl = get_planar_grid(gd)
%
%  Inputs:
%    gd		- output from spins_grid (or the grid structure from spins_gridparams)
%
%  Outputs:
%    gdpl	- structure containing the grid planes
%
%  David Deepwell, 2019

    ndims = length(fieldnames(gd));
    if ndims ~=3
        error('Grid must be 3 dimensional.')
    end

    yz.y = squeeze(gd.y(1,:,:));
    yz.z = squeeze(gd.z(1,:,:));
    xz.x = squeeze(gd.x(:,1,:));
    xz.z = squeeze(gd.z(:,1,:));
    xy.x = squeeze(gd.x(:,:,1));
    xy.y = squeeze(gd.y(:,:,1));
    gdpl.xy = xy;
    gdpl.xz = xz;
    gdpl.yz = yz;
end
