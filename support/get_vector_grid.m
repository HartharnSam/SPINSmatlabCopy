function gdvec = get_vector_grid(gd)
%  GET_VECTOR_GRID  take a full vector field and give the vectors along each dimension
%
%  Usage:
%    gdvec = get_vector_grid(gd) takes grid structure, gd, and parameters structure
%
%  Inputs:
%    gd		- output from spins_grid (or the grid structure from spins_gridparams)
%
%  Outputs:
%    gdvec	- structure containing the grid vectors
%
%  David Deepwell, 2015

    ndims = length(fieldnames(gdpar.gd));
    if ndims == 3
        gdvec.x = gd.x(:,1,1);
        gdvec.y = gd.y(1,:,1);
        gdvec.z = gd.z(1,1,:);
    elseif ndims == 2
        if isfield(gd, 'x')
            gdvec.x = gd.x(:,1);
            try
                gdvec.y = gd.y(1,:);
            catch
                gdvec.z = gd.z(1,:);
            end
        else
            gdvec.y = gd.y(:,1);
            gdvec.z = gd.z(1,:);
        end
    else
        error('Grid must be either 2 or 3 dimensional.')
    end
end
