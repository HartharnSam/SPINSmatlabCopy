function gdvec = get_vector_grid(gd, params)
% GET_VECTOR_GRID  take a full vector field and give the vectors of each dimension
%
%   gdvec = get_vector_grid(gd, params) takes grid structure, gd, and parameters structure
%
%   David Deepwell, 2015
    if params.ndims == 3
        gdvec.x = gd.x(:,1,1);
        gdvec.y = gd.y(1,:,1);
        gdvec.z = gd.z(1,1,:);
    else
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
    end
end
