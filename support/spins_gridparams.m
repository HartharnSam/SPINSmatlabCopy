function gdpar = spins_gridparams(varargin)
% SPINS_GRIDPARAMS  Parse spins.conf, read in grid and deduce other grid parameters
%
%   gdpar = spins_gridparams()  gives the (vector) grid and parameters in a structure
%   gdpar = spins_gridparams('Full')  gives the full grid and parameters in a structure
%
%   David Deepwell, 2015
    global gdpar    

    % Parser spins.conf
    params = spins_params();

    % Read in grid
    if nargin > 1
        error('Too many inputs. spins_gridparams accepts either "Full" or "Vector".')
    elseif nargin == 1
        if ~strcmpi(varargin,'Full') && ~strcmpi(varargin,'Vector')
            error('Unacceptable grid type. Must be either "Full" or "Vector".')
        else
            gd = spins_grid(varargin{1});
        end
    elseif nargin == 0
        if strcmp(params.mapped_grid,'false')
            gd = spins_grid('Vector');
        %elseif strcmp(params.mapped_grid,'true') % default to give full grid if mapped grid
        %    gd = spins_grid('Full');
        end
    end

    % Add other information into params structure
    par = add_params(gd, params);

    % Place information into another structure
    gdpar.gd = gd;
    gdpar.params = par;
end

function par = add_params(gd, params)
    % Number of dimensions
    gdnames = fieldnames(gd);
    params.ndims = length(gdnames);

    % check if grid is vectorized
    if isvector(gd.(gdnames{1}))
       gdvec = gd;
    else  % get vectorized grid
        gdvec = get_vector_grid(gd, params);
    end
    try 
        x1d = gdvec.x;
    end
    try 
        y1d = gdvec.y;
    end
    try 
        z1d = gdvec.z;
    end

    % Number of points in each dimension
    if ~isfield(params,'Nx') && isfield(gd,'x')
        Nx = length(x1d);
        params.Nx = Nx;
    elseif isfield(params,'Nx') && isfield(gd,'x')
	Nx = params.Nx;		% save for use later
    end
    if ~isfield(params,'Ny') && isfield(gd,'y')
        Ny = length(y1d);
        params.Ny = Ny;
    elseif isfield(params,'Ny') && isfield(gd,'y')
	Ny = params.Ny;
    end
    if ~isfield(params,'Nz') && isfield(gd,'z')
        Nz = length(z1d);
	params.Nz = Nz;
    elseif isfield(params,'Nz') && isfield(gd,'z')
	Nz = params.Nz;
    end

    % Check vertical expansion type
    if ~isfield(params,'type_x') && isfield(gd,'x')
        if abs((x1d(Nz/2)-x1d(Nz/2-1))/(x1d(2)-x1d(1))) > 2
            params.type_x = 'NO_SLIP';
        else
            params.type_x = 'FREE_SLIP or PERIODIC';
        end
    end
    if ~isfield(params,'type_y') && isfield(gd,'y')
        if abs((y1d(Ny/2)-y1d(Ny/2-1))/(y1d(2)-y1d(1))) > 2
            params.type_y = 'NO_SLIP';
        else
            params.type_y = 'FREE_SLIP or PERIODIC';
        end
    end
    if ~isfield(params,'type_z') && isfield(gd,'z')
        if abs((z1d(Nz/2)-z1d(Nz/2-1))/(z1d(2)-z1d(1))) > 2
            params.type_z = 'NO_SLIP';
        else
            params.type_z = 'FREE_SLIP or PERIODIC';
        end
    end

    % Grid spacing for linear grids
    if ~isfield(params,'dx') && isfield(gd,'x')
        if ~strcmp(params.type_x,'NO_SLIP') && ~strcmp(params.type_x,'CHEBY')
            dx = x1d(2) - x1d(1);
            params.dx = dx;
        end    
    end 
    if ~isfield(params,'dy') && isfield(gd,'y')
        if ~strcmp(params.type_y,'NO_SLIP') && ~strcmp(params.type_y,'CHEBY') 
            dy = y1d(2) - y1d(1);
            params.dy = dy;
        end    
    end 
    if ~isfield(params,'dz') && isfield(gd,'z')
        if ~strcmp(params.type_z,'NO_SLIP') && ~strcmp(params.type_z,'CHEBY') 
            dz = z1d(2) - z1d(1);
            params.dz = dz;
        end
    end 

    % Length of domain (assumes topography only affects the z-coordinate)
    if ~isfield(params, 'Lx') && isfield(gd, 'x')
        if strcmp(params.type_x,'NO_SLIP') || strcmp(params.type_x,'CHEBY')
            params.Lx = roundn(abs(x1d(end)-x1d(1)), -15);
        else
            params.Lx = roundn(Nx*dx, -15);
        end
    end
    if ~isfield(params,'Ly') && isfield(gd,'y')
        if strcmp(params.type_y,'NO_SLIP') || strcmp(params.type_y,'CHEBY')
            params.Ly = roundn(abs(y1d(end)-y1d(1)), -15);
        else
            params.Ly = roundn(Ny*dy, -15);
        end
    end
    if ~isfield(params,'Lz') && isfield(gd,'z')
        if strcmp(params.type_z,'NO_SLIP') || strcmp(params.type_z,'CHEBY')
            params.Lz = roundn(max(gd.z(:)) - min(gd.z(:)), -15);
        else
            params.Lz = roundn(Nz*dz, -15);
        end
    end

    % Bounds of domain
    if ~isfield(params,'xlim') && isfield(gd,'x')
        if strcmp(params.type_x,'NO_SLIP') || strcmp(params.type_x,'CHEBY')
            x1 = min(x1d(1), x1d(end));
            x2 = max(x1d(1), x1d(end));
        else
            x1 = min(x1d(1), x1d(end)) - dx/2;
            x2 = max(x1d(1), x1d(end)) + dx/2;
        end
        xL = roundn(x1, -15);
        xR = roundn(x2, -15);
        params.xlim = [xL, xR];
    end
    if ~isfield(params,'ylim') && isfield(gd,'y')
        if strcmp(params.type_y,'NO_SLIP') || strcmp(params.type_y,'CHEBY')
            y1 = min(y1d(1), y1d(end));
            y2 = max(y1d(1), y1d(end));
        else
            y1 = min(y1d(1), y1d(end)) - dy/2;
            y2 = max(y1d(1), y1d(end)) + dy/2;
        end
        yL = roundn(y1, -15);
        yR = roundn(y2, -15);
        params.ylim = [yL, yR];
    end
    if ~isfield(params,'zlim') && isfield(gd,'z')
        if strcmp(params.type_z,'NO_SLIP') || strcmp(params.type_z,'CHEBY')
            z1 = min(gd.z(:));
            z2 = max(gd.z(:));
        else
            z1 = min(z1d(1), z1d(end)) - dz/2;
            z2 = max(z1d(1), z1d(end)) + dz/2;
        end
        zL = roundn(z1, -15);
        zR = roundn(z2, -15);
        params.zlim = [zL, zR];
    end

    par = params;
end
