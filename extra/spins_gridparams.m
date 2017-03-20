function gdpar = spins_gridparams(varargin)
%  SPINS_GRIDPARAMS  Parse spins.conf, read in grid and deduce other grid parameters
%
%  Usage:
%    gdpar = spins_gridparams()  gives the grid and parameters in a structure
%
%  Inputs:
%
%    Optional arguments are:
%   First argument:
%	'Vector'    - (default) gives vector grid in output structure
%	'Full'      - gives the entire grid in output structure
%	'FastFull'  - creates the entire grid from parameters in spins.conf (equivalent to 'Full')
%   Second argument:
%   {boolean}   - check if grid is mapped or not?
%
%  Outputs:
%    gdpar	- a structure containing two structures: the grid (gd), and parameters (params)
%
%  David Deepwell, 2015
global gdpar    

    % Parse spins.conf
    params = spins_params();

    % Read in grid
    if nargin > 2
        error('Too many inputs. spins_gridparams accepts either "Full", "FastFull", or "Vector".')
    elseif nargin == 2
        if islogical(varargin{2})
            check_grid = varargin{2};
        else
            error('The second argument was not a boolean. Use false or true.')
        end
        gd = spins_grid(varargin{1});
    elseif nargin == 1
        check_grid = true;
        gd = spins_grid(varargin{1});
    elseif nargin == 0
        check_grid = true;
        gd = spins_grid('Vector');
    end

    % Add other information into params structure
    par = add_params(gd, params, check_grid, varargin);

    % Place information into output structure
    gdpar.gd = gd;
    gdpar.params = par;
end

function par = add_params(gd, params, check_grid, varargin)
    % parse varargin
    if length(varargin{1}) == 0
       vectorized = true;
    else
        if strcmpi(varargin{1},'Vector')
            vectorized = true;
        else
            vectorized = false;
        end
    end

    % get list of grid names
    gdnames = fieldnames(gd);

    % check if grid is vectorized
    if isvector(gd.(gdnames{1}))
       gdvec = gd;
    else  % get vectorized grid
        gdvec = get_vector_grid(gd);
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

    % Number of dimensions
    gdnames = fieldnames(gd);
    params.ndims = length(gdnames);
    % print error if one dimensions is accidentally written out
    if params.ndims == 3
        if Nx == 1 || Ny == 1 || Nz == 1
            warning('Grid appears to be 2 dimensional but all 3 grid fields are present. Remove the unnecessary field before proceeding.')
        end
    end

    % add gravitational constant    
    if ~isfield(params, 'g')
        params.g = 9.81;
    end

    % add reference density
    if ~isfield(params, 'rho_0')
        params.rho_0 = 1;
        warning('Reference density is not specified, default is chosen to be 1.')
    end

    % add number of outputs
    if ~isfield(params, 'noutputs')
        noutputs = length(dir('u.*'));
        if exist('u.dump','file')
            noutputs = noutputs - 1;
        end
        params.noutputs = noutputs;
    end

    % check if grid is mapped
    if isfield(gd,'z') && check_grid
        % get vector of depths at mid-depth of domain
        if isvector(gd.z)
            middepth = spins_reader('zgrid',1:10:Nx,1,round(Nz/2));
        else
            if params.ndims == 3
                middepth = gd.z(1:10:Nx,1,round(Nz/2));
            elseif params.ndims == 2
                middepth = gd.z(1:10:end,round(Nz/2));
            end
        end
        % grid is mapped if there is variation in the depth
        zratio =  min(middepth)/max(middepth);
        if zratio ~= 1
            if vectorized
                error('The grid appears to be mapped. Use the "Full" option rather than "Vector".')
            end
            if isfield(params,'mapped_grid') == false
                params.mapped_grid = 'true';
            elseif strcmp(params.mapped_grid, 'false')
                error('The grid appears to be mapped, but the config file says otherwise. Fix before proceeding.')
            end
        else
            if isfield(params,'mapped_grid') == false
                params.mapped_grid = 'false';
            elseif strcmp(params.mapped_grid, 'true')
                error('The grid appears to be unmapped, but the config file says otherwise. Fix before proceeding.')
            end
        end
    else
        % if not checking grid, assume that grid is unmapped if field doesn't exist
        if isfield(params,'mapped_grid') == false
            params.mapped_grid = 'false';
            disp('mapped_grid not found in spins.conf. Assuming it is false.')
        end
    end

    % Check vertical expansion type
    if ~isfield(params,'type_x') && isfield(gd,'x')
        midx = round(Nx/2);
        if abs((x1d(midx) - x1d(midx-1))/(x1d(2) - x1d(1))) > 2
            params.type_x = 'NO_SLIP';
        else
            params.type_x = 'FREE_SLIP or PERIODIC';
        end
    end
    if ~isfield(params,'type_y') && isfield(gd,'y')
        midy = round(Ny/2);
        if abs((y1d(midy) - y1d(midy-1))/(y1d(2) - y1d(1))) > 2
            params.type_y = 'NO_SLIP';
        else
            params.type_y = 'FREE_SLIP or PERIODIC';
        end
    end
    if ~isfield(params,'type_z') && isfield(gd,'z')
        midz = round(Nz/2);
        if abs((z1d(midz) - z1d(midz-1))/(z1d(2) - z1d(1))) > 2
            params.type_z = 'NO_SLIP';
        else
            params.type_z = 'FREE_SLIP or PERIODIC';
        end
    end

    % Grid spacing for linear grids
    if isfield(gd, 'x')
        if isfield(params, 'dx')
            warning(['Use of dx is reserved for grid spacing. ',...
                     'dx from spins.conf has been renamed as dx_conf.'])
            params.dx_conf = params.dx;
        end
        if ~strcmp(params.type_x,'NO_SLIP') && ~strcmp(params.type_x,'CHEBY')
            dx = x1d(2) - x1d(1);
            params.dx = dx;
        else
            params.dx = '-';
        end    
    end
    if isfield(gd, 'y')
        if isfield(params, 'dy')
            warning(['Use of dy is reserved for grid spacing. ',...
                     'dy from spins.conf has been renamed as dy_conf.'])
            params.dy_conf = params.dy;
        end
        if ~strcmp(params.type_y,'NO_SLIP') && ~strcmp(params.type_y,'CHEBY') 
            dy = y1d(2) - y1d(1);
            params.dy = dy;
        else
            params.dy = '-';
        end    
    end 
    if isfield(gd, 'z')
        if isfield(params, 'dz')
            warning(['Use of dz is reserved for grid spacing. ',...
                     'dz from spins.conf has been renamed as dz_conf.'])
            params.dz_conf = params.dz;
        end
        if ~strcmp(params.type_z,'NO_SLIP') && ~strcmp(params.type_z,'CHEBY') 
            dz = z1d(2) - z1d(1);
            params.dz = dz;
        else
            params.dz = '-';
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
