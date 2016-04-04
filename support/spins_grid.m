function gd = spins_grid(varargin)
%  SPINS_GRID  load the grid from SPINS.
%
%  Usage:
%    gd = spins_grid() gives the vector grid in a structure
%
%  Inputs:
%    Optional arguments:
%	'Vector'    - (default) gives vector grid in output structure
%	'Full'      - gives the entire grid in output structure
%
%  Outputs:
%    gd		- structure containing the grid
%
%  David Deepwell, 2015

if nargin > 1
    error('Too many input arguments. spins_grid accepts "Vector" or "Full".');
elseif nargin == 1 && ~strcmpi(varargin,'Vector') && ~strcmpi(varargin,'Full')
    error('Argument not understood. spins_grid accepts "Vector" or "Full".');
% make vector grid
elseif nargin == 0 || strcmpi(varargin,'Vector')
    % check if some parameters are set in spins.conf
    params = spins_params();
    eNx = isfield(params,'Nx');
    eNy = isfield(params,'Ny');
    eNz = isfield(params,'Nz');
    eLx = isfield(params,'Lx');
    eLy = isfield(params,'Ly');
    eLz = isfield(params,'Lz');
    emin_x = isfield(params,'min_x');
    emin_y = isfield(params,'min_y');
    emin_z = isfield(params,'min_z');
    etype_x = isfield(params,'type_x');
    etype_y = isfield(params,'type_y');
    etype_z = isfield(params,'type_z');

    % are all fields present?
    e3D = eNx*eNy*eNz*...
          eLx*eLy*eLz*...
          emin_x*emin_y*emin_z*...
          etype_x*etype_y*etype_z;
    e2D = eNx*eNz*...
          eLx*eLz*...
          emin_x*emin_z*...
          etype_x*etype_z;
    % find dimensionality
    if eNx*eNy*eNz == 1
        Nx = params.Nx;
        Ny = params.Ny;
        Nz = params.Nz;
        if Nx == 1 || Ny == 1 || Nz == 1
            e3D = 0;
        end
    end

    % create or write grid
    if e3D == 1 || e2D == 1
        Lx = params.Lx;
        Lz = params.Lz;
        min_x = params.min_x;
        min_z = params.min_z;
        type_x = params.type_x;
        type_z = params.type_z;
        if e3D == 1
            Ly = params.Ly;
            min_y = params.min_y;
            type_y = params.type_y;
        end

        if strcmpi(type_x,'FREE_SLIP')
            gd.x = min_x + Lx*([0:(Nx-1)]+0.5)/Nx;
        else
            gd.x = min_x + Lx*0.5*(1 - cos(pi*[0:(Nx-1)]/(Nx-1)));
        end
        if strcmpi(type_z,'FREE_SLIP')
            gd.z = min_z + Lz*([0:(Nz-1)]+0.5)/Nz;
        else
            gd.z = min_z + Lz*0.5*(1 - cos(pi*[0:(Nz-1)]/(Nz-1)));
        end
        if e3D == 1
            if strcmpi(type_y,'FREE_SLIP')
                gd.y = min_y + Ly*([0:(Ny-1)]+0.5)/Ny;
            else
                gd.y = min_y + Ly*0.5*(1 - cos(pi*[0:(Ny-1)]/(Ny-1)));
            end
        end
    else
        if (exist('xgrid_reader.m', 'file') == 2) && (exist('xgrid', 'file') == 2)
            try
                gd.x = xgrid_reader(':',1,1);
            catch
                gd.x = xgrid_reader(':',1);
            end
        end
        if (exist('ygrid_reader.m', 'file') == 2) && (exist('ygrid', 'file') == 2)
            try
                gd.y = ygrid_reader(1,':',1);
            catch
                gd.y = ygrid_reader(1,':');	% assumes 2D is never y-z plane
            end
        end
        if (exist('zgrid_reader.m', 'file') == 2) && (exist('zgrid', 'file') == 2)
            try
                gd.z = zgrid_reader(1,1,':');
            catch
                gd.z = zgrid_reader(1,':');
            end
        end
    end
% make Full grid
elseif strcmpi(varargin,'Full') && nargin == 1
    if (exist('xgrid_reader.m', 'file') == 2) && (exist('xgrid', 'file') == 2)
        gd.x = xgrid_reader();
    end
    if (exist('ygrid_reader.m', 'file') == 2) && (exist('ygrid', 'file') == 2)
        gd.y = ygrid_reader();
    end
    if (exist('zgrid_reader.m', 'file') == 2) && (exist('zgrid', 'file') == 2)
        gd.z = zgrid_reader();
    end
end

% error message if no grid_readers are in directory
if ~exist('gd','var')
    error('The grid was not found in the working directory. Do you know where you are?')
end
