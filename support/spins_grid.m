function gd = spins_grid(varargin)
% SPINS_GRID  load the grid from SPINS.
%
%   gd = spins_grid()		gives the grid in a structure
%
%   Optional arguments are:
%	'Vector'    - (default) gives vector grid in output structure
%	'Full'      -  gives the entire grid in output structure
%
% David Deepwell, 2015.

if nargin > 1
    error('Too many input arguments. spins_grid accepts "Vector" or "Full".');
elseif nargin == 1 && ~strcmpi(varargin,'Vector') && ~strcmpi(varargin,'Full')
    error('Argument not understood. spins_grid accepts "Vector" or "Full".');
elseif nargin == 0 || strcmpi(varargin,'Vector')
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
