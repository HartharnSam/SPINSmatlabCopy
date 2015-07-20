function gd = spins_grid(varargin)
% SPINS_GRID  loads the grid from SPINS.
%
%   gd = spins_grid()		gives vectors of each dimension
%   gd = spins_grid('Full')	gives full matrices
%
% David Deepwell, 2015.

if nargin > 1
    error('Too many input arguments. spins_grid accepts "Vector" or "Full".');
elseif nargin == 1 && ~strcmpi(varargin,'Vector') && ~strcmpi(varargin,'Full')
    error('Argument not understood. spins_grid accepts "Vector" or "Full".');
elseif nargin == 0 || strcmpi(varargin,'Vector')
    if exist('xgrid_reader.m', 'file') == 2
        try
            gd.x = xgrid_reader(':',1,1);
        catch
            gd.x = xgrid_reader(':',1);
        end
    end
    if exist('ygrid_reader.m', 'file') == 2
        try
            gd.y = ygrid_reader(1,':',1);
        catch
            gd.y = ygrid_reader(1,':');	% assumes 2D is never y-z plane
        end
    end
    if exist('zgrid_reader.m', 'file') == 2
        try
            gd.z = zgrid_reader(1,1,':');
        catch
            gd.z = zgrid_reader(1,':');
        end
    end
elseif strcmpi(varargin,'Full') && nargin == 1
    if exist('xgrid_reader.m', 'file') == 2
        gd.x = xgrid_reader();
    end
    if exist('ygrid_reader.m', 'file') == 2
        gd.y = ygrid_reader();
    end
    if exist('zgrid_reader.m', 'file') == 2
        gd.z = zgrid_reader();
    end
end
