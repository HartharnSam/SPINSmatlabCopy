function ar = spins_reader(varname, varargin);
% SPINS data reader with slab support and a user defined name
% Opens and reads a slab of SPINS data, optionally
% loading only a portion of the total array. This
% functionality is most useful when only a portion
% of a large, 3D array is needed.
%
% Usage:
%   Entire field:
%     slab = spins_reader(varname, seq);
%   3D field (portion):
%     slab = spins_reader(varname, seq, xrange, yrange, zrange);
%   2D field (portion):
%     slab = spins_reader(varname, seq, xrange, zrange);
%
%   Entire grid:
%     slab = spins_reader(varname);
%   3D grid (portion):
%     slab = spins_reader(varname, xrange, yrange, zrange);
%   2D grid (portion):
%     slab = spins_reader(varname, xrange, zrange);
% Where:
%   - varname is the variable to read (rho, u, v, xgrid,...)
%   - seq is the sequence number of the output
%   - xrange, yrange, and zrange are the ranges of the
%     array to be read. Empty values, [], and : each imply
%     reading the full dimension. 2D data are expected
%     to be the x-z plane.
%     1D and 2D fields can also be read using the 3D format,
%     so long as the singleton dimensions have *range=1.
%     ie. slab = spins_reader(varname, 1, 1, zrange);

% Version 2.0, Jan 2017.  Original general MATLAB code provided
% courtesy of Michael Dunphy (mdunphy@uwaterloo.ca), adapted
% for SPINS by Christopher Subich (csubich@uwaterloo.ca)
% and David Deewell (ddeepwel@uwaterloo.ca).

%% Find dimensions
params = spins_params();
Nx = params.Nx;
Nz = params.Nz;
if ~isfield(params, 'Ny')
    Ny = 1;
else
    Ny = params.Ny;
end

% a grid doesn't have an extension
var_list = dir([varname,'.*']);
if length(var_list) == 0
    is_grid = true;
else
    is_grid = false;
end

%% Get sequence if reading a field variable
if ~is_grid
    if nargin >= 2
        seq = varargin{1};
    else
        error('Output number not understood or given.')
    end
end

%% Get ranges
if is_grid && nargin==4
    % if all ranges given and varname is a grid
    xrange = varargin{1};
    yrange = varargin{2};
    zrange = varargin{3};
elseif ~is_grid && nargin==5
    % if all ranges given and varname is a field
    xrange = varargin{2};
    yrange = varargin{3};
    zrange = varargin{4};
elseif is_grid && nargin==3 && Ny==1
    % if 2D ranges given and varname is a grid
    xrange = varargin{1};
    zrange = varargin{2};
elseif ~is_grid && nargin==4 && Ny==1
    % if 2D ranges given and varname is a field
    xrange = varargin{2};
    zrange = varargin{3};
elseif (is_grid && nargin>1) || (~is_grid && nargin>2)
    fname = sprintf('%s.%d',varname,varargin{1});
    if ~(exist(fname, 'file') == 2)
        error([fname,' does not exist.']);
    else
        error('Input arguments not understood.')
    end
end

% Sanitize the ranges:
if (~exist('xrange') || isempty(xrange) || strcmp(xrange,':')) xrange = [1:Nx]; end;
if (~exist('yrange') || isempty(yrange) || strcmp(yrange,':')) yrange = [1:Ny]; end;
if (~exist('zrange') || isempty(zrange) || strcmp(zrange,':')) zrange = [1:Nz]; end;
xrange(xrange < 1) = []; xrange(xrange > Nx) = [];
yrange(yrange < 1) = []; yrange(yrange > Ny) = [];
zrange(zrange < 1) = []; zrange(zrange > Nz) = [];
% reset ranges if they were input as values outside of the extrema
if isempty(xrange) xrange = [1:Nx]; end;
if isempty(yrange) yrange = [1:Ny]; end;
if isempty(zrange) zrange = [1:Nz]; end;

% Define ranges in file-ordering
ranges = {xrange,yrange,zrange};
ranges = ranges([2,3,1]);
% Output file name:
if is_grid
    fname = sprintf('%s',varname);
else
    fname = sprintf('%s.%d',varname,seq);
end

% memmap the file for reading
m = memmapfile(fname, 'Offset',0, ...
   'Format', {'double' [Ny,Nz,Nx] 'x'}, ...
   'Writable',false);

% Extract the data and clear the memmap
ar = m.Data.x(ranges{1},ranges{2},ranges{3}); clear m

% Permute, check endianness, and return
ar = squeeze(ipermute(ar,[2,3,1]));
[~,~,endian] = computer();
if (~isequal(endian,'L'))
   ar = swapbytes(ar);
end
