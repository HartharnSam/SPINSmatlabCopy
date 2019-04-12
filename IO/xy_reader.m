function ar = xy_reader(varname, seq, varargin);
% SPINS xy reader with slab support and a user defined name
% Opens and reads a slab of SPINS data, optionally
% loading only a portion of the total array. This
% functionality is most useful when only a portion
% of a large, 3D array is needed.
%
% Usage:
%   Entire field:
%     slab = xy_reader(varname, seq);
%   2D field (portion):
%     slab = xy_reader(varname, seq, xrange, zrange);
%
% Where:
%   - varname is the file to read
%   - seq is the sequence number of the output
%   - xrange and zrange are the ranges of the
%     array to be read. Empty values, [], and : each imply
%     reading the full dimension. 2D data are expected
%     to be the x-z plane.
%     1D and 2D fields can also be read using the 3D format,
%     so long as the singleton dimensions have *range=1.
%     ie. slab = xy_reader(varname, 1, 1, zrange);

% Apr 2019.  Original general MATLAB code provided
% courtesy of Michael Dunphy (mdunphy@uwaterloo.ca), adapted
% for SPINS by Christopher Subich (csubich@uwaterloo.ca)
% and David Deewell (ddeepwel@uwaterloo.ca).

%% Find dimensions
params = spins_params();
Nx = params.Nx;
Ny = params.Ny;
Nz = 1;

%% Get ranges
if nargin==4
    % if 2D ranges given and varname is a field
    xrange = varargin{1};
    yrange = varargin{2};
end

% Sanitize the ranges:
if (~exist('xrange') || isempty(xrange) || strcmp(xrange,':')) xrange = [1:Nx]; end;
if (~exist('yrange') || isempty(yrange) || strcmp(yrange,':')) yrange = [1:Ny]; end;
zrange = [1];
xrange(xrange < 1) = []; xrange(xrange > Nx) = [];
yrange(yrange < 1) = []; yrange(yrange > Ny) = [];
% reset ranges if they were input as values outside of the extrema
if isempty(xrange) xrange = [1:Nx]; end;
if isempty(yrange) yrange = [1:Ny]; end;

% Define ranges in file-ordering
ranges = {xrange,yrange,zrange};
ranges = ranges([2,3,1]);
% Output file name:
fname = sprintf('%s.%d',varname,seq);

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
