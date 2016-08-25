function ar = spins_reader(varname, seq, xrange, yrange, zrange)
% SPINS data reader with slab support and a user defined name
% Opens and reads a slab of SPINS data, optionally
% loading only a portion of the total array. This
% functionality is most useful when only a portion
%  of a large, 3D array is needed.
%
% Usage:
%   slab = spins_reader(varname,seq,xrange,yrange,zrange);
% Where varname is the variable to read,
% seq is the sequence number of the output, and
% xrange, yrange, and zrange are the ranges
% of the array to be read.
% Empty values, [], and : each imply reading the
% full dimension.

% Version 1.1, July 09 2012.  Original general
% MATLAB code provided courtesy of Michael Dunphy
% (mdunphy@uwaterloo.ca), adapted for SPINS by 
% Christopher Subich (csubich@uwaterloo.ca).
params = spins_params();
Nx = params.Nx;
Ny = params.Ny;
Nz = params.Nz;

% Sanitize the ranges:
if (~exist('xrange') || isempty(xrange) || strcmp(xrange,':')) xrange = [1:Nx]; end;
if (~exist('yrange') || isempty(yrange) || strcmp(yrange,':')) yrange = [1:Ny]; end;
if (~exist('zrange') || isempty(zrange) || strcmp(zrange,':')) zrange = [1:Nz]; end;
xrange(xrange < 1) = []; xrange(xrange > Nx) = [];
yrange(yrange < 1) = []; yrange(yrange > Ny) = [];
zrange(zrange < 1) = []; zrange(zrange > Nz) = [];

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
