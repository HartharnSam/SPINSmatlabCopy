function field = resize_z(field_name, field, Nz_new, varargin)
% RESIZE_Z  Increase or decrease the resolution in the z-dimension.
%
%  Assumptions:
%    - if 2D, data must be in x-z plane
%    - Boundary condition must be free-slip or periodic
%      - To keep the code simple the periodic BC is not optimized
%      - since the field will still be doubled
%
%  Usage:
%    rho_new = resize_z('rho', rho, 512);
%
%  Inputs:
%    'field_name'   - the name of the field to resize
%    'field'        - the field data
%    'Nz_new'       - number of points to change to
%   Optional arguments:
%    'opts'         - a structure holding simulation parameters
%                   - these are: type_z and Lz
%
%  Outputs:
%    'field'        - the resized field
%
%  David Deepwell, 2018

% get number of points in each dimension
sz     = size(field);
Nx     = sz(1);
Nz_old = sz(end);
if length(sz) == 3
    Ny = sz(2);
else
    Ny = 1;
end
% get other parameters
if nargin == 3
    % read from spins.conf
    params = spins_params();
    type_z = params.type_z;
    Lz     = params.Lz;
elseif nargin == 4
    % read from optional input argument
    opts = varargin{1};
    type_z = opts.type_z;
    Lz     = opts.Lz;
else
    error('Incorrect number of input arguments. 3 or 4 are expected.')
end
dz_old = Lz/Nz_old;

% error check expansion type
if ~strcmp(type_z, 'FREE_SLIP') && ~strcmp(type_z, 'FOURIER')
    error('z boundary condition is not FREE_SLIP or FOURIER.')
end

% check sizes
mult = Nz_new/Nz_old;
if mult > 1
    increase_points = true;
else
    increase_points = false;
end

% reshape 2D to match the shape of a 3D field
if Ny == 1
    field = reshape(field, [Nx 1 Nz_old]);
end
% permute to put extending dimension in the first dimension
field = permute(field, [3 2 1]);

% extension types
odd_ext    = @(f) [f; -flipud(f)];
even_ext   = @(f) [f;  flipud(f)];
period_ext = @(f) [f;  f];

% double the field
if strcmp(type_z, 'FOURIER')
    field = period_ext(field);
elseif strcmp(field_name, 'w')
    field = odd_ext(field);
else
    field = even_ext(field);
end

% phase shift because grid is cell-centered, and take fft
dz_new = Lz/Nz_new;
kzs = fftfreq(2*Nz_old, dz_old);
phase_corr = exp(-1i*kzs*(dz_old/2-dz_new/2)).';
field = bsxfun(@times, fft(field), phase_corr);

% pad or truncate
if increase_points
    % pad with zeros
    field = [field(1:Nz_old,:,:); zeros(2*Nz_new-2*Nz_old,Ny,Nx); field(end-Nz_old+1:end,:,:)];
else
    % truncate high frequencies
    field = [field(1:Nz_new,:,:); field(end-Nz_new+1:end,:,:)];
end

% take ifft
field = mult*real(ifft(field));
field = field(1:Nz_new,:,:);
% permute back
field = permute(field, [3 2 1]);
% return field to 2D if it was
if Ny == 1
    field = reshape(field,[Nx Nz_new]);
end
