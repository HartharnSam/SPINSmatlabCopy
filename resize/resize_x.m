function field = resize_x(field_name, field, Nx_new, varargin)
% RESIZE_X  Increase or decrease the resolution in the x-dimension.
%
%  Assumptions:
%    - if 2D, data must be in x-z plane
%    - Boundary condition must be free-slip or periodic
%      - To keep the code simple the periodic BC is not optimized
%      - since the field will still be doubled
%
%  Usage:
%    rho_new = resize_x('rho', rho, 512);
%
%  Inputs:
%    'field_name'   - the name of the field to resize
%    'field'        - the field data
%    'Nx_new'       - number of points to change to
%   Optional arguments:
%    'opts'         - a structure holding simulation parameters
%                   - these are: type_x and Lx
%
%  Outputs:
%    'field'        - the resized field
%
%  David Deepwell, 2018

% get number of points in each dimension
sz     = size(field);
Nx_old = sz(1);
Nz     = sz(end);
if length(sz) == 3
    Ny = sz(2);
else
    Ny = 1;
end
% get other parameters
if nargin == 3
    % read from spins.conf
    params = spins_params();
    type_x = params.type_x;
    Lx     = params.Lx;
elseif nargin == 4
    % read from optional input argument
    opts = varargin{1};
    type_x = opts.type_x;
    Lx     = opts.Lx;
else
    error('Incorrect number of input arguments. 3 or 4 are expected.')
end
dx_old = Lx/Nx_old;

% error check expansion type
if ~strcmp(type_x, 'FREE_SLIP') && ~strcmp(type_x, 'FOURIER')
    error('x boundary condition is not FREE_SLIP or FOURIER.')
end

% check sizes
mult = Nx_new/Nx_old;
if mult > 1
    increase_points = true;
else
    increase_points = false;
end

% reshape 2D to match the shape of a 3D field
if Ny == 1
    field = reshape(field, [Nx_old 1 Nz]);
end

% extension types
odd_ext    = @(f) [f; -flipud(f)];
even_ext   = @(f) [f;  flipud(f)];
period_ext = @(f) [f;  f];

% double the field
if strcmp(type_x, 'FOURIER')
    field = period_ext(field);
elseif strcmp(field_name, 'u')
    field = odd_ext(field);
else
    field = even_ext(field);
end

% phase shift because grid is cell-centered, and take fft
dx_new = Lx/Nx_new;
kxs = fftfreq(2*Nx_old, dx_old);
phase_corr = exp(-1i*kxs*(dx_old/2-dx_new/2)).';
field = bsxfun(@times, fft(field), phase_corr);

% pad or truncate
if increase_points
    % pad with zeros
    field = [field(1:Nx_old,:,:); zeros(2*Nx_new-2*Nx_old,Ny,Nz); field(end-Nx_old+1:end,:,:)];
else
    % truncate high frequencies
    field = [field(1:Nx_new,:,:); field(end-Nx_new+1:end,:,:)];
end

% take ifft
field = mult*real(ifft(field));
field = field(1:Nx_new,:,:);
% return field to 2D if it was
if Ny == 1
    field = reshape(field,[Nx_new Nz]);
end
