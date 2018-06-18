function field = resize_y(field_name, field, Ny_new, varargin)
% RESIZE_Y  Increase or decrease the resolution in the y-dimension.
%
%  Assumptions:
%    - field must 3D
%    - Boundary condition must be free-slip or periodic
%      - To keep the code simple the periodic BC is not optimized
%      - since the field will still be doubled
%
%  Usage:
%    rho_new = resize_y('rho', rho, 512);
%
%  Inputs:
%    'field_name'   - the name of the field to resize
%    'field'        - the field data
%    'Ny_new'       - number of points to change to
%   Optional arguments:
%    'opts'         - a structure holding simulation parameters
%                   - these are: type_y and Ly
%
%  Outputs:
%    'field'        - the resized field
%
%  David Deepwell, 2018

% get number of points in each dimension
sz     = size(field);
Nx     = sz(1);
Ny_old = sz(2);
Nz     = sz(3);
% get other parameters
if nargin == 3
    % read from spins.conf
    params = spins_params();
    type_y = params.type_y;
    Ly     = params.Ly;
elseif nargin == 4
    % read from optional input argument
    opts = varargin{1};
    type_y = opts.type_y;
    Ly     = opts.Ly;
else
    error('Incorrect number of input arguments. 3 or 4 are expected.')
end
dy_old = Ly/Ny_old;

% error check expansion type
if ~strcmp(type_y, 'FREE_SLIP') && ~strcmp(type_y, 'FOURIER')
    error('y boundary condition is not FREE_SLIP or FOURIER.')
end

% check sizes
mult = Ny_new/Ny_old;
if mult > 1
    increase_points = true;
else
    increase_points = false;
end

% permute to put extending dimension in the first dimension
field = permute(field, [2 1 3]);

% extension types
odd_ext    = @(f) [f; -flipud(f)];
even_ext   = @(f) [f;  flipud(f)];
period_ext = @(f) [f;  f];

% double the field
if strcmp(type_y, 'FOURIER')
    field = period_ext(field);
elseif strcmp(field_name, 'v')
    field = odd_ext(field);
else
    field = even_ext(field);
end

% phase shift because grid is cell-centered, and take fft
dy_new = Ly/Ny_new;
kys = fftfreq(2*Ny_old, dy_old);
phase_corr = exp(-1i*kys*(dy_old/2-dy_new/2)).';
field = bsxfun(@times, fft(field), phase_corr);

% pad or truncate
if increase_points
    % pad with zeros
    field = [field(1:Ny_old,:,:); zeros(2*Ny_new-2*Ny_old,Nx,Nz); field(end-Ny_old+1:end,:,:)];
else
    % truncate high frequencies
    field = [field(1:Ny_new,:,:); field(end-Ny_new+1:end,:,:)];
end

% take ifft
field = mult*real(ifft(field));
field = field(1:Ny_new,:,:);
% permute back
field = permute(field, [2 1 3]);
