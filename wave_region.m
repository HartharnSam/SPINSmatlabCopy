function ax = wave_region(ii, varargin)
% WAVE_REGION    Find the region around the wave based on the values in wave_characteristics.mat
%
%  Usage:
%    ax = wave_region(4);
%    ax = wave_region(4,'mult',3,...);
%
%  Inputs:
%    'ii' - output number to use
%    Optional arguments:
%   Name:   Options              - Description
%	---------------------------------------------------------
%   file:   {string}             - file to read wave characteristics
%   type:   {'variable','fixed'} - type of horizontal spacing
%   mult:   {double, or [m1 m2]} - wavelength multiplier for right and left 
%   x:      {double, or [xL xR]} - vertical boundary locations wrt the wave center
%   z:      {double, or [zB zT]} - horizontal boundary locations
%           x requires type 'fixed', mult requires type 'variable' (default)
%
%  Outputs:
%    'ax' - vector [xL xR zB zT] defining the edges of the wave region
%
%  David Deepwell, 2016

%% get parameters
params = spins_params();
z_bot = 0 + params.min_z;
z_top = params.Lz + params.min_z;

%% Handle Optional arguments
% define expected options
exp_type = {'fixed','variable'};
% define default optional arguments
d.file = 'wave_characteristics.mat';    % file holding wave data
d.type = 'variable';    % type of spacing to use
d.mult = 3;             % horizontal multiplier (left and right)
d.x = 1;                % vertical boundary locations (wrt the wave center)
d.z = [z_bot z_top];    % horizontal boundary locations

% parse optional arguments
p = inputParser;
addParameter(p,'file', d.file, @char)
addParameter(p,'type', d.type, @(x) any(validatestring(x,exp_type)))
addParameter(p,'x', d.x, @isnumeric)
addParameter(p,'z', d.z, @isnumeric)
addParameter(p,'mult', d.mult, @isnumeric)
parse(p,varargin{:})
% put options into a shorter structure
opts = p.Results;

% adjust some options
if length(opts.mult) == 1
    opts.mult = [opts.mult opts.mult];
end
if length(opts.x) == 1
    opts.x = [-opts.x opts.x];
end

%% get wave characteristics
wave_char = load(opts.file);
ind = ii - first_output() + 1;
cntr = wave_char.wave_center(ind);
wl_r =    wave_char.wavelength_right(ind);
wl_l =    wave_char.wavelength_left(ind);

%% horizontal limits
if strcmp(opts.type, 'fixed')
    xL = cntr + opts.x(1);
    xR = cntr + opts.x(2);
elseif strcmp(opts.type, 'variable')
    xL = cntr - wl_l*opts.mult(1);
    xR = cntr + wl_r*opts.mult(2);
end

%% vertical limits
if length(opts.z) == 1
    zB = params.min_z;
    zT = opts.z + params.min_z;
else
    zB = opts.z(1);
    zT = opts.z(2);
end

% adjust if values are outside of domain
if xL < params.min_x
    xL = params.min_x;
end
if xR > params.min_x + params.Lx
    xR = params.min_x + params.Lx;
end
if zB < params.min_z
    zB = params.min_z;
end
if zT > params.min_z + params.Lz
    zT = params.min_z + params.Lz;
end

% output
ax = [xL xR zB zT];
