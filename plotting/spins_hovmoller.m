function [time, z, data] = spins_hovmoller(varname, x_loc, t_inds, varargin)
%  SPINS_HOVMOLLER  plot a space-time plot of the field 'varname' at x=loc
%                   over the times corresponding to the indices t_inds
%
%  Inputs:
%    varname - field name
%    loc     - x location
%   Optional arguments:
%    t_inds  - vector of time indices (default is all)
%
%  Outputs:
%    time    - time vector
%    z       - space vector
%    data    - field data
%
%  David Deepwell, 2018

%global gdpar

% get grid and parameters
gd.x = xgrid_reader;
gd.z = zgrid_reader;

params = spins_params;

% Optional arguments
if nargin == 2
    t_inds = first_output(varname):last_output(varname);
elseif isempty(t_inds)
    t_inds = first_output(varname):last_output(varname);
end
% define expected options
d.plot = false;     % create the plot? (bool)
p = inputParser;
addParameter(p,'plot', d.plot, @islogical)
parse(p,varargin{:})
% put options into a shorter structure
opts = p.Results;

% Assumptions (to be removed in the future):
%   2D
%   loc is along x axis (so plot will be z-t)

% Find x-index and initialize variables
x_ind = nearest_index(gd.x(:,1), x_loc);
N_t = length(t_inds);
data = zeros(params.Nz, N_t);

% read data
for nn = 1:N_t
    ii = t_inds(nn);
    data(:,nn) = spins_reader_new(varname, ii, x_ind, 1, []);
end

% create space and time variables
time = t_inds*params.plot_interval;
z = gd.z(x_ind,:);

% make plot
if opts.plot
    figure(54)
    clf
    pcolor(time, z, data)
    xlabel('t (s)')
    ylabel('z (m)')

    % make plot pretty
    shading flat
    opts.nlevels = 128;
    opts.style = 'pcolor';
    [clim, cmap] = choose_caxis(varname, data, opts);
    colormap(cmap)
    colorbar
    caxis(clim)
end
