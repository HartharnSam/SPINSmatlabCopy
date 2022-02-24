function [time, x, data] = spins_hovmoller_z(varname, z_loc, t_inds, varargin)
%  SPINS_HOVMOLLER  plot a space-time plot of the field 'varname' at z=loc
%                   over the times corresponding to the indices t_inds
%
%  Inputs:
%    varname - field name
%    loc     - z location (in units)
%   Optional arguments:
%    t_inds  - vector of time indices (default is all)
%
%  Outputs:
%    time    - time vector
%    z       - space vector
%    data    - field data
%
%  %  David Deepwell, 2018, adapted for z by Sam Hartharn-Evans
%global gdpar

% get grid and parameters
gd.x = xgrid_reader;
gd.z = zgrid_reader;

params = spins_params;

% Optional arguments
if nargin == 2 || isempty(t_inds)
    t_inds = first_output(varname):last_output(varname);
end
% define expected options
d.plot = true;     % create the plot? (bool)
p = inputParser;
addParameter(p,'plot', d.plot, @islogical)
addParameter(p,'Type','standard',@ischar);

parse(p,varargin{:})
% put options into a shorter structure
opts = p.Results;

type = opts.Type;

try
    load wave_characteristics.mat WaveStats time wave_center wavelength_right wavelength_left
catch
    characterize_wave;
    load wave_characteristics.mat WaveStats time wave_center wavelength_right wavelength_left
end
timei = time;

% Assumptions (to be removed in the future):
%   2D
%   loc is along x axis (so plot will be z-t)

% Find x-index and initialize variables
[~, z_ind] = min(abs(gd.z - z_loc), [], 2);
N_t = length(t_inds);
data = zeros(params.Nx, N_t);
z_ind_minmax = min(z_ind):max(z_ind);

% read data
for nn = 1:N_t
    ii = t_inds(nn);
    
    temp_data = spins_reader_new(varname, ii, [], z_ind_minmax);
    for jj = 1:params.Nx
        data(jj, nn) = temp_data(jj, z_ind(jj)+1-z_ind_minmax(1));
    end
    %if isempty(wave_speed(nn))
    %    wave_speed(nn) = wave_speed(end);
    %end
    if strcmpi(type, 'normalised')
        data(:, nn) = data(:, nn)./WaveStats.meanWaveSpeed;
    elseif strcmpi(type, 'WaveFramed')
        data(:, nn) = data(:, nn) - WaveStats.meanWaveSpeed;
    end
    
end

% create space and time variables
time = t_inds*params.plot_interval;
x = gd.x(:,1);

% make plot
if opts.plot
    
    pcolor(x, time, data')
    xlabel('x (m)')
    ylabel('t (s)')
    
    % make plot pretty
    shading flat
    opts.nlevels = 128;
    opts.style = 'pcolor';
    [clim, cmap] = choose_caxis(varname, data, opts);
    colormap(cmap)
    c = colorbar;
    
    % Add guiding lines
    hold on
    plot(wave_center, timei, 'k-')
    plot(wave_center + wavelength_right, timei,'k--')
    plot(wave_center - wavelength_left, timei, 'k--')
    
    % Finish making pretty
    ylabel(c, varname);
    caxis(clim)
    caxis([-.1 .1])
    xlim([x(1) x(end)])
    ylim([time(1) time(end)])
    
    
end
end