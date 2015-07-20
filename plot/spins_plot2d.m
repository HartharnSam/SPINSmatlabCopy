function pltinfo = spins_plot2d(var, t_index, varargin)
% SPINS_PLOT2D  Plot cross-sections, averages or standard deviations of variables.
%
%   pltinfo = spins_plot2d(var,t_i) plots var at t_i
%   other options
%
%    David Deepwell, 2015
global gdpar

% get grid and parameters
gd = gdpar.gd;
params = gdpar.params;

% choose correct plot type
if strcmp(params.mapped_grid,'false');
    pltinfo = spins_plot2dunmapped(var,t_index,varargin);
elseif strcmp(params.mapped_grid,'true');
    pltinfo = spins_plot2dmapped(var,t_index,varargin);
end

