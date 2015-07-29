function pltinfo = spins_plot2d(var, t_index, varargin)
% SPINS_PLOT2D  Plot cross-sections, averages or standard deviations of variables.
%
%   pltinfo = spins_plot2d(var, t_i) plots var at t_i
%   pltinfo = spins_plot2d(var, t_i, 'opt1', val1, ...) plots var at t_i with option 'opt1' as val1
%   
%   % Optional arguments:
%	Name:	Options			- Description
%       dimen:	{'X','Y','Z'}		- dimension to take cross-section
%       slice:	{double}		- location to take cross-section
%       style:  {'pcolor','contourf','contour'}
%		- type of plot
%	axis:	{[x1 x2 z1 z2]}		- domain to plot
%	xskp:	{integer}		- x-grid points to skip in plot
%	yskp:	{integer}		- y-grid     "
%	zskp:	{integer}		- z-grid     "
%       fnum:	{integer}		- figure window to make plot
%	colorbar:  {boolean}		- plot colorbar?
%	savefig:   {boolean}		- save figure in figure file?
%	visible:   {boolean}		- make figure visible?
%       ncontourf: {integer}		- contours to use for contourf plot
%       ncontour:  {integer}		- contours to use for contour plot
%	cont2:	{field name}		- secondary field to plot as contours
%	ncont2:	{integer}		- contours to use for secondary field
%
%    David Deepwell, 2015
global gdpar

% get grid and parameters
%gd = gdpar.gd;
params = gdpar.params;

% choose correct plot type
if strcmp(params.mapped_grid,'false');
    pltinfo = spins_plot2dunmapped(var,t_index,varargin);
elseif strcmp(params.mapped_grid,'true');
    pltinfo = spins_plot2dmapped(var,t_index,varargin);
end

