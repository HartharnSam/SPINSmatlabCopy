function [] = print_figure(filename, varargin) 
% PRINT_FIGURE    Prints the figure (given by fig_hand) 
%
%  Usage:
%    print_figure()
%    print_figure('rho')
%    print_figure('rho', 'opt1', val1, 'opt2', val2, ...)
%
%  Inputs:
%    'filename' - the name of the file (is optional)
%
%    Optional arguments:
%   Name:       Options                    - Description
%   ----------------------------------------------------
%   fig_hand:   {integer or figure handle} - figure to use
%   format:     {'pdf', 'png', 'eps', ...} - file format
%   units:      {'Inches', 'cm', ...}      - paper units
%   size:       {[width height]}           - vector of dimensions
%   res:        {integer}                  - image resolution
%   renderer:   {'painters', 'opengl'}     - renderer
%
%  Outputs:
%    - none
%
% David Deepwepll, 2016

%% Manage inputs
% set default inputs
if ~exist('filename', 'var')
    filename = 'fig';
end
d.fig_hand = 0;     % figure handle, 0 is placeholder
d.format = 'png';   % file format
d.units = 'Inches'; % figure and paper units
d.size = 0;         % size of figure, 0 is placeholder
d.res = 500;        % image resolution
d.renderer = 'painters';    % renderer

% parse optional arguments
p = inputParser;
addParameter(p, 'fig_hand', d.fig_hand);
addParameter(p, 'format', d.format);
addParameter(p, 'units', d.units);
addParameter(p, 'size', d.size, @isvector);
addParameter(p, 'res', d.res, @isnumeric);
addParameter(p, 'renderer', d.renderer);
parse(p, varargin{:})
% put options into a shorter structure
opts = p.Results;
fig_hand = opts.fig_hand;

%% Adjust based on given options
if fig_hand == 0
    fig_hand = gcf;
elseif isnumeric(fig_hand)
    fig_hand = figure(fig_hand);
end
% set units
set(fig_hand, 'Units', opts.units,...
              'PaperUnits', opts.units);
% set size
if opts.size == 0
    pos = get(fig_hand, 'Position');
    opts.size = pos(3:4);
end
fig_hand.PaperSize = opts.size;
% set position
fig_hand.Position = [1 1 opts.size];
fig_hand.PaperPosition = [0 0 opts.size];
fig_hand.PaperPositionMode = 'manual';

% save figure
print(fig_hand, filename,...
      ['-d',opts.format],...
      ['-r',num2str(opts.res)],...
      ['-',opts.renderer])
