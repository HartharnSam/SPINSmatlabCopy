function [] = print_all(filenames, varargin)
% PRINT_ALL    Prints all the figures in the workspace
%
%  Usage:
%    print_all()
%    print_all('figs')
%    print_all('figs', 'opt1', val1, 'opt2', val2, ...)
%
%  Inputs:
%    'filenames' - the name of the files
%
%    Optional arguments:
%   Name:       Options                    - Description
%   ----------------------------------------------------
%   format:     {pdf, png, eps, ...}       - file format
%   units:      {Inches, cm, ...}          - paper units
%   size:       {[width height]}           - vector of dimensions
%   res:        {integer}                  - image resolution
%
%  Outputs:
%    - none
%
% David Deepwepll, 2017

%% Manage inputs
% set default inputs
if ~exist('filenames', 'var')
    filenames = 'fig';
end
d.format = 'png';   % file format
d.units = 'Inches'; % figure and paper units
d.size = 0;         % size of figure, 0 is placeholder
d.res = 500;        % image resolution

% parse optional arguments
p = inputParser;
addParameter(p, 'format', d.format);
addParameter(p, 'units', d.units);
addParameter(p, 'size', d.size, @isvector);
addParameter(p, 'res', d.res, @isnumeric);
parse(p, varargin{:})
% put options into a shorter structure
opts = p.Results;

% check/make figures directory
if ~exist('figures', 'dir')
    mkdir('figures')
end
cd('figures')

%% save to directory
figs = findall(0,'type','figure');
for ii = 1:length(figs)
    % choose figure and print
    fig_hand = figs(ii);
    filename = [filenames,'_',num2str(fig_hand.Number)];
    print_figure(filename,...
                 'fig_hand', fig_hand,...
                 'format', opts.format,...
                 'units', opts.units,...
                 'size', opts.size,...
                 'res', opts.res);
    % print information
    disp(['Figure ',num2str(fig_hand.Number),' printed.'])
    completion(ii, length(figs))
end
cd('..')
