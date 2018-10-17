function [clim, cmap] = choose_caxis(var, data, opts)
%  CHOOSE_CAXIS  choose the colorbar limits and colormap of data based on the field 'var'
%
%  Inputs:
%    var  - field name
%    data - field data
%    opts - plotting options
%
%  Outputs:
%    clim   - 2 element vector of colorbar limits
%    cmap      - colormap
%
%  David Deepwell, 2015

% remove extra information from the plotted field
if strncmp(var, 'Mean', 4)
    var = strsplit(var, 'Mean ');
    var = var{end};
elseif strncmp(var, 'SD', 2)
    SD = true;
    var = strsplit(var, 'SD ');
    var = var{end};
elseif strncmp(var, 'Scaled SD', 9)
    SD = true;
    var = strsplit(var, 'Scaled SD ');
    var = var{end};
end

% default maps
nlevels = opts.nlevels;
if strcmp(opts.style,'contour')
    cmap = darkjet(nlevels);
else
    cmap = temperature(nlevels);
end
% choose color axis limits and colormap based on the field name
if exist('SD', 'var') ...
    || strcmpi(var, 'KE') ...
    || strcmpi(var, 'speed') ...
    || strcmpi(var, 'diss') ...
    || strcmpi(var, 'enst')
    clim = [0 1]*max(data(:));
    if ~strcmp(opts.style,'contour')
        cmap = cmocean('-grey',nlevels);
    end
elseif ~isempty(strfind(lower(var), 'dye')) ...
    || strcmpi(lower(var), 'tracer')
    clim = [0 1];
    cmap = cmocean('tempo', nlevels);
elseif ((strcmpi(var, 'Density') || strcmp(var, 'rho')) ...
    && mean(data(:)) < 1)
    clim = [-1 1]*max(abs(data(:)));
    cmap = flipud(cmap);
elseif strcmpi(var, 'U') ...
    || strcmpi(var, 'V') ...
    || strcmpi(var, 'W') ...
    || strncmp(var, 'vort', 4)
    clim = [-1 1]*max(abs(data(:)));
elseif strcmp(var, 'Ri')
    clim = [0 5];
    if ~strcmp(opts.style,'contour')
        cmap = cmocean('grey',nlevels);
    end
else
    clim = 'auto';
end
