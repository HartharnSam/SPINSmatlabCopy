function [colaxis, cmap] = choose_caxis(var, data, opts)
%  CHOOSE_CAXIS  choose the colorbar limits and colormap of data based on the field 'var'
%
%  Inputs:
%    var  - field name
%    data - field data
%    opts - plotting options
%
%  Outputs:
%    colaxis   - 2 element vector of colorbar limits
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

% choose color axis limits and colormap based on the field name
ncmap = opts.ncmap;
if strcmp(opts.style,'contour')
    cmap = darkjet(ncmap);    % the default colormap
else
    cmap = temperature(ncmap);    % the default colormap
end
if exist('SD', 'var') ||...
    strcmpi(var, 'KE') ||...
    strcmpi(var, 'diss')
    colaxis = [0 1]*max(data(:));
    cmap = flipud(bone(ncmap));
elseif ~isempty(strfind(var, 'Dye')) ||...
    strcmpi(var, 'Tracer')
    colaxis = [-1 1];
elseif (strcmp(var, 'Density') && mean(data(:)) < 1e-5) ||...
    strcmpi(var, 'U') || ...
    strcmpi(var, 'V') || ...
    strcmpi(var, 'W') ||...
    strncmp(var, 'vort', 4)
    colaxis = [-1 1]*max(abs(data(:)));
elseif strcmp(var, 'Ri')
    colaxis = [0 5];
else
    colaxis = 'auto';
end
