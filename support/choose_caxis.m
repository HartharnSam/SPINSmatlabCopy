function [colaxis, cmap] = choose_caxis(var, data, ncmap)
%  CHOOSE_CAXIS  choose the colorbar limits and colormap of data based on the field 'var'
%
%  Inputs:
%    var   - field name
%    data  - field data
%    ncmap - length of colormap
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
cmap = darkjet(ncmap);    % the default colormap
if exist('SD', 'var')
    colaxis = [0 1]*max(data(:));
    cmap = hot(ncmap);
elseif ~isempty(strfind(var, 'Dye')) || strcmpi(var, 'Tracer')
    colaxis = [-1 1];
elseif strcmp(var, 'Density') && mean(data(:)) < 1e-5
    colaxis = [-1 1]*max(abs(data(:)));
elseif strcmpi(var, 'U') || strcmpi(var, 'V') || strcmpi(var, 'W')
    colaxis = [-1 1]*max(abs(data(:)));
elseif strcmp(var, 'KE')
    colaxis = [0 1]*max(data(:));
    cmap = hot(ncmap);
elseif strcmp(var, 'Ri')
    colaxis = [0 5];
else
    colaxis = 'auto';
end 
