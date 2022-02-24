% spins_plotoptions parses and creates the optional arguments for spins_plot2d
%
%  David Deepwell, 2015

% define expected options 
exp_dimen = {'X','x','Y','y','Z','z'};
exp_style = {'pcolor','contourf','contour'};
% define defaults 
d.dimen = 'y';          % dimension
d.slice = [0 0];        % cross-section. opt. arg. must be single number
d.axis = 0;             % axis to plot. 0 denotes use of full domain
d.style = 'pcolor';     % plotting style
d.xskp = 1;             % x-grid points to skip
d.yskp = 1;             % y	"
d.zskp = 1;             % z	"
d.fnum = 1;             % figure window number to use
d.var2 = 'Density';     % secondary field to plot
d.nlevels2 = 6;         % number of contours of secondary field
d.ncontourf = 65;       % plotting regions in contourf style
d.ncontour = 20;        % contours in contour style
d.npcolor = 128;        % number of levels to use in colorbar for pcolor
d.nlevels = 0;          % levels of colormap (0 is placeholder, value set below)
d.clim = 0;             % colour axis limits to use -
                        % 0 uses default in choose_caxis function
d.colorbar = true;      % colorbar? (bool)
d.trim = false;         % trims values outside clim range to be within it
d.axisstyle = 'normal'; % the axis style to use
d.visible = true;       % make plot visible or not (bool)
d.speed = -1;           % wave speed for streamlines
d.streamline_heatmap= true; % Add a heatmap of the magnitude of velocity under the streamlines?
d.streamline_density= 2;    % Density of streamlines
d.savefig = false;      % save figure? (bool)
d.filename = 'filename';	% name of file to save
d.dir = 'figures';      % name of file to save

% define validation functions
fnum_error_msg = 'fnum option was not understood. Provide either an integer or "New".';
time_error_msg = sprintf(['Time option was not understood.\n Provide either an integer ',...
    '(for output number)\n or a string with time (ex. "5")']);
validation_fnum = @(x) assert(isnumeric(x) || strcmpi(x, 'New'), fnum_error_msg);
validation_time = @(x) assert(isnumeric(x) || ischar(x), time_error_msg);
validation_clim = @(x) assert(isnumeric(x) || strcmpi(x, 'auto'));

% parse time input type
validation_time(t_index)
% if it's a character, then it's passing a time in seconds
% need to convert to output number
if ischar(t_index)
    dt = params.plot_interval;
    if isfield(params, 'restart_time')
        t_0 = params.restart_time;
    else
        t_0 = 0;
    end
    if isfield(params, 'restart_sequence')
        seq_0 = params.restart_sequence;
    else
        seq_0 = 0;
    end
    time = str2double(t_index);
    % find nearest output time
    t_index = round((time-t_0)/dt + seq_0);

    % check if index is in bounds of simulation outputs
    if t_index > last_output()
        t_index = last_output();
    end
    if t_index < first_output()
        t_index = first_output();
    end
end

% parse optional arguments
p = inputParser;
addParameter(p,'dimen', d.dimen, @(x) any(validatestring(x,exp_dimen)))
addParameter(p,'slice', d.slice, @isnumeric)
addParameter(p,'fnum', d.fnum, validation_fnum)
addParameter(p,'style', d.style, @(x) any(validatestring(x,exp_style)))
addParameter(p,'var2', d.var2, @ischar)
addParameter(p,'nlevels2', d.nlevels2, @isnumeric)
addParameter(p,'speed', d.speed, @isnumeric)
addParameter(p,'streamline_heatmap', d.streamline_heatmap, @islogical)
addParameter(p,'streamline_density', d.streamline_density, @isnumeric)
addParameter(p,'xskp', d.xskp, @isnumeric)
addParameter(p,'yskp', d.yskp, @isnumeric)
addParameter(p,'zskp', d.zskp, @isnumeric)
addParameter(p,'axis', d.axis, @isnumeric)
addParameter(p,'nlevels', d.nlevels, @isnumeric)
addParameter(p,'clim', d.clim, validation_clim)
addParameter(p,'trim', d.trim, @islogical)
addParameter(p,'axisstyle', d.axisstyle, @ischar)
addParameter(p,'colorbar', d.colorbar, @islogical)
addParameter(p,'visible', d.visible, @islogical)
addParameter(p,'savefig', d.savefig, @islogical)
addParameter(p,'filename', d.filename, @ischar)
addParameter(p,'dir', d.dir, @ischar)
parse(p,varargin{:})

% put options into a shorter structure
opts = p.Results;

% choose default cross-section slice based on which dimension is plotted
if params.ndims == 3
    if length(opts.slice) == 2
        if strcmpi(opts.dimen, 'X')
            opts.slice = sum(params.xlim)/2;
        elseif strcmpi(opts.dimen, 'Y')
            opts.slice = sum(params.ylim)/2;
        elseif strcmpi(opts.dimen, 'Z')
            opts.slice = sum(params.zlim)/2;
        end
    end
end
% fix if slice is outside of domain
if params.ndims == 3
    if strcmpi(opts.dimen, 'X')
        if opts.slice > params.xlim(2)
            opts.slice = params.xlim(2);
        elseif opts.slice < params.xlim(1)
            opts.slice = params.xlim(1);
        end
    elseif strcmpi(opts.dimen, 'Y')
        if opts.slice > params.ylim(2)
            opts.slice = params.ylim(2);
        elseif opts.slice < params.ylim(1)
            opts.slice = params.ylim(1);
        end
    elseif strcmpi(opts.dimen, 'Z')
        if opts.slice > params.zlim(2)
            opts.slice = params.zlim(2);
        elseif opts.slice < params.zlim(1)
            opts.slice = params.zlim(1);
        end
    end
end

% make file name more appropriate if not given
if strcmp(opts.filename, 'filename')
    filename = var;
else
    filename = opts.filename;
end

% don't plot secondary field if it matches the primary field
var_name = strrep(var, 'Mean ', '');
if strcmp(var_name, 'rho')
    var_name = 'Density';
end
if strcmp(var_name, opts.var2) || ...
        strncmp(reverse(var), 'mottob', 6) || ... % bottom surface
        strncmp(reverse(var), 'pot', 3) % top surface
    opts.var2 = 'None';
elseif strncmp(reverse(var), 'ms', 2) || ...    % ends in sm - Spanwise Mean
        strncmp(reverse(var), 'zx', 2)          % ends in xz - a cross-section
    if strcmp(var, ['rho-',var(end-1:end)])
        opts.var2 = 'None';
    else
        opts.var2 = ['rho-',var(end-1:end)];
    end
end

% adjust if plotting bottom/top surface
if strncmp(reverse(var), 'mottob', 6) || strncmp(reverse(var), 'pot', 3)
    opts.dimen = 'z';
    if strncmp(reverse(var), 'mottob', 6)
        opts.slice = params.zlim(1);
    elseif strncmp(reverse(var), 'pot', 3)
        opts.slice = params.zlim(2);
    end
end

% change default colormap length depending on plotting style
if strcmp(opts.style, 'contourf') && opts.nlevels == 0
    opts.nlevels = d.ncontourf;
elseif strcmp(opts.style, 'contour') && opts.nlevels == 0
    opts.nlevels = d.ncontour;
elseif strcmp(opts.style, 'pcolor')
    opts.nlevels = d.npcolor;
end

% are we plotting a streamline?
plotting_streamlines1 = ~isempty(strfind(lower(var),       'streamline'));
if plotting_streamlines1
    opts.var2 = 'None';
end
plotting_streamlines2 = ~isempty(strfind(lower(opts.var2), 'streamline'));

% get indices and grid for plotting
[nx, ny, nz, xvar, yvar, zvar, plotaxis] = get_plot_points(opts);
