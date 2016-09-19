% spinsplotoptions creates the optional arguments for spins_plot2d
%
%  David Deepwell, 2015

% define expected options 
exp_dimen = {'X','x','Y','y','Z','z'};
exp_style = {'pcolor','contourf','contour'};
% define defaults 
d.dimen = 'Y';			% dimension
d.slice = [0 0];		% cross-section. opt. arg. must be single number
d.fnum = 1;		    	% figure window number to use
d.style = 'contourf';	% plotting style
d.ncontourf = 50;		% plotting regions in contourf style
d.ncontour = 20;		% contours in contour style
d.cont2 = 'Density';	% secondary field to plot
d.ncont2 = 6;			% contours of secondary field
d.speed = -1;			% wave speed for streamlines
d.xskp = 1;	    		% x-grid points to skip
d.yskp = 1;	    		% y	"
d.zskp = 1;		    	% z	"
d.axis = 0;			    % axis to plot. 0 denotes use of full domain
d.ncmap = 128;          % length of colormap (only for pcolor)
d.colaxis = 0;          % colour axis limits to use -
                        % 0 uses default in choose_caxis function
d.trim = false;         % trims values outside colaxis range to be within it
d.colorbar = true;      % colorbar? (bool)
d.visible = true;		% make plot visible or not (bool)
d.savefig = false;		% save figure? (bool)
d.filename = 'filename';	% name of file to save
d.dir = 'figures';      % name of file to save

% define validation functions
validation_fnum = @(x) assert(isnumeric(x) || strcmpi(x, 'New'),...
                'fnum option was not understood. Provide either an integer or "New".');

% parse options
p = inputParser;
addParameter(p,'dimen', d.dimen, @(x) any(validatestring(x,exp_dimen)))
addParameter(p,'slice', d.slice, @isnumeric)
addParameter(p,'fnum', d.fnum, validation_fnum)
addParameter(p,'style', d.style, @(x) any(validatestring(x,exp_style)))
addParameter(p,'ncontourf', d.ncontourf, @isnumeric)
addParameter(p,'ncontour', d.ncontour, @isnumeric)
addParameter(p,'cont2', d.cont2, @ischar)
addParameter(p,'ncont2', d.ncont2, @isnumeric)
addParameter(p,'speed', d.speed, @isnumeric)
addParameter(p,'xskp', d.xskp, @isnumeric)
addParameter(p,'yskp', d.yskp, @isnumeric)
addParameter(p,'zskp', d.zskp, @isnumeric)
addParameter(p,'axis', d.axis, @isnumeric)
addParameter(p,'ncmap', d.ncmap, @isnumeric)
addParameter(p,'colaxis', d.colaxis, @isnumeric)
addParameter(p,'trim', d.trim, @islogical)
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

% change default colormap length depending on plotting style
if strcmp(opts.style, 'contourf')
    opts.ncmap = opts.ncontourf;
elseif strcmp(opts.style, 'contour')
    opts.ncmap = opts.ncontour;
end

% get indices and grid for plotting
[nx, ny, nz, xvar, yvar, zvar, plotaxis] = get_plot_points(opts);
