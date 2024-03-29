function pltinfo = spins_plot2d(var, t_index, varargin)
%  SPINS_PLOT2D  Plot cross-sections, averages or standard deviations of variables.
%
%  Usage:
%	pltinfo = spins_plot2d(var, t_i) plots var at t_i
%	pltinfo = spins_plot2d(var, 't') plots var at nearest output to time t (in seconds)
%	pltinfo = spins_plot2d(var, t_i, 'opt1', val1, ...) plots var at t_i with option 'opt1' as val1
%   
%  Inputs:
%    'var' may be of different forms:
%	any field in the working directory ('rho','u',...)
%   'Density'       reads rho or calculates it from Salt and Temp
%   'KE'            local kinetic energy
%   'speed'         magnitude of the local velocity
%   'Ri'            gradient Richardson number
%   'Streamlines'   streamlines in the x-z plane
%   'Mean ...'      takes the spanwise mean of ...
%   'SD ...'        takes the spanwise standard deviation of ...
%   'Scaled SD ...' scales SD ... by the maximum of ...
%
%    't_index' may be:
%   an integer for a particular output
%   a vector of outputs (ie. 0:10) will plot each output successively in the same figure
%   a string containing the time (ex. '15')
%
%    Optional arguments:
%	Name:	Options			- Description (defaults are in spins_plotoptions.m)
%	---------------------------------------------------------
%   dimen:  {'X','Y','Z'}       - dimension to take cross-section
%   slice:  {double}            - location to take cross-section
%   axis:   {[x1 x2 z1 z2]}     - domain to plot
%   axes:   {handle}            - axes to plot on
%   style:  {'pcolor','contourf','contour'}     - type of plot
%   xskp:   {integer}           - x-grid points to skip in plot
%   yskp:   {integer}           - y-grid     "
%   zskp:   {integer}           - z-grid     "
%   fnum:   {integer}           - figure window to make plot
%   var2:   {field name}        - secondary field to plot as contours
%   nlevels2:  {integer}        - number of contours to use for secondary field
%   nlevels:   {integer}        - number of levels in plot
%   clim:      {[c1 c2]}        - color axis limits to use
%   colorbar:  {boolean}        - plot colorbar?
%   trim:      {boolean}        - trims values outside clim range
%   visible:   {boolean}        - make figure visible?
%   speed:     {double}         - wave speed to subtract from flow in streamlines plot
%   savefig:   {boolean}        - save figure in figure file?
%   filename:  {string}         - name of file of saved figure
%   dir:       {string}         - name of relative directory to save figure
%
%  Outputs:
%    'pltinfo'	- a structure containing the plotted fields 
%
%  David Deepwell, 2015
global gdpar
if isempty(gdpar)
    gd.x = xgrid_reader();
    gd.z = zgrid_reader();
    params = spins_params;
    gdpar.gd = gd;
    gdpar.params = params;
    gdpar.params.ndims = ndims(gd.x);

end

% get grid and parameters
split_gdpar
if ~strcmpi(params.mapped_grid,'true') && ~isvector(gd.x)
    gd = get_vector_grid(gd);
end

% check that the grid is the same size as that listed in spins.conf
par = spins_params();
if par.Nx ~= params.Nx || ...
        par.Ny ~= params.Ny || ...
        par.Nz ~= params.Nz
    error('The currently loaded grid is incorrect. Rerun spins_grid or spins_gridparams')
end

% set plotting options
spins_plotoptions

% open new figure
if strcmpi(opts.fnum, 'New')
    fighand = figure;
else
    fighand = figure(opts.fnum);
end
% figure visibility options
if opts.visible == false
    set(fighand, 'Visible', 'off')
end

for ii = t_index
    clf
    hold on
    % Title
    plot_title = var;
    is_spanwise = strncmp(var,'Mean',4) || strncmp(var,'SD',2) || strncmp(var, 'Scaled SD',9);
    if is_spanwise
        plot_title = ['Spanwise ', plot_title];
    end
    % add cross-section information
    if params.ndims == 3 && ~is_spanwise
        plot_title = [plot_title,', ',opts.dimen,'=',num2str(opts.slice),' m'];
    end
    % add time in seconds or output number to title
    if strncmp(reverse(var), 'ms', 2) || strncmp(reverse(var), 'zx', 2) || ...
        strncmp(reverse(var), 'mottob', 6) || strncmp(reverse(var), 'pot', 3)
        times = read_slice_times();
        time = times(ii);
        plot_title = [plot_title,', t=',num2str(time),' s'];
    elseif isfield(params, 'plot_interval') && isfield(params, 'restart_time') && ...
       isfield(params, 'restart_sequence')
        time = params.restart_time + (ii-params.restart_sequence)*params.plot_interval;
        plot_title = [plot_title,', t=',num2str(time),' s'];
    elseif isfield(params, 'plot_interval')
        time = ii*params.plot_interval;
        plot_title = [plot_title,', t=',num2str(time),' s'];
    else
        time = 'n/a';
        plot_title = [plot_title,', t_n=',int2str(ii)];
    end
    title(plot_title);
    % axis labels
    if strcmpi(opts.dimen, 'X')
        xlabel('y (m)'), ylabel('z (m)')
    elseif strcmpi(opts.dimen, 'Y')
        xlabel('x (m)'), ylabel('z (m)')
    elseif strcmpi(opts.dimen, 'Z')
        xlabel('x (m)'), ylabel('y (m)')
    end

    % get data to plot
    data1 = spins_readdata(var, ii, nx, ny, nz);

    % if mapped grid and taking horizontal opts.slice, then find interpolation
    if strcmpi(opts.dimen, 'Z') && strcmpi(params.mapped_grid, 'true')
        if ~(strncmp(reverse(var), 'mottob', 6) || strncmp(reverse(var), 'pot', 3)) % not bottom or top
            [xvar, yvar, data1] = get_fixed_z(xvar, yvar, zvar, data1, opts.slice);
        else
            xvar = squeeze(gd.x(:,:,1));
            yvar = squeeze(gd.y(:,:,1));
        end
    end
    % find rectilinear grid for mapped 

    if (plotting_streamlines1 || plotting_streamlines2) && strcmpi(params.mapped_grid, 'true')
        if params.ndims == 3
            gd_select.x = gd.x(nx,ny,nz);
            gd_select.y = gd.y(nx,ny,nz);
            gd_select.z = gd.z(nx,ny,nz);
            gd_rect = get_rectilinear_grid(gd_select);
            xvar_slice = squeeze(gd_rect.x)';
            yvar_slice = squeeze(gd_rect.z)';
        else
            gd_select.x = gd.x(nx,nz);
            gd_select.z = gd.z(nx,nz);
            gd_rect = get_rectilinear_grid(gd_select);
            xvar_slice = gd_rect.x';
            yvar_slice = gd_rect.z';
        end
    end

    % transpose unmapped data
    if strcmpi(params.mapped_grid, 'false') && ~plotting_streamlines1
        data1 = data1';
    end
    % remove points outside of desirable plotting range (typically from spectral aliasing)
    % this sets more contour levels into region that matters
    if opts.trim == true
        if opts.clim == 0
            error('Trim requires an axis range to trim into.')
        else
            data1(data1>opts.clim(2)) = opts.clim(2);
            data1(data1<opts.clim(1)) = opts.clim(1);
        end
    end

    % choose plotting style (contourf may take up less memory,
    % but can be slower than pcolor)
    if plotting_streamlines1
        if opts.speed == -1
            prompt = 'Provide a sensible wave speed in m/s: ';
            uwave = input(prompt);
        else
            uwave = opts.speed;
        end
        disp(['background speed = ',num2str(uwave),' m/s'])
        u1 = data1(:,:,1) - uwave;
        u2 = data1(:,:,2);
        data1(:,:,1) = u1; % update data1
        vel_mag = sqrt(u1.^2 + u2.^2);
        if strcmp(params.mapped_grid, 'true')
            if opts.streamline_heatmap
                pcolor(xvar_slice, yvar_slice, vel_mag')
                colorbar
            end
            p_hand = streamslice(xvar_slice,yvar_slice,u1',u2',opts.streamline_density,'noarrows','cubic');
        else
            if opts.streamline_heatmap
                pcolor(xvar, yvar, vel_mag')
                colorbar
            end
            p_hand = streamslice(xvar,      yvar,      u1',u2',opts.streamline_density,'noarrows','cubic');
        end
        for mm = 1:length(p_hand)
            p_hand(mm).Color = [0.8500, 0.3250, 0.0980];
        end
        if opts.streamline_heatmap
            var2col = [0.2810, 0.7250, 0.9130];
        else
            var2col = [0 0 0];
        end
    elseif strcmpi(opts.style,'pcolor')
        p_hand = pcolor(xvar,yvar,data1);
        var2col = [0 0 0];
    elseif strcmpi(opts.style,'contourf')
        [~,p_hand] = contourf(xvar,yvar,data1,opts.nlevels);
        var2col = [0 0 0];
    elseif strcmpi(opts.style,'contour')
        [~,p_hand] = contour(xvar,yvar,data1,opts.nlevels);
        var2col = [0 0 0];
    end

    % get caxis limits
    [clim, cmap] = choose_caxis(var, data1, opts);
    % use user defined caxis if specified
    if length(opts.clim) ~= 1
        clim = opts.clim;
    end

    % add extra information
    shading flat
    if strcmpi(opts.style,'contourf')
        set(p_hand,'LineColor','none')
    end
    colormap(cmap)
    if ~strcmpi(clim, 'auto')
        caxis(clim);
    end
    if opts.colorbar == true && ~plotting_streamlines1
        colorbar
    end

    % add contours of another field
    if ~strcmpi(opts.var2, 'None')
        if (strncmp(var,'Mean',4) || strncmp(var,'SD',2))
            % choose Mean of field if primary field is Mean or SD
            var2 = ['Mean ',opts.var2];
        else
            var2 = opts.var2;
        end
        if strcmp(var, var2)            % read in data only if the field is different
            data2 = data1;
        else
            data2 = spins_readdata(var2, ii, nx, ny, nz);
            if strcmpi(opts.dimen, 'Z') && strcmpi(params.mapped_grid, 'true')
                [xvar, yvar, data2] = get_fixed_z(xvar, yvar, zvar, data2, opts.slice);
            end
            if strcmpi(params.mapped_grid, 'false') && ~plotting_streamlines1
                data2 = data2';
            end
            if plotting_streamlines2
                warning('Streamlines not built for secondary variable... yet') 
                if opts.speed == -1
                    prompt = 'Provide a sensible wave speed in m/s: ';
                    uwave = input(prompt);
                else
                    uwave = opts.speed;
                end
                disp(['background speed = ',num2str(uwave),' m/s'])
                u1 = data2(:,:,1) - uwave;
                u2 = data2(:,:,2);
                data2(:,:,1) = u1;
                streamslice(xvar,yvar,u1',u2',2,'noarrows','cubic')
            else
                [~,p2_hand] = contour(xvar,yvar,data2,opts.nlevels2);
                p2_hand.LineColor = var2col;
            end
        end
    end
    if strcmp(var, 'Ri')
        contour(xvar, yvar, data1, [1 1]*0.25, 'r-');
    end

    % add contour of hill if grid is mapped
    if strcmpi(params.mapped_grid,'true') && strcmpi(opts.dimen, 'Y')
        hill_nx = nx(1):nx(end); % don't skip points when plotting the hill
        hill_nz = nz(1):nz(end); % don't skip points when plotting the hill
        % check whether the top, the bottom, or both were mapped
        if params.ndims == 2
            bottom = gd.z(round(linspace(1,params.Nx,10)),1);
            top    = gd.z(round(linspace(1,params.Nx,10)),params.Nz);
            right  = gd.x(params.Nx,round(linspace(1,params.Nz,10)));
            left   = gd.x(1,        round(linspace(1,params.Nz,10)));
        else
            bottom = gd.z(round(linspace(1,params.Nx,10)),1,1);
            top    = gd.z(round(linspace(1,params.Nx,10)),1,params.Nz);
            right  = gd.x(params.Nx,1,round(linspace(1,params.Nz,10)));
            left   = gd.x(1,        1,round(linspace(1,params.Nz,10)));
        end
        bratio =  max(bottom) - min(bottom);
        tratio =  max(top)    - min(top);
        rratio =  max(right)  - min(right);
        lratio =  max(left)   - min(left);
        % plot bottom if mapped
        if bratio ~= 0
            if params.ndims == 3
                hill   = squeeze(gd.z(hill_nx,1,1));
                hill_x = squeeze(gd.x(hill_nx,1,1));
            elseif params.ndims == 2
                hill   = squeeze(gd.z(hill_nx,1));
                hill_x = squeeze(gd.x(hill_nx,1));
            end
            plot(hill_x,hill,'k')
        end
        % plot top if mapped
        if tratio ~= 0
            if params.ndims == 3
                hill   = squeeze(gd.z(hill_nx,1,params.Nz));
                hill_x = squeeze(gd.x(hill_nx,1,params.Nz));
            elseif params.ndims == 2
                hill   = squeeze(gd.z(hill_nx,params.Nz));
                hill_x = squeeze(gd.x(hill_nx,params.Nz));
            end
            plot(hill_x,hill,'k', 'LineWidth', 1)
        end
        % if right side is mapped (or tank is rotated)
        if rratio ~= 0
            if params.ndims == 3
                hill   = squeeze(gd.z(params.Nx,1,hill_nz));
                hill_x = squeeze(gd.x(params.Nx,1,hill_nz));
            elseif params.ndims == 2
                hill   = squeeze(gd.z(params.Nx,hill_nz));
                hill_x = squeeze(gd.x(params.Nx,hill_nz));
            end
            plot(hill_x,hill,'k', 'LineWidth', 1)
        end
        % if left side is mapped (or tank is rotated)
        if lratio ~= 0
            if params.ndims == 3
                hill   = squeeze(gd.z(1,1,hill_nz));
                hill_x = squeeze(gd.x(1,1,hill_nz));
            elseif params.ndims == 2
                hill   = squeeze(gd.z(1,hill_nz));
                hill_x = squeeze(gd.x(1,hill_nz));
            end
            plot(hill_x,hill,'k')
        end
    end

    % axis options
    if (plotaxis(2)-plotaxis(1))/(plotaxis(4)-plotaxis(3)) <= 5 || ...
            strcmp(opts.axisstyle, 'image')
        axis image
    else
        axis(opts.axisstyle)
    end
    axis(plotaxis)
    set(gca,'layer','top')
    box on

    % drawnow if plotting multiple outputs
    if length(t_index) > 1, drawnow, end

    % save figure
    if opts.savefig == true
        direcs = strsplit(opts.dir, '/');
        strt_dir  = pwd;
        for jj = 1:length(direcs)
            if ~(exist(direcs{jj},'dir') == 7)
                mkdir(direcs{jj})
            end
            cd(direcs{jj})
        end
        savefig(gcf,[filename,'_',int2str(ii)]);
        cd(strt_dir)
    end
    hold off
end

% add info into opts
opts.var = var;
opts.output = ii;
opts.time = time;
opts.axis = plotaxis;
opts.clim = clim;
opts.cmap = cmap;
opts.filename = filename;
opts.fig_hand = fighand;
opts.p_hand = p_hand;
if ~strcmpi(opts.var2, 'None') && ~strcmp(var, var2)
    opts.p2_hand = p2_hand;
end
opts.nx = nx;
opts.ny = ny;
opts.nz = nz;
% output plotted data
if strcmpi(params.mapped_grid, 'false') && ~plotting_streamlines1
    data1 = data1';
    if ~strcmpi(opts.var2, 'None')
        data2 = data2';
    end
end
pltinfo.xvar = xvar;
pltinfo.yvar = yvar;
pltinfo.data1 = data1;
if ~strcmpi(opts.var2, 'None')
    pltinfo.data2 = data2;
end
if (plotting_streamlines1 || plotting_streamlines2) && strcmpi(params.mapped_grid, 'true')
    pltinfo.xvar_slice = xvar_slice;
    pltinfo.yvar_slice = yvar_slice;
end
pltinfo.opts = opts;
