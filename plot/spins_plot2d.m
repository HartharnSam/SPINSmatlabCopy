function pltinfo = spins_plot2d(var, t_index, varargin)
%  SPINS_PLOT2D  Plot cross-sections, averages or standard deviations of variables.
%
%  Usage:
%	pltinfo = spins_plot2d(var, t_i) plots var at t_i
%	pltinfo = spins_plot2d(var, t_i, 'opt1', val1, ...) plots var at t_i with option 'opt1' as val1
%   
%  Inputs:
%    'var' may be of different forms:
%	any field in the working directory ('rho','u',...)
%	'Density'	   uses the 'rho_reader' file (other labeled fields exist too, ex. Salt, ...)
%	'KE' 		   finds the local kinetic energy
%	'Mean ...'	   takes the spanwise mean of ...
%	'SD ...'	   takes the spanwise standard deviation of ...
%	'Scaled SD ...'    scales SD ... by the maximum of ...
%	'Streamline'	   plots streamlines in the x-z plane
%
%    't_index' may be:
%	an integer for a particular output
%	a vector of outputs (ie. 0:10) will plot each output successively in the same figure
%
%    Optional arguments:
%	Name:	Options			- Description (defaults are in spins_plotoptions.m)
%	---------------------------------------------------------
%	dimen:	{'X','Y','Z'}		- dimension to take cross-section
%	slice:	{double}		- location to take cross-section
%	style:  {'pcolor','contourf','contour'}
%		- type of plot
%	axis:	{[x1 x2 z1 z2]}		- domain to plot
%	xskp:	{integer}		- x-grid points to skip in plot
%	yskp:	{integer}		- y-grid     "
%	zskp:	{integer}		- z-grid     "
%	fnum:	{integer}		- figure window to make plot
%	colorbar:  {boolean}		- plot colorbar?
%	savefig:   {boolean}		- save figure in figure file?
%	visible:   {boolean}		- make figure visible?
%	ncontourf: {integer}		- contours to use for contourf plot
%	ncontour:  {integer}		- contours to use for contour plot
%	cont2:	{field name}		- secondary field to plot as contours
%	ncont2:	{integer}		- contours to use for secondary field
%
%  Outputs:
%    'pltinfo' is a structure containing the plotted fields 
%
%  David Deepwell, 2015
global gdpar

% get grid and parameters
gd = gdpar.gd;
params = gdpar.params;

% set plotting options
spins_plotoptions

% figure visibility options
if opts.visible == false
    figure('Visible','off')
else
    figure(opts.fnum)
end

for ii = t_index
    clf
    hold on
    % Title
    plot_title = var;
    if strncmp(var,'Mean',4) || strncmp(var,'SD',2)
        plot_title = ['Spanwise ', plot_title];
    end
    if params.ndims == 3	% add cross-section information
        plot_title = [plot_title,', ',opts.dimen,'=',num2str(opts.slice),' m'];
    end
    if isfield(params, 'plot_interval')	% add time in seconds or output number
        plot_title = [plot_title,', t=',int2str(ii*params.plot_interval),' s'];
    else
        plot_title = [plot_title,', t_n=',int2str(ii)];
    end
    title(plot_title);
    % axis labels
    if strcmp(opts.dimen, 'X')
        xlabel('y (m)'), ylabel('z (m)')
    elseif strcmp(opts.dimen, 'Y')
        xlabel('x (m)'), ylabel('z (m)')
    elseif strcmp(opts.dimen, 'Z')
        xlabel('x (m)'), ylabel('y (m)')
    end

    % get data to plot
    [data1,primcol,cmap] = spins_readdata(var,ii,nx,ny,nz);
    % if mapped grid and taking horizontal opts.slice, then find interpolation
    if strcmp(opts.dimen, 'Z') && strcmp(params.mapped_grid, 'true')
        [xvar, yvar, data1] = get_fixed_z(xvar, yvar, zvar, data1, opts.slice);
    end
    % transpose unmapped data
    if strcmp(params.mapped_grid, 'false')
        data1 = data1';
    end

    % choose plotting style (contourf may take up less memory,
    % but can be slower than pcolor)
    if strcmpi(var,'Streamline')
        if strcmp(params.mapped_grid, 'false')
            warning('Streamline has not been tested for mapped grids.')
        end
        prompt = 'Provide a sensible wave speed in m/s: ';
        uwave = input(prompt);
        disp(['background speed = ',num2str(uwave),' m/s'])
        u1 = data1(:,:,1) - uwave;
        u2 = data1(:,:,2);
        streamslice(xvar,yvar,u1',u2',0.75)
        cont2col = 'r-';
    elseif strcmp(opts.style,'pcolor')
        pcolor(xvar,yvar,data1)
        cont2col = 'w-';
    elseif strcmp(opts.style,'contourf')
        contourf(xvar,yvar,data1,opts.ncontourf)
        cont2col = 'w-';
    elseif strcmp(opts.style,'contour')
        contour(xvar,yvar,data1,opts.ncontour)
        cont2col = 'k-';
    end

    % add extra information
    shading flat	% need to change in version 2015
    colormap(cmap)
    caxis(primcol);
    if opts.colorbar == true
        colorbar
    end

    % add contours of another field
    if ~strcmpi(opts.cont2,'None')
        if strncmp(var,'Mean',4) || strncmp(var,'SD',2)
            % choose Mean of field if primary field is Mean or SD
            cont2 = ['Mean ',opts.cont2];
        else
            cont2 = opts.cont2;
        end
        if strcmp(cont2,var)            % read in data only if the field is different
            data2 = data1;
        else
            [data2,~,~] = spins_readdata(cont2,ii,nx,ny,nz);
            if strcmp(opts.dimen, 'Z') && strcmp(params.mapped_grid, 'true')
                [xvar, yvar, data2] = get_fixed_z(xvar, yvar, zvar, data2, opts.slice);
            end
            if strcmp(params.mapped_grid, 'false')
                data2 = data2';
            end
            contour(xvar,yvar,data2,opts.ncont2,cont2col)
        end
    end

    % add contour of hill if grid is mapped
    if strcmp(params.mapped_grid,'true') && strcmp(opts.style,'contour')
        hill_nx = nx(1):nx(end);
        if params.ndims == 3
            hill   = squeeze(gd.z(hill_nx,1,params.Nz));
            hill_x = squeeze(gd.x(hill_nx,1,params.Nz));
        elseif params.ndims == 2
            hill   = squeeze(gd.z(hill_nx,params.Nz));
            hill_x = squeeze(gd.x(hill_nx,params.Nz));
        end
        plot(hill_x,hill,'k')
    end

    % axis options
    if (plotaxis(2)-plotaxis(1))/(plotaxis(4)-plotaxis(3)) > 5
        axis normal
    else
        axis image
    end
    axis(plotaxis)
    set(gca,'layer','top')

    % drawnow if plotting multiple outputs
    if length(t_index) > 1, drawnow, end

    % save figure
    if opts.savefig == true
        if ~(exist('figures','dir') == 7)
            mkdir figures
        end
        cd figures
        saveas(gcf,[filename,'.fig'],'fig');
        cd('..')
    end
    hold off

    % output plotted data
    pltinfo.xvar = xvar;
    pltinfo.yvar = yvar;
    pltinfo.data1 = data1;
    pltinfo.var1 = var;
    try
        pltinfo.data2 = data2;
        pltinfo.var2 = cont2;
    end
    pltinfo.dimen = opts.dimen;
    pltinfo.slice = opts.slice;

end
