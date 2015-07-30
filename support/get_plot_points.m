function [nx, ny, nz, xvar, yvar, zvar, plotaxis] = get_plot_points(gd, params, cross_section, opts)
% GET_PLOT_POINTS  get the points to plot given the options (p) and the parameters (params)
%
%   used exclusively in spins_plotoptions
%
%   David Deepwell, 2015

    if strcmp(params.mapped_grid,'false')
        [nx, ny, nz, xvar, yvar, zvar, plotaxis] = unmapped_points(gd, params, cross_section, opts);
    elseif strcmp(params.mapped_grid,'true')
        [nx, ny, nz, xvar, yvar, zvar, plotaxis] = mapped_points(gd, params, cross_section, opts);
    end
end

function [nx, ny, nz, xvar, yvar, zvar, plotaxis] = unmapped_points(gd, params, cross_section, opts)

    % if full grid, give vector grid
    gdnames = fieldnames(gd);
    if isvector(gd.(gdnames{1}))
        gdvec = gd;
    else
        gdvec = get_vector_grid(gd, params);
    end

    % shorten some parameters
    if isfield(gd, 'x')
        x = gdvec.x;
        Nx = params.Nx;
        xlimits = params.xlim;
    end
    if isfield(gd, 'y')
        y = gdvec.y;
        Ny = params.Ny;
        ylimits = params.ylim;
    end
    if isfield(gd, 'z')
        z = gdvec.z;
        Nz = params.Nz;
        zlimits = params.zlim;
    end

    % find grid points to read in
    if strcmp(opts.dimen,'X')		% X dimen
        if params.ndims == 3
            nx = nearestindex(x, cross_section);
        else
            nx = 1;
        end
        if length(opts.axis) == 1		% no axis flag
            ny = 1:opts.yskp:Ny;
            nz = 1:opts.zskp:Nz;
            plotaxis=[ylimits zlimits]; % the plot area
        else					% with axis flag
            plotaxis = opts.axis;
            if plotaxis(2)<= plotaxis(1) || plotaxis(4)<=plotaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(y, plotaxis(1));
            xvarR = nearestindex(y, plotaxis(2));
            yvarB = nearestindex(z, plotaxis(3));
            yvarT = nearestindex(z, plotaxis(4));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            ny = xvarL2:opts.xskp:xvarR2;
            nz = yvarB2:opts.zskp:yvarT2;
        end
        xvar = y(ny);
        yvar = z(nz);
    elseif strcmp(opts.dimen,'Y')	% Y dimen
        if params.ndims == 3
            ny = nearestindex(y, cross_section);
        else
            ny = 1;
        end
        if length(opts.axis) == 1		% no axis flag
            nx = 1:opts.xskp:Nx;
            nz = 1:opts.zskp:Nz;
            plotaxis = [xlimits zlimits]; % the plot area
        else					% with axis flag
            plotaxis = opts.axis;
            if plotaxis(2)<= plotaxis(1) || plotaxis(4)<=plotaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(x, plotaxis(1));
            xvarR = nearestindex(x, plotaxis(2));
            yvarB = nearestindex(z, plotaxis(3));
            yvarT = nearestindex(z, plotaxis(4));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            nx = xvarL2:opts.xskp:xvarR2;
            nz = yvarB2:opts.zskp:yvarT2;
        end
        xvar = x(nx);
        yvar = z(nz);
    elseif strcmp(opts.dimen,'Z')      % Z dimen
        if params.ndims == 3
            nz = nearestindex(z, cross_section);
        else
            nz = 1;
        end
        if length(opts.axis) == 1              % no axis flag
            nx = 1:opts.xskp:Nx;
            ny = 1:opts.yskp:Ny;
            plotaxis=[xlimits ylimits]; % the plot area
        else                                        % with axis flag
            plotaxis = opts.axis;
            if plotaxis(2)<= plotaxis(1) || plotaxis(4)<=plotaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(x, plotaxis(1));
            xvarR = nearestindex(x, plotaxis(2));
            yvarB = nearestindex(y, plotaxis(3));
            yvarT = nearestindex(y, plotaxis(4));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            nx = xvarL2:opts.xskp:xvarR2;
            ny = yvarB2:opts.zskp:yvarT2;
        end
        xvar = x(nx);
        yvar = y(ny);
    end

    zvar = [];	% empty arry, since it's not used for unmapped grids
end


function [nx, ny, nz, xvar, yvar, zvar, plotaxis] = mapped_points(gd, params, cross_section, opts)

    % if vector grid, call for vector grid
    gdnames = fieldnames(gd);
    if isvector(gd.(gdnames{1}))
        error('Mapped grid requires full grid to be read in. Use spins_gridparams("Full").')
    end
    % read in vector grid for easiness with other dimensions
    gdvec = get_vector_grid(gd, params);

    % shorten some parameters
    if isfield(gd, 'x')
        x = gd.x;
        x1d = gdvec.x;
        Nx = params.Nx;
        xlimits = params.xlim;
    end
    if isfield(gd, 'y')
        y = gd.y;
        y1d = gdvec.y;
        Ny = params.Ny;
        ylimits = params.ylim;
    end
    if isfield(gd, 'z')
        z = gd.z;
        z1d = gdvec.z;
        Nz = params.Nz;
        zlimits = params.zlim;
    end

    % find grid points to read in
    if strcmp(opts.dimen,'X')		% X dimen
        if params.ndims == 3
            nx = nearestindex(x1d, cross_section);
        else
            error('get_plot_points assumes x-z plane for 2D mapped grids.')
        end
        if length(opts.axis) == 1		% no axis flag
            ny = 1:opts.yskp:Ny;
            nz = 1:opts.zskp:Nz;
            plotaxis=[ylimits zlimits]; % the plot area
        else					% with axis flag
            plotaxis = opts.axis;
            if plotaxis(2)<= plotaxis(1) || plotaxis(4)<=plotaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(y1d, plotaxis(1));
            xvarR = nearestindex(y1d, plotaxis(2));
            yvarB = nearestindex(z(nx,1,:), plotaxis(3));
            yvarT = nearestindex(z(nx,1,:), plotaxis(4));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            ny = xvarL2:opts.xskp:xvarR2;
            nz = yvarB2:opts.zskp:yvarT2;
        end
        xvar = squeeze(y(nx,ny,nz));
        yvar = squeeze(z(nx,ny,nz));
        zvar = [];
    elseif strcmp(opts.dimen,'Y')	% Y dimen
        if params.ndims == 3
            ny = nearestindex(y1d, cross_section);
        else
            ny = 1;
        end
        if length(opts.axis) == 1		% no axis flag
            nx = 1:opts.xskp:Nx;
            nz = 1:opts.zskp:Nz;
            plotaxis = [xlimits zlimits]; % the plot area
        else					% with axis flag
            plotaxis = opts.axis;
            if plotaxis(2)<= plotaxis(1) || plotaxis(4)<=plotaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(x1d, plotaxis(1));
            xvarR = nearestindex(x1d, plotaxis(2));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            if params.ndims == 3
                yvarB = max(nearestindex(squeeze(z(xvarL2:xvarR2,1,:))', plotaxis(3)));
                yvarT = min(nearestindex(squeeze(z(xvarL2:xvarR2,1,:))', plotaxis(4)));
            elseif params.ndims == 2
                yvarB = max(nearestindex(z(xvarL2:xvarR2,:)', plotaxis(3)));
                yvarT = min(nearestindex(z(xvarL2:xvarR2,:)', plotaxis(4)));
            end
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            nx = xvarL2:opts.xskp:xvarR2;
            nz = yvarB2:opts.zskp:yvarT2;
        end
        if params.ndims == 3
            xvar = squeeze(x(nx,ny,nz));
            yvar = squeeze(z(nx,ny,nz));
        elseif params.ndims == 2
            xvar = x(nx,nz);
            yvar = z(nx,nz);
        end
        zvar = [];
    elseif strcmp(opts.dimen,'Z')      % Z dimen
        if params.ndims == 3
            nz = 1:Nz;
        else
            error('get_plot_points assumes x-z plane for 2D mapped grids.')
        end
        if length(opts.axis) == 1              % no axis flag
            nx = 1:opts.xskp:Nx;
            ny = 1:opts.yskp:Ny;
            plotaxis=[xlimits ylimits]; % the plot area
        else                                        % with axis flag
            plotaxis = opts.axis;
            if plotaxis(2)<= plotaxis(1) || plotaxis(4)<=plotaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(x1d, plotaxis(1));
            xvarR = nearestindex(x1d, plotaxis(2));
            yvarB = nearestindex(y1d, plotaxis(3));
            yvarT = nearestindex(y1d, plotaxis(4));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            nx = xvarL2:opts.xskp:xvarR2;
            ny = yvarB2:opts.zskp:yvarT2;
        end
        xvar = x(nx,ny,nz);
        yvar = y(nx,ny,nz);
        zvar = z(nx,ny,nz);
    end
end
