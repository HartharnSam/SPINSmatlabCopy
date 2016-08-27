function [nx, ny, nz, xvar, yvar, zvar, plotaxis] = get_plot_points(opts)
%  GET_PLOT_POINTS  get the indices, grid, and axis domain for making
%                   a 2D cross-section plot given the prescribed options (opts)
%
%  Usage:
%    [nx, ny, nz, xvar, yvar, zvar, plotaxis] = get_plot_points(opts)
%
%  Inputs:
%    opts is a structure with the following fields:
%       Name:   Options                 - Description (defaults are in spins_plotoptions.m)
%       ---------------------------------------------------------
%       dimen:  {'X','Y','Z'}           - dimension to take cross-section
%       slice:  {double}                - location to take cross-section
%       axis:   {[x1 x2 z1 z2]}         - domain to plot
%       xskp:   {integer}               - x-grid points to skip in plot
%       yskp:   {integer}               - y-grid     "
%       zskp:   {integer}               - z-grid     "
%
%  Outputs:
%    nx, ny, nz		- indices along each dimension which will be used/read
%    xvar, yvar		- grid vectors (unmapped), or matrices (mapped) for the 2D plot
%    zvar		- only used for mapped grids when the plot will be in the x-y plane.
%			  Interpolation between the changing depth requires this to be read
%    plotaxis		- domain of the 2D plot
%
%  David Deepwell, 2015
global gdpar

    % get grid and parameters
    gd = gdpar.gd;
    params = gdpar.params;

    % Get points based upon grid type
    if strcmp(params.mapped_grid, 'false')
        [nx, ny, nz, xvar, yvar, zvar, plotaxis] = unmapped_points(gd, params, opts);
    elseif strcmp(params.mapped_grid, 'true')
        [nx, ny, nz, xvar, yvar, zvar, plotaxis] = mapped_points(gd, params, opts);
    end
end

function [nx, ny, nz, xvar, yvar, zvar, plotaxis] = unmapped_points(gd, params, opts)
% get points for an unmapped grid

    % if full grid, give vector grid
    gdnames = fieldnames(gd);
    if isvector(gd.(gdnames{1}))
        gdvec = gd;
    else
        gdvec = get_vector_grid(gd);
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
    if strcmpi(opts.dimen,'X')		% X dimen
        if params.ndims == 3
            nx = nearestindex(x, opts.slice);
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
    elseif strcmpi(opts.dimen,'Y')	% Y dimen
        if params.ndims == 3
            ny = nearestindex(y, opts.slice);
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
    elseif strcmpi(opts.dimen,'Z')      % Z dimen
        if params.ndims == 3
            nz = nearestindex(z, opts.slice);
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


function [nx, ny, nz, xvar, yvar, zvar, plotaxis] = mapped_points(gd, params, opts)
% get points for an unmapped grid

    % if vector grid, call for vector grid
    gdnames = fieldnames(gd);
    if isvector(gd.(gdnames{1}))
        error('A mapped grid requires full grid to be read in. Use spins_gridparams("Full").')
    end
    % read in vector grid for easiness with other dimensions
    gdvec = get_vector_grid(gd);

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
    if strcmpi(opts.dimen,'X')		% X dimen
        if params.ndims == 3
            nx = nearestindex(x1d, opts.slice);
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
    elseif strcmpi(opts.dimen,'Y')	% Y dimen
        if params.ndims == 3
            ny = nearestindex(y1d, opts.slice);
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
    elseif strcmpi(opts.dimen,'Z')      % Z dimen
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
