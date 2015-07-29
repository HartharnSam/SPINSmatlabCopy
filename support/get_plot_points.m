function [nx, ny, nz, xvar, yvar, primaxis] = get_plot_points(gd, params, cross_section, p)
% GET_PLOT_POINTS  get the points to plot given the options (p) and the parameters (params)
%
%   used exclusively in spins_plotoptions
%
%   David Deepwell, 2015

    if strcmp(params.mapped_grid,'false')
        [nx, ny, nz, xvar, yvar, primaxis] = unmapped_points(gd, params, cross_section, p);
    elseif strcmp(params.mapped_grid,'true')
        [nx, ny, nz, xvar, yvar, primaxis] = mapped_points(gd, params, cross_section, p);
    end
end

function [nx, ny, nz, xvar, yvar, primaxis] = unmapped_points(gd, params, cross_section, p)

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
    if strcmp(p.Results.dimen,'X')		% X dimen
        if params.ndims == 3
            nx = nearestindex(x, cross_section);
        else
            nx = 1;
        end
        if length(p.Results.axis) == 1		% no axis flag
            ny = 1:p.Results.yskp:Ny;
            nz = 1:p.Results.zskp:Nz;
            primaxis=[ylimits zlimits]; % the plot area
        else					% with axis flag
            primaxis = p.Results.axis;
            if primaxis(2)<= primaxis(1) || primaxis(4)<=primaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(y, primaxis(1));
            xvarR = nearestindex(y, primaxis(2));
            yvarB = nearestindex(z, primaxis(3));
            yvarT = nearestindex(z, primaxis(4));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            ny = xvarL2:p.Results.xskp:xvarR2;
            nz = yvarB2:p.Results.zskp:yvarT2;
        end
        xvar = y(ny);
        yvar = z(nz);
    elseif strcmp(p.Results.dimen,'Y')	% Y dimen
        if params.ndims == 3
            ny = nearestindex(y, cross_section);
        else
            ny = 1;
        end
        if length(p.Results.axis) == 1		% no axis flag
            nx = 1:p.Results.xskp:Nx;
            nz = 1:p.Results.zskp:Nz;
            primaxis = [xlimits zlimits]; % the plot area
        else					% with axis flag
            primaxis = p.Results.axis;
            if primaxis(2)<= primaxis(1) || primaxis(4)<=primaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(x, primaxis(1));
            xvarR = nearestindex(x, primaxis(2));
            yvarB = nearestindex(z, primaxis(3));
            yvarT = nearestindex(z, primaxis(4));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            nx = xvarL2:p.Results.xskp:xvarR2;
            nz = yvarB2:p.Results.zskp:yvarT2;
        end
        xvar = x(nx);
        yvar = z(nz);
    elseif strcmp(p.Results.dimen,'Z')      % Z dimen
        if params.ndims == 3
            nz = nearestindex(z, cross_section);
        else
            nz = 1;
        end
        if length(p.Results.axis) == 1              % no axis flag
            nx = 1:p.Results.xskp:Nx;
            ny = 1:p.Results.yskp:Ny;
            primaxis=[xlimits ylimits]; % the plot area
        else                                        % with axis flag
            primaxis = p.Results.axis;
            if primaxis(2)<= primaxis(1) || primaxis(4)<=primaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(x, primaxis(1));
            xvarR = nearestindex(x, primaxis(2));
            yvarB = nearestindex(y, primaxis(3));
            yvarT = nearestindex(y, primaxis(4));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            nx = xvarL2:p.Results.xskp:xvarR2;
            ny = yvarB2:p.Results.zskp:yvarT2;
        end
        xvar = x(nx);
        yvar = y(ny);
    end
end


function [nx, ny, nz, xvar, yvar, primaxis] = mapped_points(gd, params, cross_section, p)

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
    if strcmp(p.Results.dimen,'X')		% X dimen
        if params.ndims == 3
            nx = nearestindex(x1d, cross_section);
        else
            error('get_plot_points assumes 2D is in x-z plane.')
        end
        if length(p.Results.axis) == 1		% no axis flag
            ny = 1:p.Results.yskp:Ny;
            nz = 1:p.Results.zskp:Nz;
            primaxis=[ylimits zlimits]; % the plot area
        else					% with axis flag
            primaxis = p.Results.axis;
            if primaxis(2)<= primaxis(1) || primaxis(4)<=primaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(y1d, primaxis(1));
            xvarR = nearestindex(y1d, primaxis(2));
            yvarB = nearestindex(z(nx,1,:), primaxis(3));
            yvarT = nearestindex(z(nx,1,:), primaxis(4));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            ny = xvarL2:p.Results.xskp:xvarR2;
            nz = yvarB2:p.Results.zskp:yvarT2;
        end
        xvar = squeeze(y(nx,ny,nz));
        yvar = squeeze(z(nx,ny,nz));
    elseif strcmp(p.Results.dimen,'Y')	% Y dimen
        if params.ndims == 3
            ny = nearestindex(y, cross_section);
        %else
        %    ny = 1;
        end
        if length(p.Results.axis) == 1		% no axis flag
            nx = 1:p.Results.xskp:Nx;
            nz = 1:p.Results.zskp:Nz;
            primaxis = [xlimits zlimits]; % the plot area
        else					% with axis flag
            primaxis = p.Results.axis;
            if primaxis(2)<= primaxis(1) || primaxis(4)<=primaxis(3)
                error('Axis must be ordered correctly.')
            end
            xvarL = nearestindex(x1d, primaxis(1));
            xvarR = nearestindex(x1d, primaxis(2));
            xvarL2 = min(xvarL,xvarR); xvarR2 = max(xvarL,xvarR);
            if params.ndims == 3
                yvarB = nearestindex(z(xvarL2:xvarR2,1,:), primaxis(3));
                yvarT = nearestindex(z(xvarL2:xvarR2,1,:), primaxis(4));
            elseif params.ndims == 2
            end
            yvarB2 = min(yvarB,yvarT); yvarT2 = max(yvarB,yvarT);
            nx = xvarL2:p.Results.xskp:xvarR2;
            nz = yvarB2:p.Results.zskp:yvarT2;
        end
        xvar = x(nx);
        yvar = z(nz);
end
