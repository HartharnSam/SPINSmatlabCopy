
[nx, ny, nz, xvar, yvar] = get_plot_grid(gd, params, cross_section, p, primaxis)

% take full grid and give vector grid
gdnames = fieldnames(gd);
if isvector(gd.(gdnames{1}))
    gdvec = gd
else
    gdvec = get_vector_grid(gd, params)
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
    if strcmp(params.mapped_grid,'true')
        z = gd.z(:,1,:);
    else
        z = gdvec.z;
    end
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
        if isvector(z)
            yvarB = nearestindex(z, primaxis(3));
