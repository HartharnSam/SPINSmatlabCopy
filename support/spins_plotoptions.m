% spinsplotoptions.m creates the optional arguments for spinsplot

% shorten some parameters
if isfield(gd, 'x')
    x = gd.x;
    Nx = params.Nx;
    xlimits = params.xlim;
end
if isfield(gd, 'y')
    y = gd.y;
    Ny = params.Ny;
    ylimits = params.ylim;
end
if isfield(gd, 'z')
    z = gd.z;
    Nz = params.Nz;
    zlimits = params.zlim;
end


% define expected options 
exp_dimen = {'X','Y','Z'};
exp_style = {'pcolor','contourf','contour'};
% define defaults 
d.dimen = 'Y';			% dimension
d.slice = [0 0];		% cross-section. opt. arg. must be single number
d.fnum = 1;			% figure window number to use
d.savefig = false;		% save figure? (bool)
d.style = 'contourf';		% plotting style
d.ncontourf = 30;		% plotting regions in contourf style
d.ncontour = 10;		% contours in contour style
d.cont2 = 'Density';		% secondary field to plot
d.ncont2 = 6;			% contours of secondary field
d.colorbar = true;		% colorbar? (bool)
d.xskp = 1;			% x-grid points to skip
d.yskp = 1;			% y	"
d.zskp = 1;			% z	"
d.axis = 0;			% axis to plot. 0 denotes use of full domain
d.visible = true;		% make plot visible or not (bool)

% parse options
p = inputParser;
addParameter(p,'dimen',d.dimen, @(x) any(validatestring(x,exp_dimen)))
addParameter(p,'slice',d.slice,@isnumeric)
addParameter(p,'fnum',d.fnum,@isnumeric)
addParameter(p,'savefig',d.savefig,@islogical)
addParameter(p,'style',d.style, @(x) any(validatestring(x,exp_style)))
addParameter(p,'ncontourf',d.ncontourf,@isnumeric)
addParameter(p,'ncontour',d.ncontour,@isnumeric)
addParameter(p,'cont2',d.cont2,@ischar)
addParameter(p,'ncont2',d.ncont2,@isnumeric)
addParameter(p,'colorbar',d.colorbar,@islogical)
addParameter(p,'xskp',d.xskp,@isnumeric)
addParameter(p,'yskp',d.yskp,@isnumeric)
addParameter(p,'zskp',d.zskp,@isnumeric)
addParameter(p,'axis',d.axis,@isnumeric)
addParameter(p,'visible',d.visible,@islogical)
try
    parse(p,varargin{:})
catch
    parse(p,varargin{1}{:})
end

% choose default cross-section slice based on which dimension is plotted
if length(p.Results.slice) == 2
    if strcmp(p.Results.dimen, 'X')
        cross_section = sum(xlimits)/2;
    elseif strcmp(p.Results.dimen, 'Y')
        cross_section = sum(ylimits)/2;
    elseif strcmp(p.Results.dimen, 'Z')
        cross_section = sum(zlimits)/2;
    end
elseif length(p.Results.slice) == 1
    cross_section = p.Results.slice;
end

% find grid points to read in
if strcmp(p.Results.dimen,'X')		% X dimen
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
    if params.ndims == 3
        nx = nearestindex(x, cross_section);
    else
        nx = 1;
    end
    xvar = y(ny);
    yvar = z(nz);
elseif strcmp(p.Results.dimen,'Y')	% Y dimen
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
    if params.ndims == 3
        ny = nearestindex(y, cross_section);
    else
        ny = 1;
    end
    xvar = x(nx);
    yvar = z(nz);
elseif strcmp(p.Results.dimen,'Z')	% Z dimen
    if length(p.Results.axis) == 1		% no axis flag
        nx = 1:p.Results.xskp:Nx;
        ny = 1:p.Results.yskp:Ny;
        primaxis=[xlimits ylimits]; % the plot area
    else					% with axis flag
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
    if params.ndims == 3
        nz = nearestindex(z, cross_section);
    else
        nz = 1;
    end
    xvar = x(nx);
    yvar = y(ny);
end
