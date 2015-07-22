function pltinfo = spins_plot2dunmapped(var, t_index, varargin)
% SPINS_PLOT2DUNMAPPED  Plot cross-section, averages and standard deviations of variables from an unmapped grid.
%
%   pltinfo = spinsplot_unmapped(var,t_i) plots var at t_i
%   other options
%
%    David Deepwell, 2015
global gdpar

    %% get grid, parameters and set plotting options
gd = gdpar.gd;
params = gdpar.params;
spins_plotoptions

    %% figure visibility options
if p.Results.visible == false
    figure('Visible','off'), clf
else
    figure(p.Results.fnum), clf
end

for ii = t_index
    clf
    hold on
    % Title and axis labels
    if strncmp(var,'Mean',4) || strncmp(var,'SD',2)
        if isfield(params, 'plot_interval')
            title(['Spanwise ',var,', t=',int2str(ii*params.plot_interval)])
        else
            title(['Spanwise ',var,', tn=',int2str(ii)])
        end
    elseif params.ndims == 2
        if isfield(params, 'plot_interval')
            title([var,', t=',int2str(ii*params.plot_interval),' s'])
        else
            title([var,', tn=',int2str(ii)])
        end
    else
        if isfield(params, 'plot_interval')
            title([var,', ',p.Results.dimen,'=',num2str(cross_section),...
                   ' m, t=',int2str(ii*params.plot_interval),' s'])
        else
            title([var,', ',p.Results.dimen,'=',num2str(cross_section),...
                   ' m, tn=',int2str(ii)])
        end
    end
    if length(nx)>1, xlabel('x (m)'), else xlabel('y (m)'), end
    if length(nz)>1, ylabel('z (m)'), else ylabel('y (m)'), end

    % get data to plot
    [data1,primcol,cmap] = spins_readdata(var,ii,nx,ny,nz);

    % choose plotting style (contourf may take up less memory,
    % but can be slower than pcolor)
    if strcmpi(var,'Streamline')
        prompt = 'Provide a sensible wave speed in m/s: ';
        uwave = input(prompt);
        disp(['background speed = ',num2str(uwave),' m/s'])
        u1 = data1(:,:,1) - uwave;
        u2 = data1(:,:,2);
        streamslice(xvar,yvar,u1',u2',0.75)
        cont2col = 'r-';
    elseif strcmpi(p.Results.style,'pcolor')
        pcolor(xvar,yvar,data1')
        cont2col = 'w-';
    elseif strcmpi(p.Results.style,'contourf')
        contourf(xvar,yvar,data1',p.Results.ncontourf)
        cont2col = 'w-';
    elseif strcmpi(p.Results.style,'contour')
        contour(xvar,yvar,data1',p.Results.ncontour)
        cont2col = 'k-';
    end
    % add extra information
    shading flat
    colormap(cmap)
    caxis(primcol);
    if p.Results.colorbar == true && ~strcmpi(var,'Streamline')
        colorbar
    end
    if ~strcmp(p.Results.cont2,'None')      % add contours of another field
        if strncmp(var,'Mean',4) || strncmp(var,'SD',2)
            % choose Mean of field if primary field is Mean or SD
            cont2 = ['Mean ',p.Results.cont2];
        else
            cont2 = p.Results.cont2;
        end
        if strcmp(cont2,var)            % read in data only if the field is different
            data2 = data1;
        else
            [data2,~,~] = spins_readdata(cont2,ii,nx,ny,nz);
            contour(xvar,yvar,data2',p.Results.ncont2,cont2col)
        end
    end
    % axis options
    if primaxis(2)/primaxis(4)>6 && (primaxis(2)-primaxis(1)>0.75)
        axis normal
    else
        axis image
    end
    axis(primaxis)
    set(gca,'layer','top')

    % drawnow if plotting multiple outputs
    if length(t_index) > 1, drawnow, end

    % save figure
    if p.Results.savefig == true
        if ~(exist('figures','dir') == 7)
            mkdir figures
        end
        cd figures
        saveas(gcf,[filename,'.fig'],'fig');
        cd('..')
    end
    hold off
end

% output plotted data
pltinfo.xvar = xvar;
pltinfo.yvar = yvar;
pltinfo.data1 = data1;
pltinfo.var1 = var;
try
    pltinfo.data2 = data2;
    pltinfo.var2 = cont2;
end
pltinfo.dimen = p.Results.dimen;
pltinfo.slice = cross_section;
