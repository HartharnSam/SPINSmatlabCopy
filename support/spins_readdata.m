function [data, primcol, cmap] = spins_readdata(var, ii, nx, ny, nz)
%  SPINS_READDATA  read in data for the given variable
%
%  Usage:
%    [data,primcol,cmap] = spinsread(var, t_i, nx, ny, nz) reads in var at positions (nx,ny,nz) at t_i
%                          where nx, ny, nz are scalars or 1D vectors
%
%  Inputs:
%    'var' may be of different forms:
%	any field in the working directory ('rho','u',...)
%	'Density'          uses the 'rho_reader' file (other labeled fields exist too, ex. Salt, ...)
%	'KE'               finds the local kinetic energy
%	'Mean ...'         takes the spanwise mean of ...
%	'SD ...'           takes the spanwise standard deviation of ...
%	'Scaled SD ...'    scales SD ... by the maximum of ...
%	'Streamline'       plots streamlines in the x-z plane
%
%  Outputs:
%    data	- field you want
%    primcol	- extent of the colorbar
%    cmap	- colormap for the given field
%
%  David Deepwell, 2015
global gdpar

% get grid and parameters
gd = gdpar.gd;
params = gdpar.params;

% shorten some parameters
if isfield(gd, 'y')
    Ny = params.Ny;
end

% default colormap
cmap = 'darkjet';

% read data
if params.ndims == 3		% for 3D data
    % find the variable to plot (ie. remove the mean or SD from name)
    if strncmp(var,'Mean',4) || strncmp(var,'SD',2) || strncmp(var,'Scaled SD',9)
        ny = 1:Ny;
        if strncmp(var,'Mean',4)
            varorig = var;
            var = strsplit(var,'Mean ');
            var = var{end};
        elseif strncmp(var,'SD',2)
            varorig = var;
            var = strsplit(var,'SD ');
            var = var{end};
        elseif strncmp(var,'Scaled SD',9)
            varorig = var;
            var = strsplit(var,'Scaled SD ');
            var = var{end};
        end
    end
    % choose which dye field to read
    if ~isempty(strfind(var,'Dye'));		primcol = [-1 1];
        if ~isempty(strfind(var,'1b'))
            data = dye1b_reader(ii,nx,ny,nz);
        elseif ~isempty(strfind(var,'1'))
            data = dye1_reader(ii,nx,ny,nz);
        elseif ~isempty(strfind(var,'2'))
            data = dye2_reader(ii,nx,ny,nz);
        elseif ~isempty(strfind(var,'3'))
            data = dye3_reader(ii,nx,ny,nz);
        else
            try
                data = dye_reader(ii,nx,ny,nz);
            catch
                data = tracer_reader(ii,nx,ny,nz);
            end
        end
    % read in salt field
    elseif strcmp(var,'Salt') || strcmp(var,'S')
        try
            data = s_reader(ii,nx,ny,nz);	primcol = [min(data(:)) max(data(:))];
        catch    
            data = S_reader(ii,nx,ny,nz);	primcol = [min(data(:)) max(data(:))];
        end
    % read in temperature field
    elseif strcmp(var,'Temperature') || strcmp(var,'T')
        try
            data = t_reader(ii,nx,ny,nz);	primcol = [min(data(:)) max(data(:))];
        catch
            data = T_reader(ii,nx,ny,nz);	primcol = [min(data(:)) max(data(:))];
        end
    % read in density field
    elseif strcmp(var,'Density')
        try
            data = rho_reader(ii,nx,ny,nz);
        catch
            s = s_reader(ii,nx,ny,nz);
            t = t_reader(ii,nx,ny,nz);
            data = eqn_of_state(t,s);
        end
        if min(data(:)) > 0 
            primcol = [min(data(:)) max(data(:))];
        else
            primcol = [-1 1]*max(abs(data(:)));
        end
    elseif strcmp(var,'Pressure') || strcmp(var,'P')
        data = p_reader(ii,nx,ny,nz);		primcol = [-1 1]*max(abs(data(:)));
    elseif strcmp(var,'U')
        data = u_reader(ii,nx,ny,nz);		primcol = [-1 1]*max(abs(data(:)));
    elseif strcmp(var,'V')
        data = v_reader(ii,nx,ny,nz);		primcol = [-1 1]*max(abs(data(:)));
    elseif strcmp(var,'W')
        data = w_reader(ii,nx,ny,nz);		primcol = [-1 1]*max(abs(data(:)));
    % read in kinetic energy
    elseif strcmp(var,'KE')
        cmap = 'hot';
        u = u_reader(ii,nx,ny,nz);
        v = v_reader(ii,nx,ny,nz);
        w = w_reader(ii,nx,ny,nz);
        data = 0.5*(u.^2 + v.^2 + w.^2);	primcol = [0 1]*max(data(:));
        clearvars u v w
    elseif strcmp(var,'Vorticity')
        error('Vorticity not written yet.');
    % read in data for plotting streamlines
    elseif strcmp(var,'Streamline')
        ny = 1:Ny;
        u = squeeze(mean(u_reader(ii,nx,ny,nz),2));
        w = squeeze(mean(w_reader(ii,nx,ny,nz),2));
        %u = u_reader(ii,nx,ny,nz);
        %w = w_reader(ii,nx,ny,nz);
        data = ones([size(u) 2]);
        data(:,:,1) = u;
        data(:,:,2) = w;
        primcol = [-1 1];
    % read in data for given file name
    else
        try
            var_reader = str2func([var,'_reader']);
            data = var_reader(ii,nx,ny,nz);	primcol = [-1 1]*max(abs(data(:)));
        catch
            error('Variable not understood or output does not exist.');
        end
    end
    % take mean or standard deviation if asked
    if exist('varorig', 'var')
        if strncmp(varorig,'Mean',4)	 % take mean
            data = squeeze(mean(data,2));
            if ~isempty(strfind(var,'Dye'))  % put dye colorbar between [-1,1]
                primcol = [-1 1];
            elseif ~strcmp(var,'KE')
                primcol = [-1 1]*max(abs(data(:)));
            end
        elseif strncmp(varorig,'SD',2)       % take standard deviation
            data = squeeze(std(data,[],2));
            cmap = 'hot';
            primcol = [0 1]*max(data(:));
        elseif strncmp(varorig,'Scaled SD',9)       % take standard deviation
            datamax = max(data(:));
            data = squeeze(std(data,[],2))/datamax;
            cmap = 'hot';
            primcol = [0 1]*max(data(:));
        end
    end
elseif params.ndims == 2
    % can't accept Mean or SD
    if strncmp(var,'Mean',4) || strncmp(var,'SD',2) || strncmp(var,'Scaled SD',9)
        error('Mean, SD, and Scaled SD are not supported on 2D data.');
    end
    % choose which dye field to read
    if ~isempty(strfind(var,'Dye'));		primcol  =  [-1 1];
        if ~isempty(strfind(var,'1b'))
            data = dye1b_reader(ii,nx,nz);
        elseif ~isempty(strfind(var,'1'))
            data = dye1_reader(ii,nx,nz);
        elseif ~isempty(strfind(var,'2'))
            data = dye2_reader(ii,nx,nz);
        elseif ~isempty(strfind(var,'3'))
            data = dye3_reader(ii,nx,nz);
        else
            try
                data = dye_reader(ii,nx,nz);
            catch
                data = tracer_reader(ii,nx,nz);
            end
        end
    elseif strcmp(var,'Salt') || strcmp(var,'S')
        try
            data = s_reader(ii,nx,nz);		primcol = [min(data(:)) max(data(:))];
        catch
            data = S_reader(ii,nx,nz);		primcol = [min(data(:)) max(data(:))];
        end
    elseif strcmp(var,'Temperature') || strcmp(var,'T')
        try
            data = t_reader(ii,nx,nz);		primcol = [min(data(:)) max(data(:))];
        catch
            data = T_reader(ii,nx,nz);		primcol = [min(data(:)) max(data(:))];
        end
    elseif strcmp(var,'Density')
        try
            data = rho_reader(ii,nx,nz);
        catch
            s = s_reader(ii,nx,nz);
            t = t_reader(ii,nx,nz);
            data = eqn_of_state(t,s);
        end
        if min(data(:)) > 0
            primcol = [min(data(:)) max(data(:))];
        else
            primcol = [-1 1]*max(abs(data(:)));
        end
    elseif strcmp(var,'U')
        data = u_reader(ii,nx,nz);		primcol = [-1 1]*max(abs(data(:)));
    elseif strcmp(var,'W')
        data = w_reader(ii,nx,nz);		primcol = [-1 1]*max(abs(data(:)));
    elseif strcmp(var,'KE');
        data = u_reader(ii,nx,nz).^2 + w_reader(ii,nx,nz).^2;	  primcol = [0 1]*max(data(:));
        cmap = 'hot';
    elseif strcmp(var,'Vorticity')
        error('Vorticity not written yet.');
    elseif strcmp(var,'Streamline')
        u = u_reader(ii,nx,nz);
        w = w_reader(ii,nx,nz);
        data = ones([size(u) 2]);
        data(:,:,1) = u;
        data(:,:,2) = w;
        primcol = [-1 1];
    else
        try
            var_reader = str2func([var,'_reader']);
            data = var_reader(ii,nx,nz);	primcol = [-1 1]*max(abs(data(:)));
        catch
            error('Variable not understood or output does not exist.');
        end
    end
end
