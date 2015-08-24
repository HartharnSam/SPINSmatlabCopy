function data = spins_readdata(var, ii, nx, ny, nz)
%  SPINS_READDATA  read in data for the given variable
%
%  Usage:
%    data = spins_readdata(var, t_i, nx, ny, nz) reads in var at positions (nx,ny,nz) at t_i
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
%    data	- field you want, just where you want it
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

% read data
if params.ndims == 3		% for 3D data
    % Remove prefix from var (ie. remove the mean or SD from name)
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
    if ~isempty(strfind(var,'Dye'));
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
    elseif strcmp(var,'Salt') || strcmp(var,'S')
        try
            data = s_reader(ii,nx,ny,nz);
        catch    
            data = S_reader(ii,nx,ny,nz);
        end
    elseif strcmp(var,'Temperature') || strcmp(var,'T')
        try
            data = t_reader(ii,nx,ny,nz);
        catch
            data = T_reader(ii,nx,ny,nz);
        end
    elseif strcmp(var,'Density')
        try
            data = rho_reader(ii,nx,ny,nz);
        catch
            try
                s = s_reader(ii,nx,ny,nz);
            catch
                s = S_reader(ii,nx,ny,nz);
            end
            try
                t = t_reader(ii,nx,ny,nz);
            catch
                t = T_reader(ii,nx,ny,nz);
            end
            data = eqn_of_state(t, s);
        end
    elseif strcmp(var,'Pressure') || strcmp(var,'P')
        data = p_reader(ii,nx,ny,nz);
    elseif strcmp(var,'U')
        data = u_reader(ii,nx,ny,nz);
    elseif strcmp(var,'V')
        data = v_reader(ii,nx,ny,nz);
    elseif strcmp(var,'W')
        data = w_reader(ii,nx,ny,nz);
    % read in kinetic energy
    elseif strcmp(var,'KE')
        u = u_reader(ii,nx,ny,nz);
        v = v_reader(ii,nx,ny,nz);
        w = w_reader(ii,nx,ny,nz);
        data = 0.5*(u.^2 + v.^2 + w.^2);
        clearvars u v w
    % read in vorticity
    elseif strcmp(var,'Vorticity')
        error('Vorticity not written yet.');
    % read in gradient Richardson number
    elseif strcmp(var, 'Ri')
        if strcmp(params.type_z, 'NO_SLIP') || strcmp(params.type_z, 'CHEBY')
            error('Gradient Richardson number for Chebyshev grids are not possible yet.')
        elseif length(ny) > 1
            error('Gradient Richardson number must be plotted in x-z plane.')
        else
            % read in data
            rho = rho_reader(ii,nx,ny,nz)';
            u   = u_reader(ii,nx,ny,nz)';
            g = params.g;
            rho_0 = params.rho_0;
            % average nearby data
            navg = 2; %radius to average over
            filter = fspecial('disk', navg);
            rho = imfilter(rho, filter, 'replicate');
            u = imfilter(u, filter, 'replicate');
            % Diff matrix and derivatives
            Dz = FiniteDiff(gd.z(nz), 1, 2);
            Uz_sq = (Dz*u).^2;
            N_sq  = -g/rho_0*Dz*rho;
            data = (N_sq./Uz_sq)';
        end
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
    % read in data for given file name
    else
        try
            var_reader = str2func([var,'_reader']);
            data = var_reader(ii,nx,ny,nz);
        catch
            error('Variable not understood or output does not exist.');
        end
    end
    % take mean or standard deviation if asked
    if exist('varorig', 'var')
        if strncmp(varorig,'Mean',4)	 % take mean
            data = squeeze(mean(data,2));
        elseif strncmp(varorig,'SD',2)       % take standard deviation
            data = squeeze(std(data,[],2));
        elseif strncmp(varorig,'Scaled SD',9)       % take standard deviation
            datamax = max(data(:));
            data = squeeze(std(data,[],2))/datamax;
        end
    end
elseif params.ndims == 2
    % can't accept Mean or SD
    if strncmp(var,'Mean',4) || strncmp(var,'SD',2) || strncmp(var,'Scaled SD',9)
        error('Mean, SD, and Scaled SD are not supported on 2D data.');
    end
    % choose which dye field to read
    if ~isempty(strfind(var,'Dye'));
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
            data = s_reader(ii,nx,nz);
        catch
            data = S_reader(ii,nx,nz);
        end
    elseif strcmp(var,'Temperature') || strcmp(var,'T')
        try
            data = t_reader(ii,nx,nz);
        catch
            data = T_reader(ii,nx,nz);
        end
    elseif strcmp(var,'Density')
        try
            data = rho_reader(ii,nx,nz);
        catch
            try
                s = s_reader(ii,nx,nz);
            catch
                s = S_reader(ii,nx,nz);
            end
            try
                t = t_reader(ii,nx,nz);
            catch
                t = T_reader(ii,nx,nz);
            end
            data = eqn_of_state(t,s);
        end
    elseif strcmp(var,'U')
        data = u_reader(ii,nx,nz);
    elseif strcmp(var,'W')
        data = w_reader(ii,nx,nz);
    elseif strcmp(var,'KE');
        data = u_reader(ii,nx,nz).^2 + w_reader(ii,nx,nz).^2;
    elseif strcmp(var,'Vorticity')
        error('Vorticity not written yet.');
    elseif strcmp(var,'Streamline')
        u = u_reader(ii,nx,nz);
        w = w_reader(ii,nx,nz);
        data = ones([size(u) 2]);
        data(:,:,1) = u;
        data(:,:,2) = w;
    else
        try
            var_reader = str2func([var,'_reader']);
            data = var_reader(ii,nx,nz);
        catch
            error('Variable not understood or output does not exist.');
        end
    end
end
