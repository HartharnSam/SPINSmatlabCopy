function data = spins_readdata(var, ii, nx, ny, nz, dimen)
%  SPINS_READDATA  read in data for the given variable
%
%  Usage:
%    data = spins_readdata(var, t_i, nx, ny, nz) reads in var at positions (nx,ny,nz) at t_i
%                          where nx, ny, nz are scalars or 1D vectors
%
%  Inputs:
%    'var' may be of different forms:
%	any field in the working directory ('rho','u',...)
%	'Density'          searches for rho otherwise searches for t and s
%	'KE'               finds the local kinetic energy
%	'speed'            finds the magnitude of the local velocity vector 
%   'Ri'               gradient Richardson number
%	'Streamline'       plots streamlines in the x-z plane
%	'Mean ...'         takes the spanwise mean of ...
%	'SD ...'           takes the spanwise standard deviation of ...
%	'Scaled SD ...'    scales SD ... by the maximum of ...
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

% parse the variable
if strncmp(var,'Mean',4) || strncmp(var,'SD',2) || strncmp(var,'Scaled SD',9)
    % Remove prefix from var (ie. remove Mean, SD, or Scaled SD)
    if params.ndims == 3
        ny = 1:Ny;
        varorig = var;
        var = strrep(var, 'Mean ', '');
        var = strrep(var, 'Scaled SD ', '');
        var = strrep(var, 'SD ', '');
    else
        error('Mean, SD, and Scaled SD are not supported on 2D data.');
    end
end

% try different densities
if strcmpi(var,'Density')
    rhofiles = dir('rho.*');
    if ~isempty(rhofiles)
        data = spins_reader('rho',ii,nx,ny,nz);
    else
        try
            s = spins_reader('s',ii,nx,ny,nz);
            t = spins_reader('t',ii,nx,nz,nz);
            data = eqn_of_state(t,s);
        catch
            error('Density could not be computed.');
        end
    end
% read in kinetic energy
elseif strcmpi(var,'KE')
    u = spins_reader('u',ii,nx,ny,nz);
    if params.ndims == 3
        v = spins_reader('v',ii,nx,ny,nz);
    else
        v = 0;
    end
    w = spins_reader('w',ii,nx,ny,nz);
    data = 0.5*params.rho_0*(u.^2 + v.^2 + w.^2);
    clearvars u v w
% plot speed (magnitude of velocity vector)
elseif strcmpi(var,'speed')
    u = spins_reader('u',ii,nx,ny,nz);
    if params.ndims == 3
        v = spins_reader('v',ii,nx,ny,nz);
    else
        v = 0;
    end
    w = spins_reader('w',ii,nx,ny,nz);
    data = sqrt(u.^2 + v.^2 + w.^2);
    clearvars u  v w
% read in gradient Richardson number
elseif strcmp(var, 'Ri')
    if length(ny) > 1
        error('Gradient Richardson number must be plotted in x-z plane.')
    else
        % read in data
        rho_z = spins_reader('rho_z',ii,nx,ny,nz)';
        u_z   = spins_reader('u_z',ii,nx,ny,nz)';
        g = params.g;
        rho_0 = params.rho_0;
        % calculate Ri
        Uz_sq = u_z.^2;
        N_sq  = -g/rho_0*rho_z;
        data = (N_sq./Uz_sq)';
        % remove data that is too large
        if max(data(:)) > 5
            data(data>5) = 5;
            warning(['Ri>5 has been set to 5.',...
            'This enables contour and contourf to make meaningful plots.'])
        end
    end
% read in data for plotting streamlines
elseif strcmpi(var,'Streamline')
    if params.ndims == 3
        ny = 1:Ny;
        u = squeeze(mean(spins_reader('u',ii,nx,ny,nz),2));
        w = squeeze(mean(spins_reader('w',ii,nx,ny,nz),2));
    else
        u = spins_reader('u',ii,nx,ny,nz);
        w = spins_reader('w',ii,nx,ny,nz);
    end
    data = ones([size(u) 2]);
    data(:,:,1) = u;
    data(:,:,2) = w;
% read in data for given file name
else
    try
         data = spins_reader(var, ii, nx, ny, nz);
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
