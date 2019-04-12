function data = spins_readdata(varname, ii, nx, ny, nz)
%  SPINS_READDATA  Read in data for a given variable
%
%  Usage:
%    data = spins_readdata('var', t_i, nx, ny, nz) reads in var at positions (nx,ny,nz) at t_i
%                          where nx, ny, nz are scalars or 1D vectors
%
%  Inputs:
%    'varname' may be of different forms:
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
if strncmp(varname,'Mean',4) || strncmp(varname,'SD',2) || strncmp(varname,'Scaled SD',9)
    % Remove prefix from varname (ie. remove Mean, SD, or Scaled SD)
    if params.ndims == 3
        ny = 1:Ny;
        varorig = varname;
        varname = strrep(varname, 'Mean ', '');
        varname = strrep(varname, 'Scaled SD ', '');
        varname = strrep(varname, 'SD ', '');
    else
        error('Mean, SD, and Scaled SD are not supported on 2D data.');
    end
end

% try different densities
if strcmpi(varname,'Density')
    rhofiles = dir('rho.*');
    if ~isempty(rhofiles)
        data = spins_reader('rho',ii,nx,ny,nz);
    else
        try
            s = spins_reader('s',ii,nx,ny,nz);
            t = spins_reader('t',ii,nx,ny,nz);
            data = eqn_of_state(t,s);
        catch
            error('Density could not be computed.');
        end
    end
% read in kinetic energy
elseif strcmpi(varname,'KE')
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
elseif strcmpi(varname,'speed')
    u = spins_reader('u',ii,nx,ny,nz);
    if params.ndims == 3
        v = spins_reader('v',ii,nx,ny,nz);
    else
        v = 0;
    end
    w = spins_reader('w',ii,nx,ny,nz);
    data = sqrt(u.^2 + v.^2 + w.^2);
    clearvars u  v w
% rotate horizontal velocity vectors (for tilted tank cases)
elseif strcmpi(varname,'up')
    u = spins_reader('u',ii,nx,ny,nz);
    w = spins_reader('w',ii,nx,ny,nz);
    theta = params.tilt_angle;
    data = u*cosd(theta) + w*sind(theta);
    clearvars u w
% rotate horizontal velocity vectors (for tilted tank cases)
elseif strcmpi(varname,'wp')
    u = spins_reader('u',ii,nx,ny,nz);
    w = spins_reader('w',ii,nx,ny,nz);
    theta = params.tilt_angle;
    data = w*cosd(theta) - u*sind(theta);
    clearvars u w
elseif strncmp(reverse(varname), 'ms', 2) || ... % ends in sm - Spanwise Mean
        strncmp(reverse(varname), 'zx', 2)       % ends in xz - a cross-section
    suf = varname(end-2:end);
    if strncmpi(varname,'up', 2)
        u = xz_reader(['u',suf], ii, nx, nz);
        w = xz_reader(['w',suf], ii, nx, nz);
        theta = params.tilt_angle;
        data = u*cosd(theta) + w*sind(theta);
        clearvars u w
    elseif strncmpi(varname,'wp', 2)
        u = xz_reader(['u',suf], ii, nx, nz);
        w = xz_reader(['w',suf], ii, nx, nz);
        theta = params.tilt_angle;
        data = w*cosd(theta) - u*sind(theta);
        clearvars u w
    else
        data = xz_reader(varname, ii, nx, nz);
    end
elseif strncmp(reverse(varname), 'mottob', 6) || ... % bottom surface
        strncmp(reverse(varname), 'pot', 3) % top surface
    data = xy_reader(varname, ii, nx, ny);
% read in gradient Richardson number
elseif strcmp(varname, 'Ri')
    if length(ny) > 1
        error('Gradient Richardson number must be plotted in x-z plane.')
    else
        % read in data
        rho_z = spins_reader('rho_z',ii,nx,ny,nz)';
        u_z   = spins_reader('u_z',ii,nx,ny,nz)';
        g = params.g;
        if params.delta_rho < 1
            rho_0 = 1;
        else
            rho_0 = params.rho_0;
        end
        % calculate Ri
        N_sq  = -g/rho_0*rho_z;
        data = (N_sq./u_z.^2)';
        % remove data that is too large
        if max(data(:)) > 5
            data(data>5) = 5;
            warning(['Ri>5 has been set to 5.',...
            'This enables contour and contourf to make meaningful plots.'])
        end
    end
% read in data for plotting streamlines
elseif strcmpi(varname,'Streamline')
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
    data = spins_reader(varname, ii, nx, ny, nz);
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
