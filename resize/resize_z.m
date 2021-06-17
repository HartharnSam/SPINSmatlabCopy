function field = resize_z(field_name, field, Nz_new, varargin)
% RESIZE_Z  Increase or decrease the resolution in the z-dimension.
%
%  Assumptions:
%    - if 2D, data must be in x-z plane
%    - Boundary condition must be free-slip or periodic
%      - To keep the code simple the periodic BC is not optimized
%      - since the field will still be doubled
%
%  Usage:
%    rho_new = resize_z('rho', rho, 512);
%
%  Inputs:
%    'field_name'   - the name of the field to resize
%    'field'        - the field data
%    'Nz_new'       - number of points to change to
%   Optional arguments:
%    'opts'         - a structure holding simulation parameters
%                   - these are: type_z and Lz, Ly, Lx
%
%  Outputs:
%    'field'        - the resized field
%
%  David Deepwell, 2018

% get number of points in each dimension
sz     = size(field);
Nx     = sz(1);
Nz_old = sz(end);
if length(sz) == 3
    Ny = sz(2);
else
    Ny = 1;
end
% get other parameters
if nargin == 3
    % read from spins.conf
    params = spins_params();
    type_z = params.type_z;
    Lz     = params.Lz;
    Ly     = params.Ly;
    Lx     = params.Lx;    
elseif nargin == 4
    % read from optional input argument
    opts = varargin{1};
    type_z = opts.type_z;
    Lz     = opts.Lz;
    Ly     = opts.Ly;
    Lx     = opts.Lx;     
else
    error('Incorrect number of input arguments. 3 or 4 are expected.')
end


if strcmp(type_z, 'FREE_SLIP') || strcmp(type_z, 'FOURIER')
    dz_old = Lz/Nz_old;
end

% check sizes
mult = Nz_new/Nz_old;
if mult > 1
    increase_points = true;
else
    increase_points = false;
end

if ~strcmp(type_z,'NO_SLIP')
    
    % reshape 2D to match the shape of a 3D field
    if Ny == 1
        field = reshape(field, [Nx 1 Nz_old]);
    end

    % permute to put extending dimension in the first dimension
    field = permute(field, [3 2 1]);

    % extension types
    odd_ext    = @(f) [f; -flipud(f)];
    even_ext   = @(f) [f;  flipud(f)];
    period_ext = @(f) [f;  f];

    % double the field
    if strcmp(type_z, 'FOURIER')
        field = period_ext(field);
    elseif strcmp(field_name, 'w')
        field = odd_ext(field);
    else
        field = even_ext(field);
    end

    % phase shift because grid is cell-centered, and take fft
    dz_new = Lz/Nz_new;
    kzs = fftfreq(2*Nz_old, dz_old);
    phase_corr = exp(-1i*kzs*(dz_old/2-dz_new/2)).';
    field = bsxfun(@times, fft(field), phase_corr);

    % pad or truncate
    if increase_points
        % pad with zeros
        field = [field(1:Nz_old,:,:); zeros(2*Nz_new-2*Nz_old,Ny,Nx); field(end-Nz_old+1:end,:,:)];
    else
        % truncate high frequencies
        field = [field(1:Nz_new,:,:); field(end-Nz_new+1:end,:,:)];
    end

    % take ifft
    field = mult*real(ifft(field));
    field = field(1:Nz_new,:,:);
    % permute back
    field = permute(field, [3 2 1]);
    % return field to 2D if it was
    if Ny == 1
        field = reshape(field,[Nx Nz_new]);
    end
else 
    

    %BC must be NO_SLIP. 
    
    %Make new Cheb grid
    [~,z1d_new] = cheb(Nz_new-1);
    
    %Make old Cheb grid
    [~,z1d_old] = cheb(Nz_old-1);
    
    %Move to physical space
    z1d_old = Lz/2*(1 + z1d_old);
    z1d_new = Lz/2*(1 + z1d_new);
    
    %If data is 3D, we must use interp3
    if Ny ~= 1
        dx = Lx/Nx;
        dy = Ly/Ny;
        field = permute(field, [2 1 3]); %Permute so data agrees with meshgrid
        %Make old grid
        [x_old,y_old,z_old] = meshgrid(linspace(dx/2,Lx - dx/2,Nx),...
            linspace(dy/2,Ly - dy/2,Ny),z1d_old);   
        %Make new grid
        [x_new,y_new,z_new] = meshgrid(linspace(dx/2,Lx - dx/2,Nx),...
            linspace(dy/2,Ly - dy/2,Ny),z1d_new); 
        %Do interpolation
        field = interp3(x_old,y_old,z_old,field,x_new,y_new,z_new);
        %Permute back
        field = permute(field, [2 1 3]);
        %clear grids
        clear x_old y_old z_old x_new y_new z_new;
    else
        dx = Lx/Nx;
        %field = reshape(field, [Nx Nz_old]);
        %Make old grid
        [x_old,z_old] = meshgrid(linspace(dx/2,Lx - dx/2,Nx),...
            z1d_old); 
        %Make new grid
        [x_new,z_new] = meshgrid(linspace(dx/2,Lx - dx/2,Nx),...
            z1d_new); 
        %Interpolate and transpose
        field = (interp2(x_old,z_old,field',x_new,z_new))';
        %Clear grids
        clear x_old z_old x_new z_new;
        
    end

end
%
%End Function

end
