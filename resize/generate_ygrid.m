function yg = generate_ygrid(phys_size,num_size,varargin)
% generate_zgrid  Generates xgrid for use with SPINS data
%
%  Assumptions:
%       - Domain is three dimensional
%       - Cheb grids are only in vertical
%       - Domain is not mapped
%
%  Usage in 3D:
%    yg = generate_xgrid([Lx Ly Lz], [Nx Ny Nz], z_bc_type);

%
%  Inputs:
%    '[Lx Ly Lz]'      - Vector containing physical domain sizes
%    '[Nx Ny Nz]'      - Vector containing numerical grid sizes
%   Optional arguments:
%    'z_bc_type'       - String of z boundary conditions
%                      - these are: 'NO_SLIP', 'FREE_SLIP', 'FOURIER'
%                      - If no optional argument is supplied, function assumes 
%                        grid points are evenly spaced               
%
%  Outputs:
%    'yg'              - ygrid
%
%  Andrew Grace, 2021

if nargin == 3
    if ~strcmp(varargin{1},'NO_SLIP') 
        z_bc_type = '';
        if ~(strcmp(varargin{1},'FOURIER') | strcmp(varargin{1},'FREE_SLIP')) 
            warning('Expected ''NO_SLIP'', ''FOURIER'', or ''FREE_SLIP''. Ignoring argument.');
        end
    else 
        z_bc_type = varargin{1};
    end
elseif nargin > 3
    error('Too many input argument');
elseif nargin < 2
    error('Not enough input arguments');
else
    z_bc_type = '';
end

if ~( length(phys_size) == length(num_size) )
    %Dimensionality of inputs do not agree
    error('Dimensionality of physical and numerical sizes do not agree')
elseif length(phys_size) ~= 3
    %Data must be 3D
    error('Data must be in 3 dimensions')
end

Lx = phys_size(1);
Ly = phys_size(2);
Lz = phys_size(3);

Nx = num_size(1);
Ny = num_size(2);
Nz = num_size(3);


%Generate grid
    dx = Lx/Nx;
    dy = Ly/Ny;
    
    x1d = dx/2:dx:Lx;
    y1d = dy/2:dy:Ly;
    
    if ~strcmp(z_bc_type,'NO_SLIP')
        dz = Lz/Nz;
        z1d = dz/2:dz:Lz;
    else
        [~,z1d] = cheb(Nz-1);
        z1d = Lz/2*(1 - z1d);
    end

    [~,yg,~] = meshgrid(x1d,y1d,z1d); 
    yg = permute(yg, [2 1 3]);

end
