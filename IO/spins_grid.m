function gd = spins_grid(varargin)
%  SPINS_GRID  load the grid from SPINS.
%
%  Usage:
%    gd = spins_grid() gives the vector grid in a structure
%
%  Inputs:
%    Optional arguments:
%	'Vector'    - (default) gives vector grid in output structure
%	'Full'      - gives the entire grid in output structure
%	'FastFull'  - gives the entire grid, as calculated from parameters in spins.conf
%
%  Outputs:
%    gd		- structure containing the grid
%
%  David Deepwell, 2015

% check input arguments
% if nargin > 1
%     error('Too many input arguments. spins_grid accepts "Vector", "FastFull", or "Full".');
% elseif nargin == 1 && ~strcmpi(varargin,'Vector') && ...
%                       ~strcmpi(varargin,'Full') && ~strcmpi(varargin,'FastFull')
%     error('Argument not understood. spins_grid accepts "Vector", "FastFull", or "Full".');
% end
% 
% % check if parameters are set in spins.conf
% if nargin == 0 || strcmpi(varargin,'Vector') || strcmpi(varargin,'FastFull')
%     params = spins_params();
%     eNx = isfield(params,'Nx');
%     eNy = isfield(params,'Ny');
%     eNz = isfield(params,'Nz');
%     eLx = isfield(params,'Lx');
%     eLy = isfield(params,'Ly');
%     eLz = isfield(params,'Lz');
%     emin_x = isfield(params,'min_x');
%     emin_y = isfield(params,'min_y');
%     emin_z = isfield(params,'min_z');
%     etype_x = isfield(params,'type_x');
%     etype_y = isfield(params,'type_y');
%     etype_z = isfield(params,'type_z');
% 
%     % are all fields present?
%     e3D = eNx*eNy*eNz*...
%           eLx*eLy*eLz*...
%           emin_x*emin_y*emin_z*...
%           etype_x*etype_y*etype_z;
%     e2D = eNx*eNz*...
%           eLx*eLz*...
%           emin_x*emin_z*...
%           etype_x*etype_z;
%     % find dimensionality
%     if eNx*eNy*eNz == 1
%         Nx = params.Nx;
%         Ny = params.Ny;
%         Nz = params.Nz;
%         if Nx == 1 || Ny == 1 || Nz == 1
%             e3D = 0;
%         end
%     end
% end
% 
% if nargin == 0 || strcmpi(varargin,'Vector')
%     % create grid if everything is in spins.conf
%     if e3D == 1 || e2D == 1
%         gd = make_vector(params,e3D);
%     else % read grid
%         if exist('xgrid', 'file')
%             gd.x = xgrid_reader();
%             %('xgrid',[],1,1);
%         else
%             warning('xgrid not found.')
%         end
%         if exist('ygrid', 'file')
%             gd.y = ygrid_reader();
%             %gd.y = spins_reader_new('ygrid',1,[],1);
%         end
%         if exist('zgrid', 'file')
%             gd.z = zgrid_reader();
%             %gd.z = spins_reader_new('zgrid',1,1,[]);
%         else
%             warning('zgrid not found.')
%         end
%     end
% % read full grid 'fast'
% elseif strcmpi(varargin,'FastFull') && nargin == 1
%     if ~e3D && ~e2D
%         error('FastFull option not possible. Not all parameters are set in spins.conf');
%     else
%         gd_vec = make_vector(params,e3D);
%         if e3D
%             gd.x = repmat(gd_vec.x',1,Ny,Nz);
%             gd.y = repmat(gd_vec.y,Nx,1,Nz);
%             z_s = spins_reader_new('zgrid',[],1,[]);
%             z_s = reshape(z_s,[Nx, 1, Nz]);
%             gd.z = repmat(z_s,1,Ny,1);
%         elseif e2D
%             gd.x = repmat(gd_vec.x',1,Nz);
%             gd.z = zgrid_reader();
%         end
%     end
% % read full grid
% elseif strcmpi(varargin,'Full') && nargin == 1
%     if exist('xgrid', 'file')
%         gd.x = xgrid_reader;
%     else
%         warning('xgrid not found.')
%     end
%     if exist('ygrid', 'file')
%         gd.y = ygrid_reader;
%     end
%     if exist('zgrid', 'file')
%         gd.z = zgrid_reader;
%     else
%         warning('zgrid not found.')
%     end
% end
% 
% % error message if no grid files are in directory
% if ~exist('gd','var')
%     error('The grid was not found in the working directory. Do you know where you are?')
% end
% end
% 
% function gd_vec = make_vector(params,e3D);
%     Lx = params.Lx;
%     Lz = params.Lz;
%     Nx = params.Nx;
%     Nz = params.Nz;
%     min_x = params.min_x;
%     min_z = params.min_z;
%     type_x = params.type_x;
%     type_z = params.type_z;
%     if ~strcmp(type_x,'NO_SLIP') && ~strcmp(type_x,'CHEBY')
%         gd_vec.x = min_x + Lx*([0:(Nx-1)]+0.5)/Nx;
%     else
%         gd_vec.x = min_x + Lx*0.5*(1 - cos(pi*[0:(Nx-1)]/(Nx-1)));
%     end
%     if ~strcmp(type_z,'NO_SLIP') && ~strcmp(type_z,'CHEBY')
%         gd_vec.z = min_z + Lz*([0:(Nz-1)]+0.5)/Nz;
%     else
%         gd_vec.z = min_z + Lz*0.5*(1 - cos(pi*[0:(Nz-1)]/(Nz-1)));
%     end
%     if e3D == 1
%         Ly = params.Ly;
%         Ny = params.Ny;
%         min_y = params.min_y;
%         type_y = params.type_y;
%         if ~strcmp(type_y,'NO_SLIP') && ~strcmp(type_y,'CHEBY')
%             gd_vec.y = min_y + Ly*([0:(Ny-1)]+0.5)/Ny;
%         else
%             gd_vec.y = min_y + Ly*0.5*(1 - cos(pi*[0:(Ny-1)]/(Ny-1)));
%         end
%     end
% end


[gd.x gd.z] = spinsgrid2d;
if isfile('ygrid')
    gd.y = ygrid_reader;
end

%gd.x = x; 
%gd.z = z;
% gd.z1d = squeeze(z(1,:));
% gd.x1d = squeeze(x(:,1));
% sz = size(x);
% gd.NX = sz(1);
% gd.NZ = sz(2);
% gd.sz = sz;
% gd.dx = x(2,2)-x(1,1);
% gd.Lz = max(abs(z(:)));
% gd.Lx = sz(1)*dx;
% [gd.wtz, wt]=clencurt(sz(2)-1); gd.wt = wt*Lz/2;
