function field = resize_xyz(field_name, field, new_grid, varargin)
% RESIZE_XYZ  Increase or decrease the resolution in all three dimensions
%
%  Assumptions:
%    - For 3D data only
%    - Boundary condition must be free-slip or periodic
%      - To keep the code simple the periodic BC is not optimized
%      - since the field will still be doubled
%
%  Usage:
%    rho_new = resize_x('field_name', field, [Nx_new Ny_new Nz_new]);
%
%  Inputs:
%    'field_name'   - the name of the field to resize
%    'field'        - the field data
%    'new_grid'     - vector of number of points to change
%   Optional arguments:
%    'opts'         - a structure holding simulation parameters
%                   - these are: type_x and Lx
%
%  Outputs:
%    'field'        - the resized field
%
%  Authors: David Deepwell, 2018, Adapted by Andrew Grace 2021

field = resize_x(field_name,field,new_grid(1));
field = resize_y(field_name,field,new_grid(2));
field = resize_z(field_name,field,new_grid(3));
