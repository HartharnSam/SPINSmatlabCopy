%INITIAL_CONDITION_MAKER - sets up grids etc as in SPINS, based on input like spins.conf
% outputs a figure equivalent to running the first frame of SPINS
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% June-2020; Last revision: 26-Nov-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
clc; clearvars; close all;
%% Parameters to set to match spins.conf file
% Grid parameters
Lx = 7.0;
Ly = 0.1;
Lz = 0.3;
Nx = 4096;
Ny = 1;
Nz = 256;
min_x=  0;
min_y = 0;
min_z = -.3;

% Physical Parameters
g = 9.81;
rot_f = 0.0e-3;
rho_0 = 1026.0;
visco = 1e-6;
kappa_rho = 1e-7;

% Problem Parameters
delta_rho = (1046-rho_0)/rho_0;
pyc_loc = -0.15;
h_halfwidth = 0.075;
pyc_adj_loc = -0.2237;
h_adj_halfwidth = h_halfwidth;
delta_x = 0.04;
L_adj = 0.4;

% Topography Parameters
hill_height = 0.3;
hill_slope = 0.4;
hill_trans = 0.01;
hill_end_dist = -0.0;

%% Set up grid (don't change this bit)

x = min_x:Lx/(Nx-1):Lx;
y = min_y:Ly/(Ny-1):Ly;
z = 0:Lz/(Nz-1):(Lz);
[xx, yy, ~] = meshgrid(x, y, z);

xx = squeeze(xx);

ii = repmat(0:Nz-1, Nx,1);
kk = repmat(0:Nx-1, 1, Nz);

% set up Chebyshev in vertical
zz = -cos(ii.*pi./(Nz-1));

%% Set up topography and zg (z grid) based on problem topography parameters and equations in .cpp code
hill_length = hill_height/hill_slope;

a1 = Lx-hill_length-hill_end_dist;
d1 = hill_trans;
a2 = Lx-hill_end_dist;
d2 = hill_trans;

topo = hill_slope*((.5*(xx-a1+(d1*log(2*cosh((xx-a1)/d1))))) - ...
    (0.5*(xx-a2 +(d2*log(2*cosh((xx-a2)/d2))))));

zg = min_z + 0.5*Lz*(1+zz) + 0.5*(1-zz).*topo;

%% Plot these things upfigure
% Note Units/absolute numbers aren't relevant, just how they change within
% a plot is useful
figure
subplot(4, 1, 1)
plot(x, topo(:, 1), 'k-');
ylabel('z (m)');
title('Height of topography above base')
colorbar
set(gca, 'xlim', [min_x Lx]);

subplot(4, 1, 2)
pcolor(squeeze(xx), squeeze(zz), squeeze(min_z + 0.5*Lx*(1+zz))); shading flat; colorbar
ylabel('z (m)');
title('zgrid with no topography')

subplot(4, 1, 3)
pcolor(squeeze(xx), squeeze(zz), squeeze(0.5*(1-zz).*topo)); shading flat; colorbar
title('impact of topography')
ylabel('z (m)');

subplot(4, 1, 4); 
pcolor(squeeze(xx), squeeze(zz), squeeze(zg)); shading flat; colorbar
title('Final ZGrid');
xlabel('x (m)');
ylabel('z (m)');

%% Identify behind gate and flat bed profile locations
% Behind gate profile
inds_gate = find(xx(:, 1)<L_adj);
line_1 = inds_gate(round(length(inds_gate)/2));

% Flat bed profile
inds_free_wave = find(xx(:, 1)>L_adj*1.15 & xx(:, 1)<(Lx - ...
    (hill_height/hill_slope)*1.15));
line_2 = inds_free_wave(round(length(inds_free_wave)/2));

%% Calculate density profile (based on problem parameters and .cpp equations)

rho =  -0.5 * delta_rho * tanh((zg-pyc_loc)/h_halfwidth); % Without the gate

rho = rho.*0.5.*(1.0+tanh((xx-L_adj)/delta_x)); % Clears the region behind the gate
rho = rho + 0.5*(1.0-tanh((xx-L_adj)/delta_x))...
    .*(-0.5*delta_rho*tanh((zg-pyc_adj_loc)/h_adj_halfwidth)); % Adds in density around the gate

% Conversion to real units
rho = rho_0.*(1 + delta_rho/2 + rho.*(delta_rho.*(1+delta_rho) + 1));

%% Plot densities
figure(2)
subplot(2, 2, [1 2])

pcolor(xx,zg,rho),shading flat
title('Density')
%contourf(x,z,rho),shading flat % This is better for printed figures
ylabel('z (m)')
colorbar
set(gca, 'XDir', 'reverse')
hold on
plot([xx(line_1, 1) xx(line_1, 1)], [-0.3 0], '-b');
plot([xx(line_2, 1) xx(line_2, 2)], [-0.3 0], '-', 'Color', [ 0.9100 0.4100 0.1700]);

%% Plot Gate Water Column
subplot(2, 2, 4)
plot(rho(line_1, :), zg(line_1, :), '-b')
xlabel('\rho')
ylabel('z (m)')
set(gca, 'XDir', 'normal')
%% Plot Main Water Column
subplot(2, 2, 3)

plot(rho(line_2, :), zg(line_2, :), '-', 'Color', [ 0.9100 0.4100 0.1700]);
xlabel('\rho')
ylabel('z (m)')
set(gca, 'XDir', 'normal')

%% Tracer
% figure
% subplot(2, 2, [1 2])
% 
% pcolor(xx,zg,rho),shading flat
% title('Density')
% %contourf(x,z,rho),shading flat % This is better for printed figures
% ylabel('z (m)')
% colorbar
% set(gca, 'XDir', 'reverse')
% 
% 
% % Set up dye strong in top layer
% dye = 0.5 + 0.5*tanh((zg-pyc_adj_loc)/(h_adj_halfwidth));
% dye = 0.5*dye.*(1+tanh((L_adj-xx)/delta_x));
% 
% 
%  subplot(2, 2, [3 4])
% pcolor(xx, zg, dye), shading flat;
% colorbar
% title('Dye')


