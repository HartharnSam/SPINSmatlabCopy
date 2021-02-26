function plot_diss_uc_Ri(ii, xlimits)
%PLOT_DISS_UC_RI - Plots a diagnostic square of density, dissipation, u/c
%and richardson for a given time 
%
% Inputs:
%    ii - frame number to be plotted (ideally the point of breaking)
%    
%    xlimits - X limits in real coordinates (e.g. [5 7])
%
% Other m-files required: subaxis
% Subfunctions: none
% MAT-files required: none
%
% See also: single_time_pic_example,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 05-Nov-2020; Last revision: 05-Nov-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

%% Read in parameters
params = spins_params;
if nargin<2
    slope_length = params.hill_height/params.hill_slope;
    xlimits = [-slope_length 0]+(params.Lx - params.L_adj);
end
%% Set up x and y grids
[x z] = spinsgrid2d;
x = x-params.L_adj;
xlim_inds = find(x(:, 1)>xlimits(1) & x(:, 1)<xlimits(2));

x = x(xlim_inds, :);
z = z(xlim_inds, :);

%% Read in Data
rho=spins_reader_new('rho',ii, xlim_inds, []);
rho = rho_converter(rho); % convert density to real density
rhoRange = [params.rho_0 params.rho_0+(params.delta_rho*params.rho_0)];
u = spins_reader_new('u',ii, xlim_inds, []);
%w = spins_reader_new('w', ii, xlim_inds, []);
%umaxabs = 6 %max(abs(u(:)));
%diss = spins_reader_new('diss', ii, xlim_inds, []);
%drho_dz = spins_reader_new('rho_z', ii, xlim_inds, []);
%du_dz = spins_reader_new('u_z', ii, xlim_inds, []);

%% Set up figure
figure
clf
colormap darkjet

%% Plot density first 
subaxis(2,2,1)
pcolor(x,z,rho),shading flat; caxis(rhoRange);
%contourf(x,z,rho),shading flat % This is better for printed figures
c = colorbar;
ylabel(c, 'Diss')
ylabel('z (m)')
xticklabels([]);
xlim(xlimits)

%% Plot Dissipation
subaxis(2,2,2)
pcolor(x,z, log10(diss)); shading flat; caxis([-10 0]);
c = colorbar;
ylabel(c, 'log_10 (disssipation)')
ylabel('z (m)')
xlabel('x (m)')
xlim(xlimits)

%% Plot Horizontal Velocity Next
subaxis(2,2,3)
try load('wavestats.mat', 'WaveStats')
    wave_speed = WaveStats.meanWaveSpeed;
catch
    warning('No WaveStats saved, run characterize_wave if normalised velocity fields requested')
end

pcolor(x,z,u./wave_speed),shading flat,caxis([-1 1])
c = colorbar;
hold on 
contour(x, z, u./wave_speed, [-1 1], 'k-');
ylabel(c, 'U/c')
ylabel('z (m)')
xlabel('x (m)')
xlim(xlimits)

%% Plot Richardson
subaxis(2,2,4)
g = params.g;
Ri = (g./params.rho_0) *drho_dz./(du_dz.^2);

pcolor(x,z,Ri); shading flat; caxis([0 .5]);
c = colorbar;
hold on 
%contour(x, z, Ri, [0 .25], 'k-');

ylabel(c, 'Ri')
ylabel('z (m)')
xlabel('x (m)')
xlim(xlimits)

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------