%PLOT_WAVEFORM - Plots the density, velocities of wave and density profiles in wave and
%upstream
%Primarily for comparison of waveform between DJL/Lab/Numerics. Run in a
%model run directory
%
%
% Outputs:
% Single figure with 4 subplots:
%       Density field top left
%       Horizontal Velocity (u/c) top right
%       Vertical Velocity (w/c) bottom right
%       Density profiles upstream and in wave, bottom left
%
%
% Other m-files required: spinsgrid2d, cmocean
% Subfunctions: none
% MAT-files required: 
% MAT-files requested: wavestats.mat wave_characteristics.mat
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 05-Nov-2020; Last revision: 05-Nov-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
clearvars; close all; clc;
time_i = 38;

%% Load in data
% Data statistics files
if ~isfile('wave_characteristics.mat')
    characterize_wave(0,50); 
end
load wave_characteristics wave_center wavelength_left wavelength_right WaveStats 

c = WaveStats.meanWaveSpeed; % Wave Speed

% Data Field files
[x z] = spinsgrid2d; % Read in spins grid
rho = spins_reader_new('rho', time_i); 
u = spins_reader_new('u', time_i);
w = spins_reader_new('w', time_i);

%% Plot 
xlimits = [wave_center(time_i)-1 wave_center(time_i)+1]+.1; % Set xlim 1m either side of wave center

% Plot Density
subplot(2, 2, 1)
title('density')
pcolor(x, z, rho);
colormap(gca, flip(plasma));
xlim(xlimits);
shading flat
xlabel('x'); ylabel('z')

% Plot Normalised horizontal velocity
subplot(2, 2, 2)
title('u/c')
[~, h] = contourf(x, z, u./c, 25);
colormap(gca, flip(cmocean('curl')));
set(h,'LineColor','none');
xlim(xlimits);
shading flat
xlabel('x'); ylabel('z')
caxis([-1 1])

% Plot normalised vertical velocity
subplot(2, 2, 4)
title('w/c')
[~, h] = contourf(x, z, w./c, 11);
colormap(gca, flip(cmocean('curl')));
set(h,'LineColor','none')
xlim(xlimits);
shading flat
xlabel('x'); ylabel('z')
caxis([-1 1])

% Plot Density profiles
subplot(2, 2, 3)
[~, ii] = min(abs(x(:, 1) - wave_center(time_i).'))
plot(rho(ii , :), z(2048,:), 'r-')
hold on
plot(rho(ii+400, :), z(2048, :), 'b-');
set(gca, 'xdir', 'normal');
xlabel('\rho')
grid on
legend('Mid-wave', 'Upstream', 'Location', 'southwest')