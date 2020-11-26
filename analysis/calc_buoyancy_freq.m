%CALC_BOUYANCY_FREQ - Calculates + plots mid-tank buoyancy frequency profile in SPINS output
%
% Other m-files required: spinsgrid2d, spins_params, spins_reader_new
% Subfunctions: none
% MAT-files required: none
%
% See also: 
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 26-Nov-2020; Last revision: 26-Nov-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
clc; clearvars; close all;
%% Load in SPINS grid and parameters
spinsgrid2d;
params = spins_params; 
%% Load data
ii = 0; % Time = 0
Nx = params.Nx;
drho_dz = spins_reader_new('rho_z',ii); % Load SPINS Deriv file
N_2 = -(params.g/params.rho_0)*drho_dz;

figure
plot(N_2, z(Nx/2, :));
xlabel('N^2 Buoyancy Frequency (s^{-2})'); ylabel('z (m)');

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------