function [maxRe, initialRe, pos_maxRe] = calc_reynolds(n_layers, isplotting)
%CALC_REYNOLDS - Calculates how the Reynolds changes as wave shoals for 2 layer case.
% Based on Nakayama et al., 2019 definition of Re_{isw}, sets h_1 as
% pycnocline width, h_2 as bottom layer thickness.
%
% Inputs:
%    n_layers - Number of layers in system, currently should be 3 or 2, but
%    could be set up for more
%    isplotting - boolean (true or false) as to if plots should be made
%
% Outputs:
%    maxRe - Maximum Re reached during shoaling
%    initialRe - Re along flat bed of tank
%    pos_maxRE - x position of wave when maximum Re is reached
%
% Other m-files required: spinsgrid2d, spins_params
% Subfunctions: none
% MAT-files required: wavestats.mat
%
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 26-Nov-2020; Last revision: 26-Nov-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

if nargin<2
    isplotting = false;
end
%% Load in wave properties
[x, z] = spinsgrid2d;
params = spins_params;
load('wave_characteristics.mat', 'WaveStats');

rho_ratio = params.delta_rho; % Ratios of densities
g = params.g; 

rho_1 = params.rho_0; % upper layer density
rho_2 = params.rho_0*(1+params.delta_rho); % lower layer density
A = WaveStats.meanAmp;  % Wave amplitude

%% Calculate h parameters
if n_layers == 2
    h_1 = params.h_halfwidth*2; % upper layer thickness
elseif n_layers == 3
    h_1 = -params.pyc_loc + h_halfwidth; % upper layer thickness
else 
    error('This stratification set needs setting up');
end
h_2 = max(abs(z),[],  2) - h_1; % lower layer thickness
H = max(abs(z),[],  2); % Total thickness

%% Calculate linear wave speed
c_0_sqred = rho_ratio*g*(h_1(1)*h_2(1))./H(1); 
c_0 = sqrt(c_0_sqred);

%% Calculate KdV nonlinear term (alpha)
alpha = 3/2 * (c_0/(h_1*h_2)) * (rho_1*(h_2.^2) - rho_2*(h_1^2))./(rho_1*h_2 + rho_2*h_1);

%% Calculate length scale
h_dash = h_1.*h_2./(h_1+h_2);
h_dash(h_dash<0) = NaN;

%% Calculate Re
Re = alpha .* A .* h_dash/params.visco;

%% And statistics of Re
[maxRe, ind] = max(Re); 
initialRe = Re(1); % Identifies the value of Re on flat bed
pos_maxRe = x(ind, 1); %Identifies the position where Re peaks

%% Plot how Re changes
if isplotting
    figure
    plot(x(:, 1)-params.L_adj, Re*1e-3);
    xlabel('x (m)');
    set(gca, 'xlim', [x(1, 1)-params.L_adj x(end, 1)-params.L_adj]);
    ylabel('Re \times 10^{3}');
    set(gca, 'XDir', 'reverse');
end