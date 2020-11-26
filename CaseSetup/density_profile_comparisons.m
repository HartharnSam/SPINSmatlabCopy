%DENSITY_PROFILE_COMPARISONS - Compares initial density in the lab and numerics
%
% Other m-files required: spins_params, xgrid_reader, zgrid_reader,
% spins_reader_new, rho_converter, betterplots, cmocean
%
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
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
clr = ['b','k','r','g','y','c','m']; % Order of colours to plot lines in (don't think it is actually used though)

%% Load in SPINS grid params
params = spins_params;
x = xgrid_reader();
z = zgrid_reader();

% and the density field
rho = spins_reader_new('rho', 0);
rho = rho_converter(rho);

%% Identifies two points to plot (one behind gate, another midway along flat bed)
ind_gate = find(x(:, 1)<params.L_adj);
line_1 = ind_gate(round(length(ind_gate)/2));

ind_free_wave = find(x(:, 1)>params.L_adj*1.15 & x(:, 1)<(params.Lx - ...
    (params.hill_height/params.hill_slope)*1.15));
line_2 = ind_free_wave(round(length(ind_free_wave)/2));


betterplots; % Selects some formatting options
colormap(cmocean('dense'));

%% Plot SPINS Density
subplot(2, 2, [1 2])

pcolor(x,z,rho),shading flat
%contourf(x,z,rho),shading flat % This is better for printed figures

title('Density')
ylabel('z (m)')
colorbar
set(gca, 'XDir', 'reverse')
hold on
plot([x(line_1, 1) x(line_1, 1)], [-0.3 0], '-b'); % plots lines indicating where profiles are from
plot([x(line_2, 1) x(line_2, 2)], [-0.3 0], '-', 'Color', [ 0.9100 0.4100 0.1700]);

%% Plot SPINS Density Profile behind gate
subplot(2, 2, 4)
plot(rho(line_1, :), z(line_1, :), '-b')
xlabel('\rho')
ylabel('z (m)')
set(gca, 'XDir', 'normal');
grid off

%% Plot SPINS Density Profile on flat bed
subplot(2, 2, 3)
plot(rho(line_2, :), z(line_2, :), '-', 'Color', [ 0.9100 0.4100 0.1700]);
xlabel('\rho')
ylabel('z (m)')
set(gca, 'XDir', 'normal');
grid off

%% Also plot example probe output (090120.mat)
load('090120.mat', 'data'); % Load in this data
hold on
plot(data.FittedData(:, 2), data.FittedData(:, 1)-.3,'.','Color',clr(1));
set(gca, 'XDir', 'normal');

%print(['C:\Users\samha\OneDrive - Newcastle University\Project\Shoal_Core\figures\density_profile,'.png'], '-dpng');

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------