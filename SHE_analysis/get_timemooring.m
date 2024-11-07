close all; clearvars; clc;

% Define input positions for z and x
%z_pos = [.025 .05 .075 .1 .15];
z_pos = [.2 .23 .25 .28 .3];
x_pos = 5;

% Get the 2D grid for x and z, and times
[x, z] = spinsgrid2d;
[times, outputs] = get_output_times;

%% Find the nearest index in x and z to the input positions
xind = nearest_index(x(:, 1), x_pos);

z_ii = NaN.*z_pos; % Initialise z indices array with NaNs
for jj = 1:length(z_pos)
    z_ii(jj) = nearest_index(z(xind, :), z_pos(jj));
end

%% Read the density values for each output time and store in rho_ts
rho_ts = NaN(length(outputs), length(z_pos));

for ii = outputs'
    rho_ts(ii+1, :) = spins_reader_new('rho', ii, xind, z_ii);
    completion(ii, length(outputs));
end

%% Plot end density with "moorings" marked
figure
clrorder = colororder;

subplot(2, 1, 1);
rho = spins_reader_new('rho', ii);
pcolor(x, z, rho); hold on
scatter(x_pos.*ones(length(z_pos), 1), z_pos, 15, clrorder(1:length(z_pos), :), 's', 'filled');
subplot(2, 1, 2);
plot(times, rho_ts);

%% Now plot an animation of the timeseries
figure
for ii = outputs'
    rho = spins_reader_new('rho', ii);
    subplot(2, 1, 1);
    pcolor(x, z, rho); hold on
    scatter(x_pos.*ones(length(z_pos), 1), z_pos, 15, clrorder(1:length(z_pos), :), 's', 'filled');

    hold off
    %ylim([.1 .2])
    subplot(2, 1, 2)
    if ii >= 1
        plot(times(1:ii), rho_ts(1:ii, :));
    end
    xlim([min(times) max(times)])
    pause(1) % Pause to create an animation effect
end