function [Fr, times] = plot_froude(calc_type, save_plots)
%PLOT_FROUDE - Calculates the Froude number for each simulation timestep
% Fr = u_max / c 
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 22-Mar-2022; Last revision: 22-Mar-2022
% MATLAB Version: 9.12.0.1884302 (R2022a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

if nargin == 0
    save_plots = true;
end

params = spins_params;
if exist('wave_characteristics.mat', 'file') == 0 %
    characterize_wave;
end
load('wave_characteristics.mat', 'wave_speed', 'time', 'WaveStats')
[x, z] = spinsgrid2d;


if strcmpi(calc_type, 'fast')
    if exist('all_diagnos.mat', 'file') == 0 %
        plot_diagnos;
    end
    load('all_diagnos.mat');
    % Bring down resolution
    u_max = interp1(all_diagnos.diagnos.Time, all_diagnos.diagnos.Max_u, time);
    
else
    final_length = length(dir('u.*'));
    time = interp1(0:length(time)-1, time, 0:final_length-1, 'linear','extrap');
    wave_speed(end:final_length) = WaveStats.meanWaveSpeed;
    for i = 1:length(time)
        u = spins_reader_new('u',time(i));
        u_max(i) = max(u,[], 'all');
        x_ind = nearest_index(x(:, 1), wave_center);
        [~, wave_speed(i)] = calc_kdv_depthchange(max(abs(z(x_ind, :))));
    end
end

%% Now calculate Fr
Fr = u_max ./ wave_speed';
plot(time, Fr, 'k-');
yline(1);

if save_plots
    exportgraphics(gcf, 'Froude_timeseries.png');
end