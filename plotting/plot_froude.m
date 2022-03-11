function plot_stress(calc_type, save_plots)
if nargin == 0
    make_plots = true;
end

params = spins_params;
if exist('wave_characteristics.mat', 'file') == 0 %
    characterize_wave;
end
load('wave_characteristics.mat', 'wave_speed', 'time', 'WaveStats')

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
    end
end

%% Now calculate Fr
Fr = u_max ./ wave_speed';
plot(time, Fr, 'k-');
yline(1);

if save_plots
    exportgraphics(gcf, 'Froude_timeseries.png');
end