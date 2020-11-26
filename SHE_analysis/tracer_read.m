%TRACER_READ - Plot changing concentration of tracer in wave core for
%different sized waves
%
% Other m-files required: cmocean, spinsgrid2d, spins_params
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 26-Nov-2020; Last revision: 26-Nov-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

clc; close all; clearvars;
%% Set directory names of tracer runs
dir_names = {'../tc1_210720', '../tc2_220720', '../tc3_230720', '../tc4_240720'};

clrmap = cmocean('haline', length(dir_names)); % Set colours of tracer concentration lines

%% Produce empty objects for computational speed
trac_vol = nan(length(dir_names), 100);
wave_pos = trac_vol; % Create the same shaped empty variable for wave_pos
p = gobjects(length(dir_names)); % Empty graphics object for computational speed

%% Run for each directory
for ti = 1:length(dir_names)
    cd(dir_names{ti});
    [x, z] = spinsgrid2d;
    params = spins_params;

    videosave = false;
    fnm = 'Tracer.mp4';
    if videosave
        v = VideoWriter(fnm, 'MPEG-4');
        v.FrameRate = 1;
        open(v);
    end
    if ~isfile('wave_characteristics.mat')
        characterize_wave();
    end
    load wave_characteristics.mat
    volume = spinsvolume();
    wave_window = [wave_center-(wavelength_left*1.5)...
        wave_center+(wavelength_right*1.5)];
    
    %%
    for ii = 0:length(time)-1
        %xind = closest_index(ii+1, 1):closest_index(ii+1, 2);
        xind = find(x(:, 1)<wave_window(ii+1, 2) & x(:, 1)>wave_window(ii+1, 1));
        
        rho = spins_reader_new('rho', ii, xind, []);
        tracer = spins_reader_new('tracer', ii, xind, []);
        
        if ii == 0
            tracer_1 =  spins_reader_new('tracer', 0);
            tracer_tot_vol = sum(sum(tracer_1.*volume));
        end
        trac_vol(ti, ii+1) = sum(sum(tracer.*volume(xind, :)))/tracer_tot_vol;
        
    end
    if videosave
        close(v)
    end
    wave_pos(ti, time+1) = wave_center;
    figure(2)
    p(ti) = plot(wave_pos(ti, :), trac_vol(ti, :), 'Marker', 'x', 'Color', clrmap(ti, :));
    hold on
    name_ti(ti) = join(["amp = ", num2str(WaveStats.meanAmp, 2), " m "]);
    
    %% See exponential fit
    x2ind = find(wave_pos(ti, :)>(params.L_adj+(wavelength_right(1)*1.25)));
    [XOut, YOut] = prepareCurveData(wave_pos(ti, x2ind), trac_vol(ti, x2ind));
    
    f = fit(XOut, YOut, 'exp1');
    fp = plot(f);
    fp.Color = clrmap(ti, :);
    fp.LineStyle = '-';
    coefficient_values(ti, :) = coeffvalues(f)
    amp(ti) = WaveStats.meanAmp;
end
 %%
figure(2)
set(gca, 'XDir', 'normal');
xlabel('Wave Center Position')
ylabel('% of tracer in core');
xlim([wave_pos(end, x2ind(1)) max(wave_pos,[], 'all')])
ylim([0 1]);
legend(p, name_ti)
print('tracer_conc.png', '-dpng');
%%
disp(['Wave 1 extinction coefficient = ' num2str(coefficient_values(1, 2))])
disp(['Wave 2 ext. coefficient = ' num2str(coefficient_values(2, 2))])
disp(['Wave 3 ext. coefficient = ' num2str(coefficient_values(3, 2))])
disp(['Wave 4 ext. coefficient = ' num2str(coefficient_values(4, 2))])

%%
figure(3)
plot(amp, -log(2)./coefficient_values(:, 2), 'k-')
set(gca, 'XDir', 'normal')
xlabel('Wave Amplitude (m)')
ylabel('Half distance (m)')