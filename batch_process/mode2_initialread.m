%%INITIAL_READ - for initial processing of a SPINS output
% also useful for a function directory
%
% Other m-files required: loads
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
clc; clearvars; close all;

%% Set user options
isPlotDiagnostics = true; %Switch to plot diagnostics
isVideo = true; % Switch to make movie

RunDirectoryName = {'IC_asym4', 'IC_step2', 'IC_step3'}; % List of directories to process files from

%% Set location for video file to be saved to

%% Run loops
for ii = 2:length(RunDirectoryName)

    cd(['../', RunDirectoryName{ii}]);

    params = spins_params;
    disp(RunDirectoryName{ii});

    %% Calculate/Load wave characteristics
    try
        lastwarn('');
        load('wave_characteristics.mat', 'WaveStats');
        msg = lastwarn;
        if ~isempty(msg)
            error('h')
        end
    catch
        warning('off', 'MATLAB:DELETE:FileNotFound');
        delete('wave_characteristics.mat'); delete('wavestats.mat')
        WaveStats = characterize_wave(false, [5 50]); % runs the script to find wavelengths, amplitudes, wave position, speed etc.
        plot_wave_char; % Plots the outputs and saves them
        close all
    end
    [~, direct] = fileparts(cd);

    %% plot movie
    if isVideo
        if ispc
            SPINS_movie_maker({'rho', 'vorty', 'u_normalised', 'w_normalised'},...
                [0 6.6], 5, true, [direct '_full_output.mp4'])
        else
            SPINS_movie_maker({'rho', 'vorty', 'u_normalised', 'w_normalised'}, 'slopeonly', WaveStats.endSlope, true, fullfile(pathname, [params.name, '.avi']));
        end
    end
    close all

    %% Plot isopycnal movie
    [x, z] = spinsgrid2d;
    figure(1);
    clf
    set(gcf, 'Units', 'centimeters'); set(gcf, 'Position', [1 1 44 20]);
    set(gcf, 'PaperUnits','centimeters'); set(gcf, 'Position', [1 1 44 20]);
    vid = VideoWriter([direct '_isopycnals.mp4'], 'MPEG-4');
    vid.Quality = 100;
    vid.FrameRate = 1;
    open(vid);
    
    for tt = 0:100
        rho = spins_reader_new('rho', tt);
        contour(x, z, rho, [-params.delta_rho_1 params.delta_rho_2]/2, '-k');
        hold on
        plot(x(:, 1), z(:, 1), 'k-')
        hold off
        drawnow;
        figure_print_format(gcf);
        F = getframe(gcf);
        writeVideo(vid, F);
    end
    close(vid);
end