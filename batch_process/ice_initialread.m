%Ice_InitialRead
% This m file should be placed in the parent folder for simulations
clc; clearvars; close all;

% Parameters to set
list_of_sims = {'pyc1_1'};  % Names of folders to process
isPlotDiagnostics = true; % Do we want to (re)make plots of diagnostics (yes)
isCalcCharacteristics = true;
isVideo = true; % Do we want to (re)make the main video

pathname = '.\';

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%cd(list_of_sims{1});

%% Run loops
for ii = length(list_of_sims)
    cd(['../', list_of_sims{ii}]); % Move into simulation folder
    
    params = spins_params; % Load in simulation parameters
    disp(['Simulation ', list_of_sims{ii}]);
    
    %% Calculate/Load wave characteristics
    if exist('wave_characteristics.mat', 'file') == 2 && ~isCalcCharacteristics
        load('wave_characteristics.mat', 'WaveStats');
    else
        warning('off', 'MATLAB:DELETE:FileNotFound');
        delete('wave_characteristics.mat'); 
        %if exist('wave_characteristics.mat', 'file') == 0 % Has characterize_wave been run before for this experiment?
        WaveStats = characterize_wave(false); % runs the script to find wavelengths, amplitudes, wave position, speed etc.
        plot_wave_char; % Plots the outputs and saves them
        close all
    end
    
    %% Clean/Plot diagnostics
    if exist('diagnostics.mat', 'file') == 0 % checks for file created by clean_diagnostics
        clean_diagnostics(); % Cleans diagnostics files
    end
    
    if exist('all_diagnos.mat', 'file') == 0% checks for file created by plot_diagnos
        diagnos = plot_diagnos(isPlotDiagnostics, false, isPlotDiagnostics);
        close all
        plot_stress(isPlotDiagnostics, isPlotDiagnostics); % Also plot stresses
    else
        close all
        diagnos = load('all_diagnos.mat');
        diagnos = diagnos.all_diagnos;
    end
    %plot_froude('full', true);
    %% plot movie
    if isVideo
        if ispc
            spins_movie_maker({'rho', 'vorty'},...
                'slopeonly', WaveStats.endTank, true, fullfile(pathname, ['movie', '.mp4']));            
        else
            SPINS_movie_maker({'rho', 'vorty', 'u_normalised', 'w_normalised', 'tracer', 'diss'},...
                'slopeonly', WaveStats.endTank, true, fullfile(pathname, [params.name, '.avi']));
        end
    end
    close all
end
