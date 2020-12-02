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
isVideo = false; % Switch to make movie
n_layers = 3;
if n_layers == 2;
    isTwoLayer = true; % Switch for if two layer (affects isopycnal chosen for characterize wave)
else
    isTwoLayer = false;
end
RunDirectoryName = {'28_121120', 
    }; % List of directories to process files from

%% Set location for video file to be saved to
pathname = ['C:\Users\', getenv('username'), '\OneDrive - Newcastle University\Shared_Videos\Numerics\'];

%% Run loops
outputs = NaN(length(RunDirectoryName), 16);
for ii = 1:length(RunDirectoryName)
    cd(['../', RunDirectoryName{ii}]);
    
    params = spins_params;
    disp(params.name);
    disp(RunDirectoryName{ii});
    
    %% Calculate/Load wave characteristics
    if n_layers == 1 % exist('wave_characteristics.mat', 'file') == 0 % Has characterize_wave been run before for this experiment?
        if n_layers == 1
            WaveStats = wave_characteriser([10 30]); % Less robust script for use with continuous stratification
        else
            WaveStats = characterize_wave(isTwoLayer, [0 50]); % runs the script to find wavelengths, amplitudes, wave position, speed etc.
        end
        plot_wave_char; % Plots the outputs and saves them
        close all
    else % If this has already been run before, just need to load it
        load('wave_characteristics.mat', 'WaveStats');
    end
    
    %% Clean/Plot diagnostics
    if exist('diagnostics.mat', 'file') == 0 % checks for file created by clean_diagnostics
        clean_diagnostics(); % Cleans diagnostics files
    end
    
    if exist('all_diagnos.mat', 'file') == 0% checks for file created by plot_diagnos
        diagnos = plot_diagnos(isPlotDiagnostics, false, isPlotDiagnostics);
        close all
        plot_stress(isPlotDiagnostics, isPlotDiagnostics); % Also plot stresses
        close all
    else
        diagnos = load('all_diagnos.mat');
        diagnos = diagnos.all_diagnos;
    end
    
    %% plot movie
    if isVideo
        if ispc
            SPINS_movie_maker({'rho', 'vorty', 'u_normalised', 'w_normalised'}, 'slopeonly', WaveStats.endSlope, true, [pathname, params.name, '.mp4']);
        else
            SPINS_movie_maker({'rho', 'vorty', 'u_normalised', 'w_normalised'}, 'slopeonly', WaveStats.endSlope, true, [pathname, params.name, '.avi']);
        end
    end
    
    %% Calculate Reynolds statistics
    [maxRe, startRe] = calc_reynolds(2);
    
    %% Collate stats to be pasted into appropriate index sheet
    
    outputs(ii, :) = [params.hill_height, params.hill_slope, params.pyc_adj_loc, params.Nx,...
        params.Ny, params.Nz, 16, diagnos.diagnos.TotClockTime/3600,...
        diagnos.diagnos.AvgClockTimePerSimSec, WaveStats.meanAmp, ...
        WaveStats.meanWaveSpeed, WaveStats.meanWavelength...
        WaveStats.meanAmp/WaveStats.meanWavelength,...
        params.hill_slope/sqrt(WaveStats.meanAmp/WaveStats.meanWavelength), maxRe, startRe];
    %% Clean up
    clearvars -except is* RunDirectoryName pathname outputs* n_layers
    close all
end
outputs_names = {'HillHeight', 'HillSlope', 'PycAdjLoc', 'Nx', 'Ny', 'Nz', 'NumProcessors', 'TotalComputationTime',...
    'ComputationSpeed (clock time per sim second)', 'Amplitude', 'Wave Speed', 'Wavelength', 'Wave Steepness', 'Ir', 'Max Re', 'Initial Re'};

Outputs = array2table(outputs, 'VariableNames', outputs_names);
clear outputs*