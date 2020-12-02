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
outputs = NaN(length(RunDirectoryName), 21);
for ii = 1:length(RunDirectoryName)
    cd(['../', RunDirectoryName{ii}]);
    
    params = spins_params;
    disp(params.name);
    disp(RunDirectoryName{ii});
    
    %% Calculate/Load wave characteristics
    if exist('wave_characteristics.mat', 'file') == 0 % Has characterize_wave been run before for this experiment?
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
    slope_condition_tmp = num2str(round(params.hill_slope*1.5 *1000));
    if params.Lx>7
        slope_condition_tmp = [slope_condition_tmp, 'mm_b'];
    else
        slope_condition_tmp = [slope_condition_tmp, 'mm'];
    end 
    slope_condition{ii} = slope_condition_tmp;
    
    outputs(ii, :) = [params.hill_height, params.hill_slope, params.pyc_adj_loc, ...
        params.h_halfwidth, (params.Lz + (params.pyc_loc - params.h_halfwidth)),...
        WaveStats.meanAmp, WaveStats.meanWaveSpeed, WaveStats.meanWavelength... % Wave amplitude, speed, wavelength
        WaveStats.meanAmp/WaveStats.meanWavelength,... % wAve steepness
        params.hill_slope/sqrt(WaveStats.meanAmp/WaveStats.meanWavelength),... % Ir
        NaN, NaN, params.Nx, params.Ny, params.Nz, params.Lx, params.Ly, params.Lz,... % Resolution
        params.Lx./params.Nx, params.Ly./params.Ny, params.Lz./params.Nz ...% More resolution
        ];
    %% Clean up
    clearvars -except is* RunDirectoryName pathname outputs* n_layers slope_condition
    close all
end
outputs_names = {'HillHeight', 'HillSlope', 'PycAdjLoc', 'pyc_thickness', ...
    'h1', 'amplitude', 'speed', 'wavelength', 'Wave Steepness', 'Ir', 'Classification', ' ',...
    'Nx', 'Ny', 'Nz', 'Lx', 'Ly', 'Lz', 'Lx/Nx', 'Ly/Ny', 'Lz/Nz'};

Outputs = array2table(outputs, 'VariableNames', outputs_names);
% Outputs.FileNames = RunDirectoryName;
% Outputs.SlopeCondition = slope_condition;
% Outputs.InitialWaveCond = zeros(length(RunDirectoryName), 0);
% Outputs = Outputs(:, [end-2:end 1:end-3]);

clear outputs*