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
n_layers = 2;
if n_layers == 2
    isTwoLayer = true; % Switch for if two layer (affects isopycnal chosen for characterize wave)
else
    isTwoLayer = false;
end
RunDirectoryName = {'Surface_20L_Hyperviscous'}; %{'28_121120', '29_131120'   }; % List of directories to process files from
%run_directory_names('Continuous', "Full");
%runname = {'thin', 'broad', 'surface'};
nlayers = [2];
%% Set location for video file to be saved to
pathname = ['./'];
%pathname = ['D:\Sam\OneDrive - Newcastle University\Shared_Videos\Numerics\'];

%% Run loops
outputs = NaN(length(RunDirectoryName), 21);
for ii = 1:length(RunDirectoryName)
    n_layers = nlayers(ii);
    if n_layers == 2
        isTwoLayer = true; % Switch for if two layer (affects isopycnal chosen for characterize wave)
    else
        isTwoLayer = false;
    end
    cd(['../', RunDirectoryName{ii}]);

    params = spins_params;
    disp(params.name);
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
        %if exist('wave_characteristics.mat', 'file') == 0 % Has characterize_wave been run before for this experiment?
        if n_layers == 1
            WaveStats = wave_characteriser([10 45]); % Less robust script for use with continuous stratification
        else
            WaveStats = characterize_wave(isTwoLayer, [0 50]); % runs the script to find wavelengths, amplitudes, wave position, speed etc.
        end
        plot_wave_char; % Plots the outputs and saves them
        close all
        %else % If this has already been run before, just need to load it
        load('wave_characteristics.mat', 'WaveStats');
    end

    %% Clean/Plot diagnostics
    if exist('diagnostics.mat', 'file') == 0 % checks for file created by clean_diagnostics
        clean_diagnostics(); % Cleans diagnostics files
    end

    if exist('all_diagnos.mat', 'file') == 0% checks for file created by plot_diagnos
        diagnos = plot_diagnos(isPlotDiagnostics, false, isPlotDiagnostics);
        close all
        %        plot_stress(isPlotDiagnostics, isPlotDiagnostics); % Also plot stresses
    else
        close all

        diagnos = load('all_diagnos.mat');
        diagnos = diagnos.all_diagnos;
    end

    %% plot movie
    if isVideo
        if ispc
            SPINS_movie_maker({'rho', 'vorty', 'u_normalised', 'w_normalised'},...
                'slopeonly', WaveStats.endTank, true, 'output.mp4');

            %SPINS_movie_maker({'rho', 'vorty', 'u_normalised', 'tracer'}, 'slopeonly', WaveStats.endSlope, true, fullfile(pathname, [params.name, '.mp4']));
            %SPINS_movie_maker({'rho', 'vorty', 'u_normalised', 'w_normalised'}, 'slopeonly', WaveStats.endSlope, true, fullfile(pathname, ['numerics_example', '.mp4']));

        else
            SPINS_movie_maker({'rho', 'vorty', 'u_normalised', 'w_normalised'}, 'slopeonly', WaveStats.endSlope, true, fullfile(pathname, [params.name, '.avi']));
        end
    end

    %% Calculate Reynolds statistics
    [maxRe, startRe(ii)] = calc_reynolds(n_layers);

    %% Collate stats to be pasted into appropriate index sheet
    slope_condition_tmp = num2str(round(params.hill_slope*1.5 *1000));
    if params.Lx>7
        slope_condition_tmp = [slope_condition_tmp, 'mm_b'];
    else
        slope_condition_tmp = [slope_condition_tmp, 'mm'];
    end
    slope_condition{ii} = slope_condition_tmp;

    %% Clean up
    clearvars -except is* RunDirectoryName pathname outputs* n_layers slope_condition start* nlayers runname
    close all
end
outputs_names = {'HillHeight', 'HillSlope', 'PycAdjLoc', 'pyc_thickness', ...
    'h1', 'amplitude', 'speed', 'wavelength', 'Wave Steepness', 'Ir', 'Classification', 'empty ',...
    'Nx', 'Ny', 'Nz', 'Lx', 'Ly', 'Lz', 'LxNx', 'LyNy', 'LzNz'};

Outputs = array2table(outputs, 'VariableNames', outputs_names);
% Outputs.FileNames = RunDirectoryName;
% Outputs.SlopeCondition = slope_condition;
% Outputs.InitialWaveCond = zeros(length(RunDirectoryName), 0);
% Outputs = Outputs(:, [end-2:end 1:end-3]);

clear outputs*
