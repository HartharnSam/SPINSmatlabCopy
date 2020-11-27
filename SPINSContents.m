%% List of key SPINS files

% initialread    -  Script to call for initial processing of a run: 
%                        - Process diagnostics
%                        - Calculate Wave Statistics
%                        - Plot movie of density, u, w, vorticity
% spins_params
%% Case Setup Scripts
% Scripts to help produce spins.conf files
%
% initial_condition_maker - Like hillreader, but also sorts out density
%                       and tracer fields 
%
% density_profile_comparisons - Plots initial density field behind gate,
%                       and in main tank to compare with probes
%
%% Initial Analysis Scripts
% Scripts to produce wave metrics
% characterize_wave     -   Script to calculate wave characteristics
%                           (amplitude, wavelengths, time to hit slope,
%                           speed)
% plot_wave_char        -   Script to plot how wave characteristics from
%                           characterize_wave change over time
% plot_diagnos          -   Plot diagnostics recorded in model, including
%                           change of energy, max speeds and vorticities
% plot_stress           -   Plot stresses recorded in model. 
% plot_waveform         -   Plot showing waveform for comparison with
%                               DJL/lab
%% Plotting scripts
% basic_single_plot     - produces simplest single frame plot of real 
%                           density and horizontal velocity.
% plot_waveform
% SPINS_movie_maker     - Main script to plot the wave propagating over a
%                           fixed window
% make_movie            - SPINS built movie maker, which tracks the wave,
%                          and can be plotted with any spins_plot2d
%                          compatable variable
% spins_plot2d          - Plots a range of variables, including things like
%                           KE, Ri, , fairly complex options 
% 


%% In depth analysis




