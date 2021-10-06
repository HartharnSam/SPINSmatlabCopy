%% List of key SPINS files

% initialread    -  Script to call for initial processing of a run: 
%                        - Process diagnostics
%                        - Calculate Wave Statistics
%                        - Plot movie of density, u, w, vorticity
% spins_params  -   Outputs structure of parameters passed in spins.conf
%                   file

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
%% Analysis
% calc_bot_stress.m     
% calc_richardson.m     
% clean_diagnostics.m   
% find_position.m       
% calc_buoyancy_freq.m
% characterize_wave.m
% find_contour.m
% find_wave_max.m
% FiniteDiff.m
% calc_reynolds.m
% chebyshev_volumes.m   
% find_half_max.m       

%% CaseSetup
% density_profile_comparisons.m
% initial_condition_maker.m  

%% Extra
% fftfreq.m 
% get_planar_grid.m       
% nearest_index.m          
% subaxis.m               
% filters.m
% get_plot_points.m
% par2var.m
% vekLeg.m
% check_make_dir.m
% find_fields.m
% get_rectilinear_grid.m
% parseArgs.m
% vekplot2.m
% clear_except.m
% first_output.m
% get_vector_grid.m
% rotate_refframe.m
% wave_region.m
% completion.m
% get_fixed_z.m
% interp_onto_rect_grid.m
% rotation_transform.m
% contour_data.m
% get_lab_refframe.m
% last_output.m
% spins_plotoptions.m      
% eqn_of_state.m
% get_output_times.m
% max_grid_spacing.m
% split_gdpar.m            

%% figs
% figure_defaults.m
% print_figure.m          
% figure_print_format.m   
% rgb.m                   
% betterplots.m           
% inferno.m               
% rotate_axes.m           
% bluewhitered.m          
% magma.m                 
% shift_axis.m            
% choose_caxis.m          
% newbluewhitered.m       
% subplot_labels.m        
% cmocean.m               
% plasma.m                
% temperature.m           
% crop_cmap.m             
% plotboxpos.m            
% tick_chooser.m          
% darkjet.m               
% print_all.m             
% viridis.m               
% default_line_colours.m  
% print_clr_bkgrd.m       

%% IO

%% plotting

%% resize





