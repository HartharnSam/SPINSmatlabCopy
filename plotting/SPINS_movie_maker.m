function SPINS_movie_maker(parameters, xlimits, min_t, savevideo, savefnm)
%SPINS_MOVIE_MAKER - Plots a video of the SPINS outputs outlined by
%parameters over the time and x space specified
%
% Syntax:  SPINS_movie_maker(savevideo,parameters, xlimits, min_t, savefnm)
%
% Inputs:
%    parameters - [REQUIRED] List of fields to be plotted - in order of
%                 subplots
%           options:
%               rho
%               u
%               u_normalised
%               w
%               w_normalised
%               v
%               diss
%               vorty
%               tracer
%    xlimits - [OPTIONAL] Horizontal region to be plotted - can be either [lowerx higherx]m (in wcs), or 
%              'slopeonly' to plot the slope (as identified within the
%              wavestats slopehit parameter. Defaults to whole tank. If set
%              to [], goes to default
%    savevideo - [OPTIONAL] Boolean if movie will save or not. Default
%              false
%    min_t - [OPTIONAL] Start time of animation. Default 0
%    savefnm - [OPTIONAL] Filename of save file. Default 'v2.mp4'
%
% Example:
%    SPINS_movie_maker(true, {'rho', 'vorty', 'u', 'w'}, 'slopeonly', WaveStats.endSlope);
%    SPINS_movie_maker(false, {'rho', 'u'}, [4.5 7], 15);
%
% Other m-files required: spins_params, spinsgrid2d, rho_converter,
% spins_reader_new, subaxis, bluewhitered
% Subfunctions: none
% MAT-files required: none
%
% See also: simple_time_pic_example
% Author: Sam Hartharn-Evans
% School Maths Stats and Physics, Newcastle University, UK
% email address: s.hartharn-evans2@newcastle.ac.uk
% August 2020; Last revision: 28-Sept-2020
% MATLAB Version: 9.8.0.1323502 (R2020a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%% Load initial global SPINS parameters
params = spins_params;

rho_0 = params.rho_0;
delta_rho = params.delta_rho;

[x, z] = spinsgrid2d;
x = x-params.L_adj; 
if nargin<2 || isempty(xlimits)
    xlimits = [min(min(x)) max(max(x))];
end
if nargin<4
    savevideo = false;
end
if nargin<5
    if ispc
        savefnm = 'v2.mp4';
    else
        savefnm = 'v2.avi';
    end
end

if exist('diagnostics.mat', 'file') == 2
    load diagnostics.mat Max_u Max_w
    umaxabs = max(Max_u);
    wmaxabs = max(Max_w);
    
end
try load('wave_characteristics.mat', 'WaveStats', 'wave_center')
    wave_speed = WaveStats.meanWaveSpeed;
catch
    warning('No WaveStats saved, run characterize_wave if normalised velocity fields requested')
end

%% Identify which parameters to be used & their positions
n_plots = length(parameters);
n_columns = 2;
n_rows = ceil(n_plots/n_columns);
% List of possible parameters, identifies the order
positioning.rho = find(strcmp(parameters, 'rho'));
positioning.u = find(strcmp(parameters, 'u'));
positioning.u_normalised = find(strcmp(parameters, 'u_normalised'));
positioning.v = find(strcmp(parameters, 'v'));
positioning.w_normalised = find(strcmp(parameters, 'w_normalised'));
positioning.w = find(strcmp(parameters, 'w'));
positioning.diss = find(strcmp(parameters, 'diss'));
positioning.vorty = find(strcmp(parameters, 'vorty'));
positioning.tracer = find(strcmp(parameters, 'tracer'));

%% Begin Video and set figure settings
if savevideo % Open videowriter class
    if ispc
        vid = VideoWriter(savefnm, 'MPEG-4');
    else
        vid = VideoWriter(savefnm);
    end
    vid.FrameRate = 1;
    open(vid);
end

% Set spacing between and outside of subplots
SpacingVert = 0.09;
Margin = 0.09;

% Open figure and set sizes
figure(1);
clf
set(gcf, 'Units', 'centimeters'); set(gcf, 'Position', [1 1 44 20]);
set(gcf, 'PaperUnits','centimeters'); set(gcf, 'Position', [1 1 44 20]);

if ispc
    sgtitle(params.name, 'Interpreter','none');
end
%% Set spatial extent of plots
if strcmp(xlimits, 'slopeonly')
    if isfield(params, 'hill_height')
        slope_length = (params.hill_height/params.hill_slope)+params.hill_end_dist;
    elseif isfield(params, 'ice_length')
        slope_length = params.ice_length;
    end
    
    inds_free_wave = find(x(:, 1)>params.L_adj*1.15 & x(:, 1)<(params.Lx - slope_length*1.15)); % Identifies region where wave on flat bed and outside of gate influence
    
    inds =  [inds_free_wave(round(length(inds_free_wave)/6)), ...
        inds_free_wave(round(length(inds_free_wave)/2)),...
        inds_free_wave(round((length(inds_free_wave)/6)*5.5))]; % identifies three key points over flat bed
    
    xlimits_ind = inds(3):params.Nx;
    xlimits = [x(inds(3), 1) x(params.Nx, 1)];
else
    xlimits_ind = find(x(:, 1)>xlimits(1) & x(:, 1)<xlimits(2));
end

xspac = params.Nx/128;
zspac = params.Nz/8;
Xuv = x(xlimits_ind(1:xspac:end), 1:zspac:params.Nz);
Zuv = z(xlimits_ind(1:xspac:end), 1:zspac:params.Nz); 

x = x(xlimits_ind, :);
z = z(xlimits_ind, :);
%% Set up time system
max_t = length(dir('u.*'))-1;
time = round(get_output_times());

if nargin<4
    if exist('wave_characteristics.mat', 'file') == 2
        %load('wave_characteristics.mat', 'wave_center')
        loc_ind = find(wave_center>xlimits(1));
        if ~isempty(loc_ind)
            min_t = time(loc_ind(1));
        else
            min_t = time(end); 
        end
    else
        min_t = 0;
    end
end

%% Run video code

for ii=min_t:max_t
    ax = gobjects(1, n_plots);
    %% Rho
    if ~isempty(positioning.rho)
        rho = rho_converter(spins_reader_new('rho',ii, xlimits_ind, [])); % Read in density and convert to real units
        rholimits = [rho_0 rho_0*(1+delta_rho)];
        
        ax(positioning.rho) = subaxis(n_rows,n_columns, positioning.rho, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,rho),shading flat
        colormap darkjet
        caxis(gca, rholimits);
        title('Density')
                
    end
    %% U Velocity
    if ~isempty(positioning.u)
        u = spins_reader_new('u', ii, xlimits_ind, []);
        %umaxabs=max(abs(u(:)));
        
        ax(positioning.u) = subaxis(n_rows,n_columns, positioning.u, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,u),shading flat
        colormap darkjet
        caxis(gca, umaxabs*[-1 1]);
        title('u (m/s)')
        
    end
    %% Normalised U Velocity
    if ~isempty(positioning.u_normalised)
        u = spins_reader_new('u', ii, xlimits_ind, []);
        c = wave_speed;
        u_normalised = u./c;
        ax(positioning.u_normalised) = subaxis(n_rows, n_columns, positioning.u_normalised, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x, z, u_normalised), shading flat
        colormap darkjet
        caxis(gca, [-1.1 1.1]);
        title('u/c')
        
    end
    %% V Velocity
    if ~isempty(positioning.v)
        v = spins_reader_new('v', ii, xlimits_ind, []);
        %vmaxabs=max(abs(v(:)));
        
        ax(positioning.v) = subaxis(n_rows,n_columns, positioning.v, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,v),shading flat
        colormap darkjet
        caxis(gca, vmaxabs*[-1 1]);
        title('v (m/s)')
        
    end
    %% W Velocity
    if ~isempty(positioning.w)
        w = spins_reader_new('w', ii, xlimits_ind, []);
        %wmaxabs=max(abs(w(:)));
        
        ax(positioning.w) = subaxis(n_rows,n_columns, positioning.w, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,w),shading flat
        colormap darkjet
        caxis(gca, wmaxabs*[-1 1]);
        title('w (m/s)')
        
    end
        %% Normalised W Velocity
    if ~isempty(positioning.w_normalised)
        w = spins_reader_new('w', ii, xlimits_ind, []);
        c = wave_speed;
        w_normalised = w./c;
        
        ax(positioning.w_normalised) = subaxis(n_rows, n_columns, positioning.w_normalised, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x, z, w_normalised), shading flat
        colormap darkjet
        caxis(gca, [-1.1 1.1]);
        title('w/c')
    end
    %% Dissipation
    if ~isempty(positioning.diss)
        diss = spins_reader_new('diss', ii, xlimits_ind, []);
        disslimits = [-12 0];
        
        ax(positioning.diss) = subaxis(n_rows,n_columns, positioning.diss, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,log10(diss)),shading flat
        colormap darkjet
        caxis(gca, disslimits);
        title('dissipation')
        
        
    end
    %% Vorticity
    if ~isempty(positioning.vorty)
        vorty = spins_reader_new('vorty', ii, xlimits_ind, []);
        u_small = spins_reader_new('u',ii, xlimits_ind(1:xspac:end), (1:zspac:params.Nz));
        w_small = spins_reader_new('w', ii, xlimits_ind(1:xspac:end), (1:zspac:params.Nz));
        
        ax(positioning.vorty) = subaxis(n_rows,n_columns, positioning.vorty, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,vorty),shading flat; hold on
        if isfield(params, 'hill_height')
            plot([params.Lx, params.Lx-params.hill_height/params.hill_slope]-params.L_adj, params.min_z+[params.hill_height 0], 'k-', 'linewidth', .2);
        end
        caxis([-1 1]*6);
        colormap(gca, bluewhitered);
        q_scale = 1.2;
        linespec = struct('Color', 'k', 'LineWidth', .2);
        vekplot2(Xuv, Zuv, u_small, w_small, q_scale, linespec);
        if ispc
            vekLeg('SouthWest', q_scale, .1, linespec);
        end
        title('Vorticity/Velocity')
        
    end
    %% Tracer
    if ~isempty(positioning.tracer)
        tracer = spins_reader_new('tracer', ii, xlimits_ind, []);
        
        ax(positioning.tracer) = subaxis(n_rows,n_columns, positioning.tracer, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,tracer),shading flat
        colormap darkjet
        caxis(gca, [0 1]);
        title('tracer concentration')
    end
    
    %% Correct axis titles etc.
    
    for ai = 1:n_plots
        if mod(ai, n_columns) ~= 1
            
            ax(ai).YLabel = [];
            ax(ai).YTickLabel = {};
        else
            ax(ai).YLabel.String = 'z (m)';
        end
        
        if ceil(ai/n_rows) ~= n_rows
            ax(ai).XLabel = [];
            ax(ai).XTickLabel = {};
        else
            ax(ai).XLabel.String = 'x (m)';
        end
        ax(ai).YLim = [params.min_z params.min_z+params.Lz];
        set(ax(ai), 'XDir', 'reverse');
        colorbar(ax(ai));
        ax(ai).XLim = xlimits;
        %daspect(ax(ai), [1 1 1]);
        
        
    end
    if savevideo
        F = getframe(gcf);
        writeVideo(vid, F);
    end
end
if savevideo
    close(vid)
end

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------