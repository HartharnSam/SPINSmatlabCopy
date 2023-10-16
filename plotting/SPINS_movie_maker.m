function spins_movie_maker(parameters, xlimits, min_t, savevideo, savefnm)
%SPINS_MOVIE_MAKER - Plots a video of the SPINS outputs outlined by
% parameters over the time and x space specified
%
% Syntax:  SPINS_movie_maker(parameters, xlimits, min_t, savevideo, savefnm)
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
%               bot_stresses
%               timeseries_dissipation
%               timeseries_tracer
%
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
%x = x-params.L_adj;
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

% if exist('diagnostics.mat', 'file') == 2
%     load diagnostics.mat Max_u Max_w
%     umaxabs = max(Max_u);
%     wmaxabs = max(Max_w);
%
% end
try load('wave_characteristics.mat', 'WaveStats', 'wave_center')
    wave_speed = WaveStats.meanWaveSpeed;
catch
    warning('No WaveStats saved, run characterize_wave if normalised velocity fields requested')
end

%% Identify which parameters to be used & their positions
n_plots = length(parameters);
if n_plots > 2
    n_columns = 2;
else
    n_columns = 1;
end
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
positioning.bot_stresses = find(strcmp(parameters, 'bot_stresses'));
positioning.timeseries_dissipation = find(strcmp(parameters, 'timeseries_dissipation'));
positioning.timeseries_tracer = find(strcmp(parameters, 'timeseries_tracer'));
positioning.dye1 = find(strcmp(parameters, 'dye1'));
positioning.dye2 = find(strcmp(parameters, 'dye2'));

pcolor_plots = zeros(n_columns, n_rows);

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
set(gcf, 'Units', 'centimeters'); set(gcf, 'Position', [1 1 26 11]);
set(gcf, 'PaperUnits','centimeters'); set(gcf, 'Position', [1 1 26 11]);

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

xspac = round(1:params.Nx/128:length(xlimits_ind));
zspac = linspace(-1, 1, 8);

for yi = 1:length(zspac)
    yspaceinds(yi) = nearest_index(clencurt(params.Nz-1), zspac(yi));
end
zspac = yspaceinds;
Xuv = x(xlimits_ind(xspac), yspaceinds);
Zuv = z(xlimits_ind(xspac), yspaceinds);

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
bot_stresses = nan(length(x(:, 1)), length(time));
%% Run video code
ax = gobjects(1, n_plots);
for ii = min_t:max_t

    %% Rho
    if ~isempty(positioning.rho)
        rho = rho_converter(spins_reader_new('rho',ii, xlimits_ind, [])); % Read in density and convert to real units
        rholimits = [rho_0 rho_0*(1+delta_rho)];

        ax(positioning.rho) = subaxis(n_rows,n_columns, positioning.rho, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,rho),shading flat
        hold on
        plot(x(:, 1), z(:, end), 'k-');
        %colormap darkjet
        c = colorbar(ax(positioning.rho));
        ylabel(c, '$\rho (kg m^{-3})$');
        cmocean('dense');
        caxis(gca, rholimits);
        %title('Density')
        title(subplot_labels(positioning.rho, 2))
        pcolor_plots(positioning.rho) = 1;
    end
    %% U Velocity
    if ~isempty(positioning.u)
        u = spins_reader_new('u', ii, xlimits_ind, []);
        %umaxabs=max(abs(u(:)));
        ax(positioning.u) = subaxis(n_rows,n_columns, positioning.u, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,u),shading flat
        cmocean('delta');
        c = colorbar(ax(positioning.u));
        ylabel(c, '$u$');
        umaxabs = 0.15;
        caxis(gca, umaxabs*[-1 1]);
        title('u (m/s)')
        pcolor_plots(positioning.u) = 1;

    end
    %% Normalised U Velocity
    if ~isempty(positioning.u_normalised)
        u = spins_reader_new('u', ii, xlimits_ind, []);
        c = wave_speed(1);
        u_normalised = u./c;
        ax(positioning.u_normalised) = subaxis(n_rows, n_columns, positioning.u_normalised, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x, z, u_normalised), shading flat
        %colormap darkjet
        cmocean('delta');
        c = colorbar(ax(positioning.u_normalised));
        ylabel(c, '$u / c$');
        caxis(gca, [-.1 .1]);
        %title('u/c')
        title(subplot_labels(positioning.u_normalised, 2))
        pcolor_plots(positioning.u_normalised) = 1;

    end
    %% V Velocity
    if ~isempty(positioning.v)
        v = spins_reader_new('v', ii, xlimits_ind, []);
        %vmaxabs=max(abs(v(:)));

        ax(positioning.v) = subaxis(n_rows,n_columns, positioning.v, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,v),shading flat

        caxis(gca, vmaxabs*[-1 1]);
        title('v (m/s)')
        pcolor_plots(positioning.v) = 1;

    end
    %% W Velocity
    if ~isempty(positioning.w)
        w = spins_reader_new('w', ii, xlimits_ind, []);
        %wmaxabs=max(abs(w(:)));

        ax(positioning.w) = subaxis(n_rows,n_columns, positioning.w, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,w),shading flat
        colormap(ax(positioning.w), darkjet)
        cmocean('delta');
        wmaxabs = 0.15;
        caxis(gca, wmaxabs*[-1 1]);
        c = colorbar(ax(positioning.w));
        ylabel(c, '$w$');
        title('w (m/s)')
        pcolor_plots(positioning.w) = 1;

    end
    %% Normalised W Velocity
    if ~isempty(positioning.w_normalised)
        w = spins_reader_new('w', ii, xlimits_ind, []);
        c = wave_speed;
        w_normalised = w./c(1);

        ax(positioning.w_normalised) = subaxis(n_rows, n_columns, positioning.w_normalised, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x, z, w_normalised), shading flat
        %colormap darkjet
        cmocean('delta');
        c = colorbar(ax(positioning.w_normalised));
        ylabel(c, '$w / c$');
        caxis(gca, [-0.1 0.1]);
        %title('w/c')
        title(subplot_labels(positioning.w_normalised, 2))
        pcolor_plots(positioning.w_normalised) = 1;

    end
    %% Dissipation
    if ~isempty(positioning.diss)
        diss = spins_reader_new('diss', ii, xlimits_ind, []);
        disslimits = [-12 0];

        ax(positioning.diss) = subaxis(n_rows,n_columns, positioning.diss, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,log10(diss)),shading flat
        cmocean('amp');
        caxis(gca, disslimits);
        title('dissipation')
        pcolor_plots(positioning.diss) = 1;
        c = colorbar(ax(positioning.diss));
        ylabel(c, '$\epsilon (m^2s^{-3})$');
    end
    %% Vorticity
    if ~isempty(positioning.vorty)
        vorty = spins_reader_new('vorty', ii, xlimits_ind, []);
        u_small = spins_reader_new('u',ii, xlimits_ind(xspac), (zspac));
        w_small = spins_reader_new('w', ii, xlimits_ind(xspac), (zspac));

        ax(positioning.vorty) = subaxis(n_rows,n_columns, positioning.vorty, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,vorty),shading flat; hold on
        plot(x(:, 1), z(:,1), 'k-');
        plot(x(:, end), z(:,end), 'k-');

        caxis([-1 1]*6);
        colormap(gca, newbluewhitered);
        q_scale = 1.2;
        linespec = struct('Color', 'k', 'LineWidth', .2);
        vekplot2(Xuv, Zuv, u_small, w_small, q_scale, linespec);
        if ispc
            vekLeg('SouthWest', q_scale, .1, linespec);
        end
        c = colorbar(ax(positioning.vorty));
        ylabel(c, '$\omega (s^{-1})$');
        %title('Vorticity/Velocity')
        title(subplot_labels(positioning.vorty, 2))
        pcolor_plots(positioning.vorty) = 1;

    end
    %% Tracer
    if ~isempty(positioning.tracer)
        tracer = spins_reader_new('tracer', ii, xlimits_ind, []);

        ax(positioning.tracer) = subaxis(n_rows,n_columns, positioning.tracer, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,tracer),shading flat
        cmocean('turbid');
        caxis(gca, [0 1]);
        c = colorbar(ax(positioning.tracer));
        ylabel(c, 'Tracer');
        title(subplot_labels(positioning.tracer, 2))
        pcolor_plots(positioning.tracer) = 1;

    end

    %% Dye 1
    if ~isempty(positioning.dye1)
        dye1 = spins_reader_new('dye1', ii, xlimits_ind, []);

        ax(positioning.dye1) = subaxis(n_rows,n_columns, positioning.dye1, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,dye1),shading flat
        cmocean('turbid');
        caxis(gca, [0 1]);
        c = colorbar(ax(positioning.dye1));
        ylabel(c, 'Tracer 1');
        title(subplot_labels(positioning.dye1, 2))
        pcolor_plots(positioning.dye1) = 1;

    end
    if ~isempty(positioning.dye2)
        dye2 = spins_reader_new('dye2', ii, xlimits_ind, []);

        ax(positioning.dye2) = subaxis(n_rows,n_columns, positioning.dye2, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x,z,dye2),shading flat
        cmocean('turbid');
        caxis(gca, [0 1]);
        c = colorbar(ax(positioning.dye2));
        ylabel(c, 'Tracer 2');
        title(subplot_labels(positioning.dye2, 2))
        pcolor_plots(positioning.dye2) = 1;

    end
    %% Dissipation Timeseries
    if ~isempty(positioning.timeseries_dissipation)
        load('all_diagnos', 'all_diagnos')
        [~, timeind] = min(abs(ii - all_diagnos.EnergyBudget.Time));
        dissip = all_diagnos.EnergyBudget.KE2Int_tot(1:timeind);
        mixing = all_diagnos.EnergyBudget.APE2BPE_tot(1:timeind);
        ax(positioning.timeseries_dissipation) = subaxis(n_rows,n_columns, positioning.timeseries_dissipation, 'SpacingVert', SpacingVert, 'Margin', Margin);
        plot(all_diagnos.EnergyBudget.Time(1:timeind), dissip, 'k-');
        hold on
        plot(all_diagnos.EnergyBudget.Time(1:timeind), mixing, 'b-');
        xlim([min_t max_t]);
        ylim([min(all_diagnos.EnergyBudget.KE2Int_tot) max(all_diagnos.EnergyBudget.KE2Int_tot)])

        xlabel('Time (s)'); ylabel('Energy (J)');
        c = colorbar;
        ax_pos = get(ax(positioning.timeseries_dissipation), 'Position');
        delete(c);
        set(ax(positioning.timeseries_dissipation),'Position', ax_pos);
        set(gca, 'XDir', 'normal');
        pcolor_plots(positioning.timeseries_dissipation) = 0;

    end
    %% Tracer Concentration Timeseries
    %     if ~isempty(positioning.timeseries_tracer)
    %         rho = spins_reader_new('rho', ii);
    %         x_start = nearest_index(z(:, 1), params.pyc_loc);
    %
    %         [tracer_volume(ii) weightings] = density_tracer(rho, x(x_start, 1));
    %
    %         ax(positioning.timeseries_tracer) = subaxis(n_rows,n_columns, positioning.timeseries_tracer, 'SpacingVert', SpacingVert, 'Margin', Margin);
    %         plot(1:ii, tracer_volume)
    %         xlim([min_t max_t]);
    %         ylim([0 sum(weightings(:))/2])
    %
    %         xlabel('Time (s)'); ylabel('Tracer beyond 6.5m ($m^3$)');
    %         c = colorbar;
    %         ax_pos = get(ax(positioning.timeseries_tracer), 'Position');
    %         delete(c);
    %         set(ax(positioning.timeseries_tracer),'Position', ax_pos);
    %         set(gca, 'XDir', 'normal');
    %         pcolor_plots(positioning.timeseries_tracer) = 0;
    %
    %     end
    %% Hovmoller of bed stress
    if ~isempty(positioning.bot_stresses)
        tmp_bot_stresses = calc_bot_stress(ii);
        bot_stresses(:, ii) = tmp_bot_stresses(xlimits_ind);

        ax(positioning.bot_stresses) = subaxis(n_rows,n_columns, positioning.bot_stresses, 'SpacingVert', SpacingVert, 'Margin', Margin);
        pcolor(x(:, 1), time, bot_stresses'); shading flat;
        cmocean('curl');
        caxis(ax(positioning.bot_stresses), [-.1 .1]);
        c = colorbar(ax(positioning.bot_stresses));
        ylabel(c, '$t_x$');
        title(subplot_labels(positioning.bot_stresses, 2))
        pcolor_plots(positioning.bot_stresses) = 1;
    end
    %% Correct axis titles etc.

    for i = 1:n_columns
        [~, last_in_column(i)] = find((pcolor_plots(i,:) == 1), 1, 'last');
    end
    for ai = 1:n_plots
        if pcolor_plots(ai)
            if mod(ai, n_columns) ~= 1 && n_columns ~= 1

                ax(ai).YLabel = [];
                ax(ai).YTickLabel = {};
            else
                ax(ai).YLabel.String = 'z (m)';
            end
            ax(ai).YLim = [params.min_z params.min_z+params.Lz];
           % set(ax(ai), 'XDir', 'reverse');
            ax(ai).XLim = xlimits;
            col_num = mod(ai-1, n_columns)+1;
            if ceil(ai/n_columns) == last_in_column(col_num)
                ax(ai).XLabel.String = 'x (m)';
            else
                ax(ai).XLabel = [];
                ax(ai).XTickLabel = {};
            end
        end
        %daspect(ax(ai), [1 1 1]);
    end

    if ~isempty(positioning.bot_stresses)
        ax(positioning.bot_stresses).YLim = [min_t max_t];
        ax(positioning.bot_stresses).YLabel.String = 't (s)';
        ax(positioning.bot_stresses).YTickLabelMode = 'auto';
        ax(positioning.bot_stresses).XLabel.String = 'x (m)';

    end
    drawnow;
    sgtitle(['t = ', num2str(ii), 's'])
    %dark_figure(gcf, [20 20 20]);
    if savevideo
        figure_print_format(gcf, 14);
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