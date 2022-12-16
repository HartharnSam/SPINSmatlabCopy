function plot_timeslices(ColumnProperties, c_range, labelled_columns, ylims)
%PLOT_LAB_VS_NUMERICS - Plots panelled figure showing lab vs numerics at a chosen
%number of comparable times.
%
% Syntax:  plot_LabVsNumerics(ColumnProperties,output_fnm,input3)
%
% Inputs:
%    ColumnProperties - MATLAB Structure containing all the variables
%    needed, for each column. See fission_plotter.m for example of setup
%    outputFilename - Name of file where figure saved to
%    y1 - Y axis limits
%    c_range - Limits of the c axis
%    isFinal - If final true (default false), formats plot properly and saves image output
%    labelled_columns - Numbers of columns to have a colourbar associated
%                       with them
%
% Outputs:
%
% Other m-files required: newbluewhitered, cmocean, betterplots,
% figure_defaults, spinsgrid2d, spins_params, subaxis, subplot_labels,
% find_contour, spins_reader_new, rho_converter, dfireadvel, plotPIV
%
% Subfunctions: plot_model, plot_lab
% MAT-files required: none
%
% See also:
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 26-Feb-2021; Last revision: 26-Feb-2021
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
% Set up variables
n_columns = length(ColumnProperties);
n_times = length(ColumnProperties{1}.Times);
if nargin < 4
    ylims = [];
end

if nargin < 2
    c_range = 6*[-1 1];
    labelled_columns = n_columns;
end
if ~isfield(ColumnProperties, 'clrmap')
    clrmap = 'newbluewhitered'; % Colormap for the background image
else
    clrmap = ColumnProperties.clrmap;
end
contour_clrmap = brighten(cmocean("grey", 6), -.5); % Colours of the contour lines
fig = figure;
tlayout = tiledlayout(n_times, n_columns,'TileSpacing','compact', ...
    'Padding','compact'); %, 'TileIndexing', 'columnmajor');
betterplots(fig);
%% Load in data

orig_dir = cd; % Save initial folder
for jj = 1:n_columns
    s{jj} = gobjects(1, length(ColumnProperties{jj}.Times));
    for ii = 1:length(ColumnProperties{jj}.Times)
        plot_num = ((n_columns)*(ii-1))+(jj);
        s{jj}(ii) = nexttile(plot_num);
    end
    
    if strcmpi(ColumnProperties{jj}.Type, 'Model')
        cd(ColumnProperties{jj}.DirectoryName);
        
        plot_model(ColumnProperties{jj}, jj, clrmap, c_range, contour_clrmap, labelled_columns, s{jj}, ylims);
    else
        warning('Type must be Model');
    end
    drawnow;
end


cd(orig_dir);

end


function plot_model(ColumnProperties, column_number, ~, c_range, contour_clrmap, labelled_columns, s, ylims)
%% Load in model key parameters
[x, z] = spinsgrid2d; % Load SPINS grid
params = spins_params;
rho_0 = params.rho_0;
delta_rho = params.delta_rho;
x = x -.3 - params.L_adj; % Adjust to tank coordinates system

%% Crop x and z grids
xlimits = ColumnProperties.x1;
xinds = find(x(:, 1)>xlimits(1) & x(:, 1)<xlimits(2));
x = x(xinds, :);
z_midline = z(params.Nx/2,:); % Take a profile at the mid-tank before cropping
z = z(xinds, :);
z = z*100;
%% Pre-allocate the axes array
n_times = length(ColumnProperties.Times);
%s = gobjects(n_times, 1);
% Sort out timings
[actual_times, actual_otpt_nums] = get_output_times(false);
actual_times = round(actual_times,1);
times = round(ColumnProperties.Times, 1);
%iis = ColumnProperties.Times;%actual_otpt_nums(any(actual_times == round(ColumnProperties.Times, 2), 2));
% for i = 1:length(times)
%     iis(i) = actual_otpt_nums(actual_times == times(i));
% end
iis = times;

%% Start plotting
for i=1:n_times
    ii = iis(i);
    %% Set up Axis
    title(s(i), subplot_labels(n_times.*(column_number-1) +i  , 2));
    s(i).TitleHorizontalAlignment = 'left';
    %axes(s(i));
    %% Identify which field to plot
    if length(ColumnProperties.Field)>1
        Field = ColumnProperties.Field{i};
    else
        Field = ColumnProperties.Field{1};
    end
    %% Plot
    switch lower(Field)
        case 'rho' % Plot density field
            rho = spins_reader_new('rho',ii, xinds, []); %read in density
            rho = rho_converter(rho); % convert to real densities
            params = spins_params;
            pcolor(s(i), x, z, rho-1000),shading(s(i),'flat')
            caxis(s(i), [params.rho_0 params.rho_0.*(1+params.delta_rho)]-1000);
            colormap(s(i), cmocean('dense'))
            if any(column_number == labelled_columns) && i == 1
                c = colorbar(s(i));
                c.Location = 'northoutside';
                ylabel(c, '$\sigma_t (kg m^{-3})$')
            end
        case 'vorty' % plot vorticity field
            vorty = spins_reader_new('vorty', ii, xinds, []); % Read in data
            
            pcolor(s(i), x, z, vorty), shading(s(i),'flat'); % Plot data
            c_range = [-1 1].*6;
            caxis(s(i), c_range); % Sort color range
            colormap(s(i), cmocean('balance'))

            if any(column_number == labelled_columns) && i  == 1
                c = colorbar(s(i));
                c.Location = 'northoutside';
                ylabel(c, '$\omega (s^{-1})$')
                set(c, 'YTick', [min(c_range) 0 max(c_range)]);
            end
            
        case 'vortyrho' % Plot Vorticity (background) with isopycnals
            vorty = spins_reader_new('vorty', ii, xinds, []); % Read in vorticity
            rho = spins_reader_new('rho',ii, xinds, []); %read in density
            rho = rho_converter(rho); % convert to real densities
            pcolor(s(i), x, z, vorty), shading(s(i),'flat'); % Plot background
            caxis(s(i), c_range);
            colormap(s(i), cmocean('balance'))

            hold(s(i), 'on')
            %% Add Contour on top
            strat = spins_reader_new('rho', 0, params.Nx/2, 1, []); % Read a single vertical profile
            % Calculate the contour positions
            if params.h_halfwidth < 0.01
                isopyc_loc = params.pyc_loc + (params.h_halfwidth .* [-1.5 -.5 0 .5 1.5]);
                contour_clrmap = [contour_clrmap; flip(contour_clrmap)];
                contour_clrmap = contour_clrmap([4 5 6 7 8 9], :);
            elseif params.h_halfwidth < .07
                isopyc_loc = params.pyc_loc - params.h_halfwidth + (0.015:0.01125:.07);
                
            else
                isopyc_loc = params.pyc_loc + linspace(-params.h_halfwidth/2, params.h_halfwidth/2, 5);
                contour_clrmap = [flip(contour_clrmap); (contour_clrmap)];
                contour_clrmap = contour_clrmap([4 5 6 7 8 9], :);
                
            end
            contval = rho_converter(interp1(z_midline, strat, isopyc_loc));
            contval = contval(~isnan(contval));
            % Plot the contour
            for nn = 1:length(contval)
                [cont_x, cont_y] = find_contour(x, z, rho, contval(nn));
                if ~isempty(cont_x)
                    %plot(cont_x, cont_y, '-', 'Color', contour_clrmap(nn,:), 'LineWidth', .5);%);
                    plot(s(i), cont_x, cont_y, '-', 'Color', [.4 .4 .4], 'LineWidth', .5);%);
                end
            end
            % Sort the colourbar if needed
            if any(column_number == labelled_columns) && i == 1
                %originalSize = get(s(i), 'Position');
                c = colorbar(s(i));
                c.Location = 'northoutside';
                ylabel(c, '$\omega (s^{-1})$')
                %set(s(i), 'Position', originalSize);
                set(c, 'YTick', [min(c_range) 0 max(c_range)]);
            end
        case 'urho'
            u = spins_reader_new('u', ii, xinds, []); % Read in vorticity
            rho = spins_reader_new('rho',ii, xinds, []); %read in density
            rho = rho_converter(rho); % convert to real densities
            pcolor(s(i), x, z, u), shading(s(i),'flat'); % Plot background
            caxis(s(i), c_range);
            newbluewhitered(255, 0, s(i));
            
            hold(s(i), 'on');
            %% Add Contour on top
            strat = spins_reader_new('rho', 0, params.Nx/2, 1, []); % Read a single vertical profile
            % Calculate the contour positions
            if params.h_halfwidth < 0.01
                isopyc_loc = params.pyc_loc + (params.h_halfwidth .* [-1.5 -.5 0 .5 1.5]);
                contour_clrmap = [contour_clrmap; flip(contour_clrmap)];
                contour_clrmap = contour_clrmap([4 5 6 7 8 9], :);
            elseif params.h_halfwidth < .07
                isopyc_loc = params.pyc_loc - params.h_halfwidth + (0.015:0.0125:.075);
                
            else
                isopyc_loc = params.pyc_loc + linspace(-params.h_halfwidth/2, params.h_halfwidth/2, 5);
                contour_clrmap = [flip(contour_clrmap); (contour_clrmap)];
                contour_clrmap = contour_clrmap([4 5 6 7 8 9], :);
                
            end
            contval = rho_converter(interp1(z_midline, strat, isopyc_loc));
            contval = contval(~isnan(contval));
            % Plot the contour
            for nn = 1:length(contval)
                [cont_x, cont_y] = find_contour(x, z, rho, contval(nn));
                if ~isempty(cont_x)
                    %plot(cont_x, cont_y, '-', 'Color', contour_clrmap(nn,:), 'LineWidth', .5);%);
                    plot(s(i), cont_x, cont_y, '-', 'Color', [.4 .4 .4], 'LineWidth', .5);%);
                end
            end
            % Sort the colourbar if needed
            if any(column_number == labelled_columns) && i == 1
                %originalSize = get(s(i), 'Position');
                c = colorbar(s(i));
                ylabel(c, 'u (ms$^{-1})$')
                c.Location = 'northoutside';

                %set(s(i), 'Position', originalSize);
                set(c, 'YTick', [min(c_range) 0 max(c_range)]);
            end
        case 'ri' % Plot richardson number
            % Read in Ri
            try
                rho_z = spins_reader_new('rho_z',ii, xinds, [])';
                u_z   = spins_reader_new('u_z',ii, xinds, [])';
                
                g = params.g;
                if params.delta_rho < 1
                    rho_0 = 1;
                else
                    rho_0 = params.rho_0;
                end
                % calculate Ri
                N_sq  = -g/rho_0*rho_z;
                data = (N_sq./u_z.^2)';
            catch
                
                data = SPINS_derivs('ri', ii, true);
                data = data(xinds, :);
            end
            axes(s(i));
            ogdata = data;
            data = data.*NaN;
            data(ogdata>=0.25) = 2.5;
            data(ogdata< 0.25 & ogdata>0) = 1.5;
            data(ogdata<=0) = 0.5;
            contourf(s(i), x, z, data, [1 2]);
            shading flat;
            caxis(s(i), [0 3]);
            c = colorbar(s(i));
            colormap(s(i), cmocean('-dense', 3));
            if any(column_number == labelled_columns) && i == 1
                stops = linspace(0,3, 3+3+1); 
                c.YTick = stops(2:2:end-1);               
                c.Location = 'northoutside';
                c.TickLabels = {'$<0$', '$0-0.25$', '$>0.25$'};
                set(c, 'FontSize', 8);
                ylabel(c, '$Ri$');
            end
                  
            
        case 'diss'
            diss = spins_reader_new('diss', ii, xinds, []);
            pcolor(s(i), x,z,log10(diss)),shading(s(i),'flat')
            colormap(s(i), cmocean('amp'));
            caxis(s(i), [-7 -1]); % Sort color range
            if any(column_number == labelled_columns) && i == 1
                %originalSize = get(s(i), 'Position');
                c = colorbar(s(i));
                c.Location = 'northoutside';

                ylabel(c, '$log_{10}(\epsilon (m^2s^{-3}))$')
                %set(s(i), 'Position', originalSize);
            end
        case 'fr'
            if exist('KdV.mat', 'file') == 0 %
                calc_kdv_depthchange;
                error('re-run script')
            end
            load KdV.mat KdV
            % Calculate Fr = U/c_lw
            u = spins_reader_new('u',ii, xinds, []);
            data = u/KdV.c_0;
            axis(s(i));
            % remove data that is too large
            if max(data(:)) > 2
                data(data>2) = 2;
                warning(['Fr>2 has been set to 2.',...
                    'This enables contour and contourf to make meaningful plots.'])
            end
            %TODO _ Fix this and the colorbar
            % Plot data
            pcolor(s(i), x, z, data); shading(s(i), 'flat');
            hold(s(i), 'on');
            contour(s(i), x, z, data, [-.25 0 1 2], 'LineColor', 'none');
            clim(s(i), [-2 2]);
            colormap(s(i), cmocean('delta'));
                        
            if any(column_number == labelled_columns) && i == 1
                %originalSize = get(s(i), 'Position');
                c = colorbar(s(i));        
                c.Location = 'northoutside';
                set(c, 'FontSize', 8);
               % c.TickLabels = {'$<0$', '$0-1$', '$>1$'};
                ylabel(c, '$Fr$');
                %set(s(i), 'Position', originalSize);
            end
            
            
            
    end
    if ~(strcmpi(Field, 'Ri') || strcmpi(Field, 'rho') || strcmpi(Field, 'diss') || strcmpi(Field, 'Fr'))
        %eval(['colormap(s(i),' clrmap, ');']);
        %caxis(s(i), c_range);
    end
    %% Sort out axes
    % Y Axis & Colorbar
    if column_number ~= 1
        yticklabels(s(i), []);
    else
        ylabel(s(i), '$z (cm)$');
    end
    
    % X Axis
    if i == n_times
        xlabel(s(i),'$x (m)$')
    else
        xticklabels(s(i), []);
    end
    xlim(s(i), xlimits);
    if isempty(ylims)
        ylim(s(i), 'tight')
    else
        ylim(s(i), ylims);
    end
    
    daspect(s(i), [1 100 1]);
    tick_chooser('XTick', s(i));
    set(s(i), 'XDir', 'normal');
    %% Draw on slope, and outside
    hold(s(i), 'on')
    plot(s(i), x(:, 1), z(:, 1), 'k-', 'LineWidth', 1);
    plot(s(i), x(:, end), z(:, end), 'k-', 'LineWidth', 1);
    drawnow

    hold(s(i), 'off')
    box(s(i), 'off')
    set(s(i), 'TickDir', 'in');
    
    %% Add Time stamp
    %text(xlimits(1)*.1 + xlimits(2)*.9, -0.25, ['t = ', num2str(actual_times(actual_otpt_nums == ii)),' s'], 'FontSize', 14); % Writes the time of each plot
    
end

end
