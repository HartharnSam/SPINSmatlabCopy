%Fission_Lab_VsNumerics
clearvars; close all; clc;
orig_dir = cd; % Save initial folder
isFinalfinal = false;
%% Settings for whole plot
time_spacing = 10; % dt between panels

output_fnm = 'F1_Fission'; % save filename
left_col.x1 = [10 12.5]-.6; % x limits for the left hand column
right_col.x1 = left_col.x1;% x limits for right hand column
%y1 = [0.0188 0.1020]; % y limits
isFinal = true; % Switch to format and save the figure, false for a quick draft
c_range = .1*[-1 1]; % color axes limits on the background plot

%% Settings for left column
left_col.DirectoryName = ['..\..\02_Raw_Data\Surface_10L_270720_33'];
left_col.Type = 'Model'; % 'Model' or 'Lab'
left_col.Field = {'urho'}; % Field to plot - for model only, choice of 'rho', 'vortyrho', 'vorty', 'ri'
left_col.Times = (0:4)*time_spacing +130; % Time for each panel in LHS

%% Settings for right column
right_col.DirectoryName = left_col.DirectoryName;
right_col.Type = 'Model';
right_col.Field = {'rho'};
right_col.Times = (0:4)*time_spacing +130;

%% run the script/plot data
ColumnProperties{1} = left_col;
ColumnProperties{2} = right_col;

plot_timeslices(ColumnProperties, c_range, [1 2], [-0.3 0])
fig = gcf;
figaxes = findall(fig, 'Type', 'Axes');
for i = 1:size(figaxes, 1)
    daspect(figaxes(i), [2 1 1])
end

figure_print_format(fig);
set(gcf, 'PaperUnits', 'centimeters', 'Units', 'centimeters');
set(gcf, 'PaperPosition', [0, 0, 19, 15], 'Position', [0 0 19 15]);

if isFinal
    %figure_print_format(fig);
    cd(orig_dir);
    
    set(gcf, 'PaperUnits', 'centimeters', 'Units', 'centimeters');
    set(gcf, 'PaperPosition', [0, 0, 19, 15], 'Position', [0 0 19 15]);
    
    %print(['../../05_Output/PaperFigures_Raster/', output_fnm, '.png'], '-dpng');
    exportgraphics(fig, ['../../04_Output/PaperFigures_Raster/', output_fnm, '.png'], 'Resolution', 300);
    if isFinalfinal
        exportgraphics(fig, ['../../04_Output/PaperFigures_Vector/', output_fnm, '.eps'], 'ContentType', 'vector');
    else
        exportgraphics(fig, ['../../04_Output/PaperFigures_Vector/', output_fnm, '.eps']);
    end
    %print(['../../05_Output/PaperFigures_Vector/', output_fnm, '.eps'], '-depsc', '-opengl');
end
