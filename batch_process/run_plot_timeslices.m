%Fission_Lab_VsNumerics
clearvars; close all; clc;
orig_dir = cd; % Save initial folder
%% Settings for whole plot
output_fnm = 'F1_Fission'; % Filename for saved images (saves as both .png and .eps)
isFinal = true; % Switch to format and save the figure, false for a quick draft

c_range = .1*[-1 1]; % color axes limits on the background plot. 
% List of good options dependent on plotting variable:
%   vorty / vortyrho - 6*[-1 1];
%   rho - [1026 1046];
%   u / v / urho / vrho - .1*[-1 1];

%% Settings for left column
left_col.DirectoryName = ['./'];
left_col.Type = 'Model'; % 'Model' or 'Lab'
left_col.Field = {'Ri'}; % Field to plot - for model only, choice of :
% 'rho' (density color background only) 
% 'vortyrho' (vorticity color background overlaid by isopycnals) 
% 'vorty' (vorticity color background only) 
% 'ri' (Richardson number color background only) 
% 'urho' (horizontal velocity color background overlaid by isopycnals) 
% 'diss' (dissipation color background only) 
left_col.Times = [40 50 60 70]; % Time for each panel in LHS
left_col.x1 = [3 6]; % x limits for the left hand column

%% Settings for right column
right_col.DirectoryName = ['./'];
right_col.Type = 'Model';
right_col.Field = {'Fr'};
right_col.Times = [40 50 60 70];
right_col.x1 = [3 6];% x limits for right hand column

%% run the script/plot data
ColumnProperties{1} = left_col;
ColumnProperties{2} = right_col;

plot_timeslices(ColumnProperties, c_range, [1 2], [-0.405 0])
fig = gcf;
figaxes = findall(fig, 'Type', 'Axes');
for i = 1:size(figaxes, 1)
    daspect(figaxes(i), [1 1 1])
end

figure_print_format(fig);
set(gcf, 'PaperUnits', 'centimeters', 'Units', 'centimeters');
set(gcf, 'PaperPosition', [0, 0, 19, 15], 'Position', [0 0 19 15]);

if isFinal
    %figure_print_format(fig);
    cd(orig_dir);
    
    
    % IF YOUR MATLAB VERSION IS BEFORE R2020a UNCOMMENT THESE 2 LINES
    %print(['../Figures/', output_fnm, '.png'], '-dpng');
    %print(['../Figures/', output_fnm, '.eps'], '-depsc', '-opengl');

    % IF YOUR MATLAB VERSION IS BEFORE R2020a COMMENT THESE LINES
    exportgraphics(fig, ['../Figures/', output_fnm, '.png'], 'Resolution', 300);
    exportgraphics(fig, ['../Figures/', output_fnm, '.eps']);
end
