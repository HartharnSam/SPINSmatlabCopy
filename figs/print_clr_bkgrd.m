function print_clr_bkgrd(bkgrd_clr, fname)
%PRINT_LIGHTONDARK - Changes plot background to be the colour, and all
%lines/axes text to be white. Saves using export_fig
%
% Inputs:
%    bkgrd_clr - Color of background, either as [R, G, B] or colour string
%    fname     - Filename to save to (with extension)

[~, ~, ext] = fileparts(fname); % Identify file type

% Collect axes handles (including colorbars)
gfigs = findall(gcf,'type','ColorBar', '-or', 'type', 'axes', '-or', 'type', 'title'); 

% Change axes properties
set(gfigs, 'YColor', 'w'); %Set x axis and labels etc to white
set(gfigs, 'XColor', 'w'); % set 
set(gcf, 'Color', bkgrd_clr);
% Export
export_fig(fname, ['-d', ext]);
