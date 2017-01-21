function [] = print_figure(fig_hand, filename, type) 
% PRINT_FIGURE    Prints the figure (given by fig_hand) 
%
%  Usage:
%    print_figure(gcf,'filename')
%
%  Inputs:
%    'fig_hand' - a figure handle
%    'filename' - the name of the file
%    'type'     - the file type
%
%  Outputs:
%    - none
%
% David Deepwepll, 2016

    % default file name
    if ~exist('type', 'var')
        type = 'png';
    if ~exist('filename', 'var')
        filename = 'fig';
    end

    % set-up figure size
    set(fig_hand, 'Units', 'Inches')
    pos = get(fig_hand, 'Position');
    set(fig_hand, 'PaperPositionMode','Auto',...
                  'PaperUnits', 'Inches',...
                  'PaperSize', [pos(3) pos(4)]);
    % save figure
    print(fig_hand, filename, ['-d',type], '-r500')
end
