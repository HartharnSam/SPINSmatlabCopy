function [] = print_figure(fig_hand,filename) 
% PRINT_FIGURE    Prints the figure (given by fig_hand) 
%
%  Usage:
%    print_figure(gcf,'filename')
%
%  Inputs:
%    'fig_hand' - a figure handle
%    'filename' - the name of the file
%
%  Outputs:
%    - none
%
% David Deepwepll, 2016

    % default file name
    if nargin<2
        filename = 'test';
    end

    % set-up figure size
    set(fig_hand, 'Units', 'Inches')
    pos = get(fig_hand, 'Position');
    set(fig_hand, 'PaperPositionMode','Auto',...
                  'PaperUnits', 'Inches',...
                  'PaperSize', [pos(3) pos(4)]);
    % save figure
    print(fig_hand, filename, '-dpdf', '-r500')
end
