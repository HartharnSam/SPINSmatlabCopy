function [] = figure_defaults(fig_hand)
% FIGURE_DEFAULTS       Makes the figure (given by fig_hand) to be prettier (more latex)
%
%  Usage:
%    figure_defaults(gcf)
%
%  Inputs:
%    'fig_hand' - a figure handle
%
%  Outputs:
%    - none
%
% David Deepwell, 2016

    % default figure
    if ~exist(fig_hand,'var')
        fig_hand = gcf;
    end

    % settings
    fontsize = 10;

    % get figure children
    childs = allchild(fig_hand);

    % Change axis
    ax = findobj(childs,'Type','Axes');
    if ~isempty(ax)
        ax.TickLabelInterpreter = 'Latex';
        ax.FontSize = fontsize;
        ax.XLabel.Interpreter = 'Latex';
        ax.YLabel.Interpreter = 'Latex';
    end

    % Change colorbar
    cbar = findobj(childs,'Type','ColorBar');
    if ~isempty(cbar)
        cbar.TickLabelInterpreter = 'Latex';
    end

    % Change legend
    leg = findobj(childs,'Type','Legend');
    if ~isempty(leg)
        leg.Interpreter = 'Latex';
    end
end
