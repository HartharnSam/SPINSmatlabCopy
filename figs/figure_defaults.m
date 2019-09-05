function [] = figure_defaults(fig_hand)
% FIGURE_DEFAULTS       Makes the figure to be prettier (more like Latex).
%  Run this function after making the plot to improve it's look.
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
    if nargin == 0
        fig_hand = gcf;
    end

    % settings
    fontsize = 12;

    % get figure children
    childs = allchild(fig_hand);

    % Change axis
    ax = findobj(childs,'Type','Axes');
    if ~isempty(ax)
        for ii = 1:length(ax)
            if ~strcmp(ax(ii).Type, 'axestoolbar')
                % change font size and font
                ax(ii).FontSize = fontsize;
                ax(ii).TickLabelInterpreter = 'Latex';
                ax(ii).XLabel.Interpreter = 'Latex';
                ax(ii).YLabel.Interpreter = 'Latex';

                % change axis outline
                set(ax(ii),'layer','top');
                ax(ii).Box = 'on';

                % remove line colour on contourf plots
                for jj = 1:length(ax(ii).Children)
                    if strcmp(ax(ii).Children(jj).Type, 'contour') ...
                            && strcmp(ax(ii).Children(jj).Fill, 'on')
                        ax(ii).Children(jj).LineStyle = 'none';
                    end
                end
            end
        end
    end

    % Change colorbar
    cbar = findobj(childs,'Type','ColorBar');
    if ~isempty(cbar)
        for ii = 1:length(cbar)
            cbar(ii).TickLabelInterpreter = 'Latex';
            cbar(ii).Label.Interpreter = 'Latex';
        end
    end

    % Change legend
    leg = findobj(childs,'Type','Legend');
    if ~isempty(leg)
        for ii = 1:length(leg)
            leg(ii).Interpreter = 'Latex';
        end
    end

    % Change text size
    text_hand = findall(fig_hand, 'Type', 'Text');
    %texthands = texthands(~strcmp(get(texthands,'string'),''));
    set(text_hand, 'Interpreter', 'Latex',...
                   'Fontsize',fontsize)
end
