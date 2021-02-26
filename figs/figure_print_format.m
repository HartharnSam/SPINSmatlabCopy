function figure_print_format(fig_hand)
%FIGURE_PRINT_FORMAT - Formats the figure nicely ready for printing
%Changes to LaTeX like text, at appropriate text size, sets a nice line
%width, box outline, ticks, etc.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans, based on SPINSmatlab figure_defaults (David
% Deepwell, 2016)
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 15-Dec-2020; Last revision: 15-Dec-2020
% MATLAB Version: 9.9.0.1467703 (R2020b)

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
            
            % Get rid of the annoying exponent thing
            for i = 1:length(ax(ii).YAxis)
                ax(ii).YAxis(i).Exponent = 0;
            end
            for i = 1:length(ax(ii).XAxis)
                ax(ii).XAxis(i).Exponent = 0;
            end
            
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

%% Set some other things
settings = struct('DefaultLineLineWidth',1, 'DefaultTextFontWeight','normal',...
    'DefaultAxesFontWeight','normal', 'defaultfigurecolor', 'w',...
    'DefaultAxesTickDir', 'in');

fields = fieldnames(settings);
for i = 1:length(fields)
    str = strjoin(['settings.', fields(i),';'], '');
    set(fig_hand, fields{i}, eval(str));
end




