function vekLeg(Location, scale, arrow_value, line,mark)
%VEKLEG - Produces a "legend" like scale arrow as produced by vekplot2.m
%
% Syntax: vekLeg(Location, scale)
%
% Inputs: 
%   Location - Location of the vector legend in format like legend
%               'Location'
%   scale    - Length of arrow in the figure when u^2 + v^2 = 1. Should be
%               same as vekplot2 scale
%   arrow_value - Arrow representing this velocity (in m/s)
%
% Outputs:
%   N/A
% 
% Example:
%   vekLeg('NorthEast', 1, 1)
% 
% Other m-files required: vekplot2.m 
% Subfunctions: n/a
% MAT-files required: n/a
% 
% See also: vekplot2.m

% Author: Sam Hartharn-Evans
% Newcastle University - s.hartharn-evans2@ncl.ac.uk
% April 2020; Last Revision 07/04/2020
% Version 1.0

% -------------------------------------------------------------------------
%% Check Variables Set
% -------------------------------------------------------------------------
if nargin < 4 || isempty(line)
    line = 'b';
end
if nargin <5 || isempty(mark)
    mark = '-';
end

% -------------------------------------------------------------------------
%% Begin Code
% -------------------------------------------------------------------------
% Open a legend, record it's position and close it
f = gca; % Call current axis, open one if not present
l = legend('Location', Location);
Position_leg = l.Position;
old_axesPos = f.Position; % Also record the revised position of the axis
                          % so that changes in this due to placement of 
                          % the legend can be retained
delete(l); % Delete the legend again

f.Position = old_axesPos; % Reposition main axes as if a legend had been placed there


% Calculate width of figure in x units
f_xunit_width = diff(f.XLim); % Width of plot in data units
f_axes_width = f.InnerPosition(3); % Width of plot in figure units
%xUnitChangeInFigure = f_xunit_width/f_axes_width; % Change in data units for 1 change in figure units
DataToFigure = f_axes_width/f_xunit_width;
FigureToData = f_xunit_width/f_axes_width; % 

% Calculate required width of box
t = text(6 , .15, strcat([num2str(arrow_value), ' m/s']), 'Horizontalalignment', 'center');

box_width_dataUnits = max(scale*arrow_value, t.Extent(3));
box_width = box_width_dataUnits*DataToFigure*1.3;

%% Produce the vector legend
g = axes('InnerPosition', [Position_leg(1), Position_leg(2), box_width, Position_leg(4)*.5], 'Box', 'on'); % Plot in the space previously identified by the legend creation

relative_height = Position_leg(4)/f.Position(4)*diff(f.YLim);
vekplot2(0, relative_height/3, arrow_value, 0, scale, line, mark); % Produce an arrow 1 unit in size, scaled as per "scale"

AxesXLim = [0 box_width_dataUnits]-(box_width_dataUnits-scale*arrow_value)/2; %scale*arrow_value]

xlim(g, AxesXLim);
ylim(g, [0 relative_height])
delete(t);
t = text(g, mean(AxesXLim) , relative_height*.66, strcat([num2str(arrow_value), ' m/s']), 'Horizontalalignment', 'center');
% Remove axis like formatting
xticks([]); yticks([]); box on
g.Position = [Position_leg(1), Position_leg(2), Position_leg(3), Position_leg(4)*.5];
% Restore current axes to the main axes just in case this is needed
set(gcf, 'CurrentAxes', f);
end