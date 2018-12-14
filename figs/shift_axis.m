function [] = shift_axis(xshift, yshift)
% SHIFT_AXES    Shift the axis horizontally and vertically
%
% David Deepwell, 2018

% get position
pos = get(gca, 'Position');
pos(1) = pos(1) + xshift;
pos(2) = pos(2) + yshift;
% set position
set(gca, 'Position', pos);
