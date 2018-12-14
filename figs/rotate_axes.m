function [] = rotate_axes(fig_hand, theta);
% Rotate the axes object clockwise by some angle theta
%
% David Deepwell, 2018

% get the correct figure
if nargin > 0
    figure(fig_hand); % use provided figure handle
end
% if not, use the current figure
ax = gca;

% If no angle specified, use the one in the spins.conf
if nargin < 2
    par = spins_params();
    theta = par.tilt_angle;
end

% rotate axes
UpVector = [-sind(theta), cosd(theta), 0];
DAR = get(ax, 'DataAspectRatio');
set(ax, 'CameraUpVector', DAR.*UpVector);
