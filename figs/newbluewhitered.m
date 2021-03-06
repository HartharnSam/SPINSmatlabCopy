function newmap = newbluewhitered(m, pivot)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.


if nargin < 1
   m = size(get(gcf,'colormap'),1);
end


bottom = [0 0 0.5];
middle = [1 1 1];
top = [0.5 0 0];

RGB_weight = [.299, .587, .114]';

% Find middle
lims = get(gca, 'CLim');
if nargin<2
    pivot = 0;
end

% Find ratio of negative to positive
if (lims(1) < pivot) & (lims(2) > pivot)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = (pivot - lims(1)) / diff(lims);
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    % Just negative
    new = [bottom; middle];
    len = size(new, 1);
    
    %Fix around brightness
    map_brightness = sqrt(new*RGB_weight);

    newsteps = linspace(map_brightness(1), map_brightness(2), neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(map_brightness, new(:,i), newsteps)', 0), 1);
        
    end
    map_brightness_2 = sqrt(newmap1*RGB_weight);
    for i = 1:3
        % Interpolate over RGB spaces of colormap again to become
        % perpetually uniform with more sample points
        newmap1(:, i) = min(max(interp1(map_brightness_2, newmap1(:, i), newsteps)', 0), 1);
    end
    
    % Just positive
    new = [middle; top];
    len = length(new);
    map_brightness = sqrt(new*RGB_weight);
    newsteps = linspace(map_brightness(1), map_brightness(2), poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(map_brightness, new(:,i), newsteps)', 0), 1);
        
    end
    map_brightness_2 = sqrt(newmap*RGB_weight);
    for i = 1:3
        % Interpolate over RGB spaces of colormap again to become
        % perpetually uniform with more sample points
        newmap(:, i) = min(max(interp1(map_brightness_2, newmap(:, i), newsteps)', 0), 1);
    end
    
    % And put 'em together
    newmap = [newmap1; newmap];
    
elseif lims(1) >= pivot
    % Just positive
    new = [middle; top];
    len = length(new);
    map_brightness = sqrt(new*RGB_weight);
    newsteps = linspace(map_brightness(1), map_brightness(2), m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(map_brightness, new(:,i), newsteps)', 0), 1);
        
    end
    map_brightness_2 = sqrt(newmap*RGB_weight);
    for i = 1:3
        % Interpolate over RGB spaces of colormap again to become
        % perpetually uniform with more sample points
        newmap(:, i) = min(max(interp1(map_brightness_2, newmap(:, i), newsteps)', 0), 1);
    end    
else
    % Just negative
    new = [bottom; middle];
    len = length(new);
    map_brightness = sqrt(new*RGB_weight);
    newsteps = linspace(map_brightness(1), map_brightness(2), m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(map_brightness, new(:,i), newsteps)', 0), 1);
        
    end
    map_brightness_2 = sqrt(newmap*RGB_weight);
    for i = 1:3
        % Interpolate over RGB spaces of colormap again to become
        % perpetually uniform with more sample points
        newmap(:, i) = min(max(interp1(map_brightness_2, newmap(:, i), newsteps)', 0), 1);
    end
    
end
% 
% m = 64;
% new = [bottom; botmiddle; middle; topmiddle; top];
% % x = 1:m;
% 
% oldsteps = linspace(0, 1, 5);
% newsteps = linspace(0, 1, m);
% newmap = zeros(m, 3);
% 
% for i=1:3
%     % Interpolate over RGB spaces of colormap
%     newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
% end
% 
% % set(gcf, 'colormap', newmap), colorbar