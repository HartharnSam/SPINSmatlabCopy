function [pos, ind] = find_half_max(x, y)
% FIND_HALF_MAX     Find the first location (while moving away from the location of maximum)
%    that y drops to half it's maximum. It assumes that the maximum is at the last index.
%
%  Usage:
%    [pos, ind] = find_half_max(cont_x, cont_y)
%
%  Inputs:
%    'x'   - horizontal component of contour
%    'y'   - height of component of contour
%
%  Outputs:
%    'pos'      - the interpolated position of first half maximum
%    'ind'      - the nearest index (of x or y) to this position
%
%  David Deepwell, 2016

    % set-up for loop
    len = length(y);
    ii = len+1;
    max_val = y(len);
    val = max_val;

    % loop from right to left
    while val > max_val/2 && ii>1
        ii = ii-1;
        val = y(ii);
    end

    if ii == 1 || ii >= len
        % didn't reach max_val/2, so take minimum in domain
        [~, ind] = min(y);
        pos = x(ind);
    else
        % interpolate to find location of half max
        inds = [ii-1 ii ii+1];
        [pos, ~] = find_position(x(inds), y(inds), max_val/2);
        ind = nearest_index(x, pos);
    end
end
