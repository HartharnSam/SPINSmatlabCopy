function [pos, ind] = find_position(x, y, val)
% FIND_POSITION     Find the position at which the variable y is equal to 'val'.
%
%  Usage:
%    [pos, ind] = find_position(gd.x, rho, rho_0)
%
%  Inputs:
%    'x'   - a monotonically increasing vector (typically a grid vector)
%    'y'   - a monotonic vector to search for where y=val
%    'val' - a number (should be in range of y)
%
%  Outputs:
%    'pos' - the interpolated position at which y=val
%    'ind' - the nearest index (of x or y) to where y=val
%
%  David Deepwell, 2016

    ind = nearestindex(y, val);
    if val < min(y) || val > max(y)
        pos = x(ind);
    else
        % find points for interpolation
        if val > y(ind) && ind < length(x)
            ind1 = ind;
            ind2 = ind + 1;
        elseif val < y(ind) && ind > 1
            ind1 = ind -1;
            ind2 = ind;
        end
        x1 = x(ind1);
        x2 = x(ind2);
        y1 = y(ind1);
        y2 = y(ind2);
        % do interpolation
        pos = x1 + (val-y1)*(x2-x1)/(y2-y1);
    end
end
