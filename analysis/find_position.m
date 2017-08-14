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

    ind = nearest_index(y, val);
    if val < min(y) || val > max(y)
        pos = x(ind);
    else
        % find points for quadratic interpolation
        if ind > 1 && ind < length(x)
            ind1 = ind - 1;
            ind2 = ind;
            ind3 = ind + 1;
        elseif ind == 1
            ind1 = ind;
            ind2 = ind + 1;
            ind3 = ind + 2;
        elseif ind == length(x)
            ind1 = ind - 2;
            ind2 = ind - 1;
            ind3 = ind;
        end
        inds = [ind1 ind2 ind3];
        % do interpolation
        xs = x(inds)';
        ys = y(inds);
        mat = [xs.^2, xs, ones(3,1)];
        quad_fit = mat\ys;
        a = quad_fit(1);
        b = quad_fit(2);
        c = quad_fit(3);
        pos1 =  sqrt(val/a + (b/2/a)^2 - c/a) - b/2/a;
        pos2 = -sqrt(val/a + (b/2/a)^2 - c/a) - b/2/a;
        % find which root is the proper one
        posn = [pos1 pos2];
        dist = abs([1 1]*x(ind) - posn);
        [~,p_ind] = min(dist);
        pos = posn(p_ind);
    end
end
