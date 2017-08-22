function [max_val, pos, ind] = find_wave_max(x, y)
% FIND_WAVE_MAX     Find the positions at which the variable y has local maximums.
%
%  Usage:
%    [max_val, pos, ind] = find_wave_max(gd.x, rho_contour)
%
%  Inputs:
%    'x'   - horizontal component of contour
%    'y'   - height of component of contour
%
%  Outputs:
%    'max_val'  - value of the local maxima of the contour given by (x,y)
%    'pos'      - the interpolated positions giving the local maxima
%    'ind'      - the nearest index (of x or y) to these maxima
%
%  David Deepwell, 2016

% set parameters for findpeaks
[y_max, max_ind] = max(y);
min_height = y_max/3;
if max_ind < length(x)
    [y_sep, ind] = find_position(x(max_ind:end), y(max_ind:end), 0.5*y_max);
    if ~isreal(y_sep)
        y_sep = x(max_ind-1+ind);
    end
    min_pk_dist = (y_sep - x(max_ind))/3;
    % remove locations where grid doubles back
    ni = 1;
    xmax = x(end);
    inds(1) = length(x);
    for nn = inds(1)-1:-1:1
        if x(nn) < xmax
            ni = ni+1;
            xmax = x(nn);
            inds(ni) = nn;
        end
    end
    xc = x(fliplr(inds));
    yc = y(fliplr(inds));

    %x_diff = x(2:end) > x(1:end-1);
    %x = x(x_diff);
    %y = y(x_diff);
    % find peaks if x is monotonic
    if all(diff(xc)>0)
        [pks, locs, width, prom] = findpeaks(yc, xc, 'MinPeakHeight', min_height,...
            'MinPeakProminence',min_height/2,'NPeaks', 15, 'MinPeakDistance', min_pk_dist);
        %findpeaks(y, x, 'MinPeakHeight', min_height,...
        %                        'SortStr', 'descend', 'NPeaks', 5,...
        %                        'MinPeakDistance',min_pk_dist,...
        %                        'WidthReference', 'halfheight',...
        %                        'Annotate', 'extents');
    else
        pks = [];
    end

    if isempty(pks)
        % if not using findpeaks, just take the absolute maxima
        max_val = y_max;
        ind = nearest_index(y, y_max);
        pos = x(ind);
        warning('No peaks found, using largest value as peak.')
    else
        % improve location and max given by pks
        % by fitting three points near peak with a quadratic (y = ax^2 + bx + c)
        for ii = 1:length(pks)
            loc_ind = nearest_index(x, locs(ii));
            if loc_ind > 1 && loc_ind < length(x)
                x1 = loc_ind - 1;
                x2 = loc_ind;
                x3 = loc_ind + 1;
                inds = [x1 x2 x3];
                % make vectors of the points
                xs = x(inds);
                ys = y(inds);
                szx = size(xs);
                if szx(1) == 1
                    xs = xs';
                end
                szy = size(ys);
                if szy(1) == 1
                    ys = ys';
                end
                % do interpolation
                mat = [xs.^2, xs ones(3,1)];
                quad_fit = mat\ys;
                a = quad_fit(1);
                b = quad_fit(2);
                c = quad_fit(3);
                % find max value and position
                max_val(ii) = -b^2/4/a + c;
                pos(ii) = -b/2/a;
            else
                max_val(ii) = pks(ii);
                pos(ii) = locs(ii);
            end
            ind(ii) = loc_ind;
        end
    end

    % sort peaks to be from right to left
    if length(pks) > 1
        [pos, order] = sort(pos, 'descend');
        max_val = max_val(order);
        ind = ind(order);
    end
else
    % wave has reached the tank end
    max_val = 0;
    pos = x(end);
    ind = length(x);
end
