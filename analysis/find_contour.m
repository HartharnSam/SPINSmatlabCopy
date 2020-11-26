function [cont_x, cont_y] = find_contour(x, y, field, val, cont_ind)
% FIND_CONTOUR     Find the contour where field == val.
%
%  Usage:
%    [cont_x, cont_y] = find_contour(gd.x, gd.z, rho, rho_0)
%
%  Inputs:
%    'x'        - the x-grid
%    'y'        - the y-grid
%    'field'    - the field given on the grid
%    'val'      - the contour value
%    'cont_ind' - contour index (optional argument)
%
%  Outputs:
%    'cont_x' - the x positions of the contour
%    'cont_y' - the y positions of the contour
%
% Example:
%    [cont_x, cont_y] = find_contour(x, y, rho, 1026); % would require
%    readjusted rho to be used
%
%
%  David Deepwell, 2016

    % set contour index to zero if not provided
    % (to use the longest contour)
    if nargin < 5
        cont_ind = 0;
    end

    cont = contour_data(x, y, field, [1 1]*val);
    % fix for when there are multiple contours in cont
    len_c = length(cont);
    if len_c == 1 % if a single contour
        cont_x = cont.X;
        cont_y = cont.Y;
    elseif len_c > 1 % if multiple contours
        if cont_ind == 0
            % find which contour has most elements (this is the default)
            cont_size = zeros(len_c,1);
            for kk = 1:len_c
                cont_size(kk) = cont(kk).Length;
            end
            [~, cont_ind] = max(cont_size);
        end
        % select that contour
        cont_x = cont(cont_ind).X;
        cont_y = cont(cont_ind).Y;
    else
        disp('No contours found for top isopycnal'), return
    end
end
