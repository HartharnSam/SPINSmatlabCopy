function [cont_x, cont_y] = find_contour(x, y, field, val)
% FIND_CONTOUR     Find the contour where field=val.
%
%  Usage:
%    [cont_x, cont_y] = find_contour(gd.x, gd.z, rho, rho_0)
%
%  Inputs:
%    'x'     - the x-grid
%    'y'     - the y-grid
%    'field' - the field given on the grid
%    'val'   - the contour value
%
%  Outputs:
%    'cont_x' - the x positions of the contour
%    'cont_y' - the y positions of the contour
%
%  David Deepwell, 2016

    cont = contour_data(x, y, field, [1 1]*val);
    % fix for when there are multiple contours in cont
    len_c = length(cont);
    if len_c == 1 % if a single contour
        cont_x = cont.X;
        cont_y = cont.Y;
    elseif len_c > 1 % if multiple contours
        % find which contour has most elements (it's likely the proper one)
        cont_size = zeros(len_c,1);
        for kk = 1:len_c
            cont_size(kk) = cont(kk).Length;
        end
        [~, cont_ind] = max(cont_size);
        % select that contour
        cont_x = cont(cont_ind).X;
        cont_y = cont(cont_ind).Y;
    else
        disp('No contours found for top isopycnal'), return
    end
end
