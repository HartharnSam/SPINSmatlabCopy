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
%    'val'   - a number to find contour level
%
%  Outputs:
%    'cont_x' - the x positions of the contour
%    'cont_y' - the y positions of the contour
%
%  David Deepwell, 2016

    cont = contourcs(x, y, field, [1 1]*val);
    % fix for when there are multiple contours in cont
    if length(cont) == 1 % if a single contour
        cont_x = cont.X;
        cont_y = cont.Y;
    elseif length(cont) > 1 % if multiple contours
        % find which contour has most elements (it's likely the proper one)
        cont_size = zeros(length(cont),1);
        for kk = 1:length(cont)
            cont_size(kk) = length(cont(kk).Y);
        end
        [~, cont_ind] = max(cont_size);
        cont_x = cont(cont_ind).X;    % select that contour
        cont_y = cont(cont_ind).Y;
    else
        disp('No contours found for top isopycnal'), return
    end
end
