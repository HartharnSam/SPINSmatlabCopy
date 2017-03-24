function Cout = contour_data(x, y, field, val)
% CONTOUR_DATA  Wrapper of CONTOUR to obtain a structured output
%
%  Usage:
%    S = contour_data(x, y, rho, 0.01);
%
%  Inputs:
%    'x'     - horizontal grid. Can be a vector or matrix
%    'y'     - vertical grid. Can be a vector or matrix
%    'field' - matrix of values at grid points
%    'val'   - value of contour
%
%  Outputs:
%    'Cout' - structure containing organized contours
%       Sub-elements are:
%         Level  - contour line value
%         Length - number of contour line points
%         X      - X coordinate array of the contour line
%         Y      - Y coordinate array of the contour line
%
%   See also contour.
%
% David Deepwell, 2017. This code is based off contourcs.m
% written by Takeshi Ikuma (see Mathworks.com).

% adjust if val is a single number
if length(val) == 1
    val = [1 1]*val;
end

% make contour
fig_hand = figure('Visible', 'Off');
[c_data, ~] = contour(x, y, field, val);
close(fig_hand)

% Count the number of contour segments
c_num = 0;  % number of contours
ii = 1;
while ii <= size(c_data,2)
    c_num = c_num + 1;
    ii = ii + c_data(2,ii) + 1;
end

% fill the output structure with contour data
ii = 1;
for nn = 1:c_num
    idx = (ii+1):(ii+c_data(2,ii));
    Cout(nn).Level  = c_data(1, ii);
    Cout(nn).Length = c_data(2, ii);
    Cout(nn).X =      c_data(1, idx);
    Cout(nn).Y =      c_data(2, idx);
    % re-order if necessary
    if Cout(nn).X(end) < Cout(nn).X(1)
        Cout(nn).X = fliplr(Cout(nn).X);
        Cout(nn).Y = fliplr(Cout(nn).Y);
    end
    ii = idx(end) + 1; % next starting index
end
