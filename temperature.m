function c = temperature(m)
%  TEMPERATURE  divergent colour map with nice red and blue end values
%   TEMPERATURE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with cool blue, range through shades of
%   blue to white, and then through shades of red to a beautiful red.
%   TEMPERATURE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(temperature)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT, DARKJET.
%
%  David Deepwell, 2015

if nargin < 1
    m = size(get(gcf,'colormap'),1);
elseif nargin > 1
    error('Too many input arguments. Only one is expected.')
end

c = crop_cmap('balance',0.1,m);
