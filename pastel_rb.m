function c = pastel_rb(m, gamma)
%  PASTEL_RB  divergent colour map with nice red and blue end values
%   PASTEL_RB(M, gamma), is an M-by-3 matrix that defines a colormap.
%   The width of the gradient around the middle colour can be
%   changed by adjusting gamma.
%   The colors begin with light blue, range through shades of
%   blue to black, and then through shades of red to a light red.
%   PASTEL_RB, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(pastel_rb)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT, DARKJET, TEMPERATURE.
%
%  David Deepwell, 2015

if nargin < 1
    m = size(get(gcf,'colormap'),1);
    gamma = 1; % exponential value
elseif nargin == 1
    gamma = 1; % exponential value
elseif nargin > 2
    error('Too many input arguments. Only one is expected.')
end

% Nice blue
rb = 181;
gb = 214;
bb = 255;

% Nice red
rr = 255;
gr = 155;
br = 137;

% From nice blue to black to nice red
if (mod(m,2) == 0)
    m1 = m*0.5+1;
else
    m1 = ceil(m*0.5);
end
lr = linspace(0, 1, m1);
ll = linspace(1, 0, m1);
r1 = (rb/255)*ll'.^gamma;
g1 = (gb/255)*ll'.^gamma;
b1 = (bb/255)*ll'.^gamma;
r2 = (rr/255)*lr'.^gamma;
g2 = (gr/255)*lr'.^gamma;
b2 = (br/255)*lr'.^gamma;
if (mod(m,2) == 0)
    r = [r1(1:end-1); r2(2:end)];
    g = [g1(1:end-1); g2(2:end)];
    b = [b1(1:end-1); b2(2:end)];
else
    r = [r1(1:end); r2(2:end)];
    g = [g1(1:end); g2(2:end)];
    b = [b1(1:end); b2(2:end)];
end

c = [r g b]; 

