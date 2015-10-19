function c = pastel_pg(m, gamma)
%  PASTEL_PG  divergent colour map with pastel purple and green end values
%   PASTEL_PG(M, gamma), is an M-by-3 matrix that defines a colormap.
%   The width of the gradient around the middle colour can be
%   changed by adjusting gamma.
%   PASTEL_PG, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(pastel_pg)
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

% starting colour
rb = 182;
gb = 255;
bb = 157;

% end colour
rr = 231;
gr = 157;
br = 255;

% From starting colour to black to end colour
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

