function c = temperature(m, gamma)
%  TEMPERATURE  divergent colour map with nice red and blue end values
%   TEMPERATURE(M), is an M-by-3 matrix that defines a colormap.
%   The width of the gradient around the middle colour can be
%   changed by adjusting gamma.
%   The colors begin with cool blue, range through shades of
%   blue to white, and then through shades of red to a beautiful red.
%   TEMPERATURE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%   Gamma specifies the strength of gradient around the white middle 
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
    gamma = 1; % exponential value
elseif nargin == 1
    gamma = 1; % exponential value
elseif nargin > 2
    error('Too many input arguments. Only one is expected.')
end

% Nice blue
rb = 10;
gb = 23;
bb = 139;

% Nice red
rr = 128; % 193
gr = 6;   % 5
br = 9;   % 6

% From nice blue to white,
% then white to nice red
if (mod(m,2) == 0)
    m1 = m*0.5+1;
else
    m1 = ceil(m*0.5);
end
ll = linspace(0, 1, m1);
lr = linspace(1, 0, m1);
r1 = (1-rb/255)*ll'.^gamma + rb/255; 
g1 = (1-gb/255)*ll'.^gamma + gb/255; 
b1 = (1-bb/255)*ll'.^gamma + bb/255; 
r2 = (1-rr/255)*lr'.^gamma + rr/255;
g2 = (1-gr/255)*lr'.^gamma + gr/255;
b2 = (1-br/255)*lr'.^gamma + br/255;
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

