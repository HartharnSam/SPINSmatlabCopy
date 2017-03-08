function c = crop_cmap(map_name, perc, varargin)
%  crop_cmap  function to crop the extrema of a cmocean colormap
%   CROP_CMAP, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(temperature)
%
%  David Deepwell, 2015

% defaults
NLevels = 128;
autozero = false;

% Does the user want to center a diverging colormap on zero?  
tmp = strncmpi(varargin,'zero',3);
if any(tmp)
    autozero = true;
    varargin = varargin(~tmp);
end

% Has user requested a specific number of levels? 
tmp = isscalar(varargin);
if any(tmp)
    NLevels = varargin{tmp};
end

% get colormap and it's indices
inds = 1:256;
if autozero
    cmap = cmocean(map_name, 'zero');
else
    cmap = cmocean(map_name);
end
% find new colormap
frst = perc*256;
last = (1-perc)*256;
new_inds = linspace(frst, last, NLevels);
for ii=1:3
    c(:,ii) = interp1(inds, cmap(:,ii), new_inds);
end
