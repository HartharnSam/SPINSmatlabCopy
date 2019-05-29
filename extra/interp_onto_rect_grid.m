function [gd_rect, field_rect] = interp_onto_rect_grid(gd, field)
%  INTERP_ONTO_RECT_GRID  Interpolate the field onto a rectilinear grid
%
%  Usage:
%    [gd_rect, field_rect] = interp_onto_rect_grid(field)
%
%  Inputs:
%    field      - Matrix of a data field
%
%  Outputs:
%    gd_rect    - Rectilinear grid
%    field_rect - Field on the rectilinear grid
%
%  David Deepwell, 2019

% find rectilinear grid
gd_rect = get_rectilinear_grid(gd);

% Interpolate field onto the rectilinear grid
F = scatteredInterpolant(gd.x(:), gd.z(:), field(:),'linear','none');
vq = F(gd_rect.x(:), gd_rect.z(:));
field_rect = reshape(vq, size(gd_rect.x));

% remove points outside initial domain
% (scatteredInterpolant makes a convex shape, not a concave one)
Nx = size(gd.x, 1);
Nx_rect = size(gd_rect.x, 1);
if Nx == Nx_rect
    % remove points above initial domain
    above_top = gd_rect.z > gd.z(:,end);
    field_rect(above_top) = NaN;
    % remove points below initial domain
    below_bot = gd_rect.z < gd.z(:,1);
    field_rect(below_bot) = NaN;
end
