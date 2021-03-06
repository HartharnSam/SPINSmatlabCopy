function gdpar_rot = rotate_refframe(full_trnfrm);
% rotate grid into laboratory reference frame
%
% David Deepwell, 2019
global gdpar
split_gdpar

% shorten some parameters
Lx    = params.Lx;
theta = params.tilt_angle;

if nargin == 0
    full_trnfrm = true;
end

% apply the rotation transformation
if full_trnfrm
    [x_lab, z_lab] = rotation_transform(gd.x, gd.z, -theta);
else
    [x_lab, z_lab] = rotation_transform(squeeze(gd.x(:,1,:)), squeeze(gd.z(:,1,:)), -theta);
end
z_lab = z_lab + Lx * sind(theta);   % set bottom of tank to zero

% create new grid-parameter structure
gd.x = x_lab;
gd.z = z_lab;
gdpar_rot.gd = gd;
gdpar_rot.params = params;
