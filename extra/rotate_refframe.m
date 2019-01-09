function gdpar_rot = rotate_frame(full_trnfrm);
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
    [x_lab, z_lab] = arrayfun(@(x, z) rotation_transform(x,z,-theta), gd.x, gd.z);
else
    [x_lab, z_lab] = arrayfun(@(x, z) rotation_transform(x,z,-theta), ...
        squeeze(gd.x(:,1,:)), squeeze(gd.z(:,1,:)));
end
z_lab = z_lab + Lx * sind(theta);   % set bottom of tank to zero

% create new grid-parameter structure
gd.x = x_lab;
gd.z = z_lab;
gdpar_rot.gd = gd;
gdpar_rot.params = params;
