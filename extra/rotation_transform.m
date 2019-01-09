function [x_new, z_new] = rotation_transform(x, z, theta)
% Rotate the the point (x,z) clockwise by some angle theta (in degrees)
% Vectors and Arrays don't work, use arrayfun with rotation_transform
%
% David Deepwell, 2019

    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    v = [x; z];
    v_new = R*v;
    x_new = v_new(1);
    z_new = v_new(2);
end
