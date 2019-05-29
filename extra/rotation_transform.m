function [x_new, z_new] = rotation_transform(x, z, theta)
    % Rotate the the point (x,z), or an array of points,
    % counter-clockwise by some angle theta (in degrees)
    %
    % David Deepwell, 2019

    if ~isscalar(x)
        [x_new, z_new] = arrayfun(@(xx,zz) rotate_element(xx, zz, theta), x, z);
    else
        [x_new, z_new] = rotate_element(x, z, theta);
    end
end

function [x_new, z_new] = rotate_element(x, z, theta)
    R = [cosd(theta) -sind(theta); ...
         sind(theta) cosd(theta)];
    v = [x; z];
    v_new = R*v;
    x_new = v_new(1);
    z_new = v_new(2);
end
