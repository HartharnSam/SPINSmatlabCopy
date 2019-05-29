function gdpar = get_lab_refframe(full_transform)
% Return the laboratory reference frame
% where the grid is in horizontal/vertical axis,
% not the along slope/across slope axis
%
% David Deepwell, 2019

% Optional arguments
if ~exist('full_transform','var')
    full_transform = true;
end

% get grid
global gdpar
gdpar = spins_gridparams('FastFull');

% rotate grid
if gdpar.params.ndims == 3
    if full_transform
        gdpar = rotate_refframe(true);
    else
        gdpar = rotate_refframe(false);
    end
else
    gdpar = rotate_refframe();
end
gdpar.params.rotated_grid = true;
