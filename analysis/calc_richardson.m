function Ri = calc_richardson(ii)
%CALC_BOT_STRESS Calculates bottom stresses of mapped SPINS output
% Offline version of bottom_stress_x.cpp for use in MATLAB
% Sam Hartharn-Evans 
% 09/08/2021

[x, z] = spinsgrid2d;
u = spins_reader_new('u', ii);
rho = spins_reader_new('rho', ii);
params = spins_params;
g = params.g;
dx = diff(x(:, 1));

% get -du/dx
[~, dudz] = get_grad2(u);

% get dw/dz 
[~, drhodz] = get_grad2(rho);

N_squared = -(g/params.rho_0).*drhodz;

Ri = N_squared./(dudz.^2);

end