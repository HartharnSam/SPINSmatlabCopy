function stress_x = calc_top_stress(ii)
%CALC_TOP_STRESS Calculates top stresses of mapped SPINS output
% Offline version of top_stress_x.cpp for use in MATLAB
% Sam Hartharn-Evans 
% 09/08/2021

[x, z] = spinsgrid2d;
u = spins_reader_new('u', ii);
w = spins_reader_new('w', ii);
params = spins_params;
dx = diff(x(:, 1));
Hprime = gradient(z(:, end), dx(1));
% get -du/dx
[dudx, dudz] = get_grad2(u);

% get dw/dz 
[dwdx, dwdz] = get_grad2(w);

stress_x = 2*Hprime .* (-dudx + dwdz);

temp = (1-(Hprime.^2)).*(dudz + dwdx);
stress_x = stress_x + temp;

stress_x = params.visco.*params.rho_0./(1+Hprime.^2) .*stress_x;
stress_x = stress_x(:, end);
end