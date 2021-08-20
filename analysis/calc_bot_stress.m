function stress_x = calc_bot_stress(ii)
%CALC_BOT_STRESS Calculates bottom stresses of mapped SPINS output
% Offline version of bottom_stress_x.cpp for use in MATLAB
% Sam Hartharn-Evans 
% 09/08/2021

[x, z] = spinsgrid2d;
u = spins_reader_new('u', ii);
w = spins_reader_new('w', ii);
params = spins_params;
dx = diff(x(:, 1));
Hprime = gradient(z(:, 1), dx(1));
% get -du/dx
[dudx, dudz] = get_grad2(u);

% get dw/dz 
[dwdx, dwdz] = get_grad2(w);

stress_x = 2*Hprime .* (-dudx + dwdz);

temp = (1-(Hprime.^2)).*(dudz + dwdx);
stress_x = stress_x + temp;

coeff =params.visco./(sqrt(1+(Hprime.^2)));

stress_x = params.visco.*params.rho_0./(1+Hprime.^2) .*stress_x;
stress_x = stress_x(:, 1);
end