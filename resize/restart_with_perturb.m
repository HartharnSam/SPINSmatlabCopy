function restart_with_perturb(u_prime, w_prime, rho_prime, intended_time, new_dir)
%RESTART_WITH_PERTURB - Outputs various fields with added perturbation
%ready for restarting with e.g. wave_reader.cpp
%
% Inputs:
%    u_prime - Perturbation field on u
%    w_prime - " " 
%    rho_prime - " "
%    intended_time - Simulation time for initial restart time (in s)
%
% Other m-files required: get_output_times, nearest_index,
% spins_reader_new, spins_writer
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% Department of Geography & Environmental Sciences, Northumbria University
% email address: sam.hartharn-evans@northumbria.ac.uk
% GitHub: https://github.com/HartharnSam
% 17-Mar-2025; Last revision: 17-Mar-2025
% MATLAB Version: 23.2.0.2668659 (R2023b) Update 9

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%intended_time = 210;

[times, outputs] = get_output_times;
ii = outputs(nearest_index(times, intended_time));

u = mean(spins_reader_new('u', ii));
w = mean(spins_reader_new('w', ii));
rho = mean(spins_reader_new('rho', ii));
[x, z] = spinsgrid2d;

% Example for periodic perturbation in velocity:
%eps = 1e-3; 
%kx = 3;
%kz = 3;

%alpha = kx * 2*pi;
%beta  = kz * 2*pi;
%u_prime = eps.*sin(alpha*x).*cos(beta*z);
%w_prime = (-eps.*cos(alpha*x).*sin(beta*z));

u = u + u_prime;
w = w + w_prime;
rho = rho + rho_prime;

figure
subplot(3, 2, 1);
pcolor(x, z, u);
subplot(3, 2, 2);
pcolor(x, z, u-mean(u));

subplot(3, 2, 3);
pcolor(x, z, w);
subplot(3, 2, 4);
pcolor(x, z, w-mean(w));

subplot(3, 2, 5);
pcolor(x, z, rho);
subplot(3, 2, 6);
pcolor(x, z, rho-mean(rho));


mkdir("../"+new_dir)
cd("../"+new_dir);
spins_writer('u2d', u);
spins_writer('w2d', w);
spins_writer('rho2d', rho);
spins_writer('x2d', x);
spins_writer('z2d', z);