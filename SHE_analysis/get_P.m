function p = get_P(u, du_dx, du_dz, dw_dx, dw_dz, rho)
% This is meant to calculate the pressure field using the poisson equation
% for pressure. I'm really not confident of it, so it is not production
% value. 
% A first check may be comparing it to djles_pressure, which IS
% production value!
% Poisson equation:
% nabla^2 p = rho_0((du/dx)^2 + (dw/dz)^2 + 2 * du/dz * dw/dx + ...
%    nu * nabla^2 (nabla . u) ) - d2rho/dz2*g


%% Read in variables from spins.conf
params = spins_params;
nu = params.visco; 
g = params.g;
rho0 = params.rho_0;
dx = params.Lx/params.Nx;

% % Calculate the laplacians
lap_u = get_grad2(du_dx) + get_grad2(du_dz);
lap_w = get_grad2(dw_dx) + get_grad2(dw_dz);

lap2p = rho0 * (du_dx.*du_dx + 2 * du_dz.*dw_dx + dw_dz .* dw_dz + ...
    nu*(lap_u+lap_w)) - get_grad2(rho)*g;

% Solve the Poisson equation for pressure (p)
% Preallocate the pressure field
p = zeros(size(u));

% Parameters for relaxation solver
maxIter = 5000;
tolerance = 1e-4;
omega = 1;  
% Relaxation factor, omega
% omega < 1 underrelaxation is stable, but slow 
% omega = 1 the Gauss-Seidel method is guaranteed if problem well posed
% omega > 1 over-relaxation is fast but can become unstable

% Iterate to solve the Poisson equation using relaxation
for iter = 1:maxIter
    p_old = p;
    
    % Relaxation loop
    p(2:end-1, 2:end-1) = (1-omega) * p(2:end-1, 2:end-1) + omega * 0.25 * ...
                           (p(3:end, 2:end-1) + p(1:end-2, 2:end-1) + ...
                            p(2:end-1, 3:end) + p(2:end-1, 1:end-2) - ...
                            lap2p(2:end-1, 2:end-1) * dx^2);
    p(:, 1) = 0;    % Bottom boundary (z = 0)
    p(:, end) = 0;  % Top boundary (z = Lz)
    % Check for convergence
    if max(max(abs(p - p_old))) < tolerance
        disp(['Pressure Poisson solver converged in ' num2str(iter) ' iterations']);
        break;
    end
end

if iter == maxIter
    disp('Warning: Pressure solver did not converge');
end