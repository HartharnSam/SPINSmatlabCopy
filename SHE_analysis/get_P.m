function p = get_P(U, adv_u, buoy_u, rot_u, dt, idx_u, idx_w, Dx, Dz)
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
%nu = params.visco; 
g = params.g;
%rho0 = params.rho_0;
dx = params.Lx/params.Nx;

% % Calculate the laplacians
%lap_u = get_grad2(du_dx) + get_grad2(du_dz);
%lap_w = get_grad2(dw_dx) + get_grad2(dw_dz);

%du_dx = get_grad2(u);
%[~, dw_dz] = get_grad2(w);

%if params.Ny == 1
%    dv_dy = du_dx *0;
%else
%    error('3-D undefined');
%end
U(:, 1, :) = 0; 
U(:, end, :) = 0;
lap2p = -get_div(adv_u, idx_u, idx_w, Dx, Dz) + get_div(buoy_u, idx_u, idx_w, Dx, Dz) ...
    + get_div(rot_u, idx_u, idx_w, Dx, Dz) + 1/dt *get_div(U, idx_u, idx_w, Dx, Dz);
%lap2p = get_div(U, idx_u, idx_w);
%du_dx + dv_dy + dw_dz;
%rho0 * (du_dx.*du_dx + 2 * du_dz.*dw_dx + dw_dz .* dw_dz + ...
    %nu*(lap_u+lap_w)) - get_grad2(rho)*g;

% Solve the Poisson equation for pressure (p)
% Preallocate the pressure field
u = U(idx_u{:});
p = zeros(size(u));

% Parameters for relaxation solver
maxIter = 5000;
tolerance = 1e-6;
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

end

function f_grad = get_div(F, idx_u, idx_w, Dx, Dz)
dFx_dx = Dx*F(idx_u{:});
dFz_dz = (Dz*(F(idx_w{:})'))';
f_grad = dFx_dx+dFz_dz;


end