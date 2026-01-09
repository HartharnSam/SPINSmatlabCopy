function q_out = calc_heatflux(params, x, outputs, times)
%CALC_HEATFLUX - Calculates the interfacial "heat" flux across the upper
%boundary of the tank
%
% Inputs:
%    params - output of spins_params
%    x - x grid (output of xgrid_reader)
%    outputs - output numbers to calculate for
%    times - corresponding output times
%
% Other m-files required: 	calc_Nu_turbulent, calc_Ra, cmocean,
% figure_print_format, get_grad2, spins_reader_new
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% Department of Geography & Environmental Sciences, Northumbria University
% email address: sam.hartharn-evans@northumbria.ac.uk
% GitHub: https://github.com/HartharnSam
% 03-Apr-2024; Last revision: 03-Apr-2024
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%% Calculate Heat/Density Fluxes at surface
H = params.Lz; % wall separation
delta = abs(params.rho_top-params.rho_bot);% density difference
q_hat_ts = NaN(size(x, 1), length(outputs));
q_ts = NaN(size(x, 1), length(outputs));

tau = q_hat_ts;
Nu = calc_Nu_turbulent(calc_Ra);
% pre-allocate before parfor to minimise comms
kappa = params.kappa;
Lx = params.Lx;  Nx = params.Nx;
visco = params.visco ; rho_0 = params.rho_0;

warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
pool = gcp;

parfor ii = 1:length(outputs)
    ti = outputs(ii);
    rho = spins_reader_new('rho', ti);
    u = spins_reader_new('u', ti);

    [~, drho_dz] = get_grad2(rho);
    [~, du_dz] = get_grad2(u);

    % Nondimensionalised q (Howland et al 2021 -
    % doi:10.1017/jfm.2021.952 )
    q_hat = -H/delta.*(drho_dz);
    q_hat_ts(:, ii) = q_hat(:, end);

    % Dimensionalised density flux, q (Holland & Jenkins 1999)
    q = -kappa.*Nu.*(drho_dz);
    dx = Lx/Nx;
    q_ts(:, ii) = q(:, end).*dx;

    % Boundary Shear Stress
    tau(:, ii) = visco.*rho_0.*du_dz(:, end);
    %completion(ii, length(outputs), .1, 'Heat Flux loop');
end


if nargout == 0
    figure;
    tiledlayout(3, 1);
    nexttile;
    pcolor(x(:, end), times, q_ts');
    c = colorbar;
    cmocean('thermal');
    ylabel(c, '$q (kgm^{-1}s^{-1})$');
    xlabel('$x (m)$');
    ylabel('$t (s)$');
    axis tight
    nexttile;

    pcolor(x(:, end), times, abs(tau)');
    c = colorbar;
    cmocean('amp');
    axis tight
    ylabel(c, '$\tau$');
    xlabel('$x (m)$');
    ylabel('$t (s)$');
    caxis([0 .2])

    nexttile;
    plot(times, sum(q_ts));
    xlabel('$t (s)$');
    ylabel('$q (kgm^{-1}s^{-1})$');
    axis tight

    save('heat_flux.mat', 'x', 'times', 'q_ts', 'tau');

    figure_print_format(gcf);
    exportgraphics(gcf, 'HeatFlux.png', 'Resolution', 300);
else
    q_out = sum(q_ts)*params.Ly; % Total heat flux timeseries

end
fprintf('\n Heat Fluxes Calculated \n ------------- \n')

end
%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------