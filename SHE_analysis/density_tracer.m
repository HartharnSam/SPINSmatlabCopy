function tracer_volume = density_tracer(rho, x_loc)

%% First normalise the density
params = spins_params;
if range(rho, 'all') < 10
    rho_normalised = (rho+params.delta_rho/2)/params.delta_rho;
else
    rho_normalised = rho;
end
[x, z] = spinsgrid2d;

%% Then calculate weightings
weightings = chebyshev_volumes();
weighted_rho = weightings.*rho_normalised;

%% Find x > x_loc
x_inds = find(x(:, 1) > x_loc);

tracer_volume = sum(weighted_rho(x_inds, :), 'all');
