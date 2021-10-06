function [tracer_volume weightings] = density_tracer(rho, x_loc)

%% First normalise the density
params = spins_params;
if range(rho(:)) < 10
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
weighted_rho = weighted_rho(x_inds, :);
weightings = weightings(x_inds,:);
tracer_volume = sum(weighted_rho(:));
