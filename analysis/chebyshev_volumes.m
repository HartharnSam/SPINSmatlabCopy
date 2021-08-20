function weights = chebyshev_volumes()
%CHEBYSHEV_VOLUMES - Calculates the volume linked to each grid point in
%chebyshev (z) grid, including over slopes. 
%%
[x, z] = spinsgrid2d;
params = spins_params;
x_weight = params.Lx/params.Nx;
[xx, w] = clencurt(params.Nz-1);
y_weight = w./sum(w) .* max(abs(z), [], 2);
weights = y_weight*x_weight;

%% Test
% calculated_vol = sum(sum(weight));
% actual_volume = params.Lx*params.Lz;
% actual_volume = actual_volume - (.5*params.hill_height*params.hill_height/params.hill_slope);
% if calculated_vol == actual_vol
% disp('Correct');
% end
    end
%%
