function volume = spinsvolume()


%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
params = spins_params;
[x, z] = spinsgrid2d;
params.ndims = size(size(x));

Vol = params.Lx * params.Ly * params.Lz;
if strcmp(params.mapped_grid, 'true')
    if params.ndims == 3
        bot = z(:,1,1);
        top = z(:,1,params.Nz);
    else
        bot = z(:,1);
        top = z(:,params.Nz);
    end
    if strcmp(params.type_x, 'NO_SLIP')
        warning('Volume calculation is not setup for Cheb grid in x.')
    else
        volume = Vol - params.Ly * trapz(x(:, 1), bot+top);
    end
end
%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------