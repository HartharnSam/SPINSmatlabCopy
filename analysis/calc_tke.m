function [k, k_tot] = calc_tke(ii, xinds, zinds)

% Load in raw data
u = spins_reader_new('u', ii, xinds, zinds);
w = spins_reader_new('w', ii, xinds, zinds);

u_var = var(u);
w_var = var(w);
k = 0.5*(u_var + w_var);

if nargout > 1
    params = spins_params; 
    if ~isempty(zinds)
        params.Nz = min(zinds, params.Nz);
    end
    if ~isempty(xinds)
        warning("need to configure variable x for k_tot")
    end

    [~,wci] = clencurt(params.Nz-1);
    wi = wci*(params.Lz)*0.5;
    k_tot = sum(wi.*k).*params.Lx;
end
