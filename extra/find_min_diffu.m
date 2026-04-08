function kappa_min = find_min_diffu(params)
% FIND_MIN_DIFFU
% From plot_diagnos, identifies the minimum diffusivity across all tracers
% 

diffu_types = {'kappa','kappa_rho','kappa_tracer',...
    'kappa_dye','kappa_dye1','kappa_dye2',...
    'kappa_T','kappa_S','kappa_t','kappa_s'};
kappa = nan(1, length(diffu_types));
for ii = 1:length(diffu_types)
    if isfield(params, diffu_types{ii})
        kappa(ii) = params.(diffu_types{ii});
    else
        kappa(ii) = NaN;
    end
end
kappa_min = min(kappa(:));
end