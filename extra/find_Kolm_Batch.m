function Scales = find_Kolm_Batch(diagnos, gdpar, kappa_min)
% Prints the Kolmogorov and Batchelor Scales and their relation to
% dx & dz
% Kolmogorov Scale: eta = (rho_0 nu^3/epsilon)^(1/4)
% Batchelor Scale:  lambda_B = eta * sqrt(kappa/nu)
% rho_0 is included to make epsilon the dissipation per unit mass (makes dimensions work)
% Adapted from plot_diagnos
% 
% EXAMPLE:
% diagnos = load('diagnostics.mat'); gdpar = 
% params = spins_params(); [gd.x, gd.z] = spinsgrid2d(); gdpar.gd = gd;
% gdpar.params = params;
% kappa_min = find_min_diffu(params); % find minimum diffusivity
% Scales = find_Kolm_Batch(diagnos, gdpar, kappa_min);


% shorten parameters
split_gdpar
rho_0 = params.rho_0;
visco = params.visco;

if isfield(diagnos, 'Max_diss')
    % find max diss. ignoring the first 100 points containing the random perturbations
    if length(diagnos.Diss_tot) >= 100
        max_diss = max(diagnos.Max_diss(100:end));
    else
        max_diss = max(diagnos.Max_diss);
    end
    
    % Kolmogorov and Batchelor scales
    Kolm = (rho_0*visco^3/max_diss)^(1/4);
    Batch = Kolm * sqrt(kappa_min/visco);
    
    % compare max grid size to these scales
    max_dxyz = max(max_grid_spacing(gdpar));
    dx_Kolm  = max_dxyz/Kolm;
    dx_Batch = max_dxyz/Batch;
    
    % print out the info
    fprintf('\n')
    disp('---- Kolmogorov and Batchelor Scales ----')
    disp(['dx/eta =      ',num2str(dx_Kolm)])
    disp(['dx/lambda_B = ',num2str(dx_Batch)])
    
    Scales.Kolm  = Kolm;
    Scales.Batch = Batch;
    Scales.dx_Kolm  = dx_Kolm;
    Scales.dx_Batch = dx_Batch;
else
    Scales = 'N/A';
end
end