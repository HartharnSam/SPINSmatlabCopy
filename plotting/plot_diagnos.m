function all_diagnos = plot_diagnos(make_plots, do_filter, save_plots)  
%  PLOT_DIAGNOS  Plot diagnostics (timing, energy budget, etc) of a SPINS case
%
%  Info: Diagnostics are in one of two diagnostic files
%           1) diagnostics.txt - raw data
%           2) diagnostics.mat - cleaned data
%                  run clean_diagnostics to remove time overlaps
%
%        The Energy diagnostics are only true for a closed system with no external forcing
%        ie. no tidal forcing, or forced shear flow
%        Those energy transfer terms must be computed on their own
%
%  Usage:
%   diagnos = plot_diagnos();
%
%  Inputs:
%    Optional argument:
%   make_plots      - {boolean} Do you want to make the plots? (default: true)
%   do_filter       - {boolean} Do you want to filter BPE_tot? (default: false)
%                               This smooths the energy budget plots
%
%  Outputs:
%   'all_diagnos'   - a structure containing multiple diagnostic structures:
%                       * diagnos       - diagnostics directly from the diagnostic file
%                       * EnergyBudget  - Energy in each form (KE, APE, etc) and amounts transfered between them
%                       * EnergyRates   - Rates of change of different energy components (KE, APE, etc)
%                                         and the rates of transfer between them
%                       * Mixing        - Mixing efficiencies
%                       * Scales        - Kolmogorov and Batchelor length scales
%
%  David Deepwell, 2019

%% set default plotting option
if nargin == 0
    make_plots = true;
    do_filter = false;
    save_plots = true;  
elseif nargin == 1
    do_filter = false;
    save_plots = true; 
elseif nargin == 2
    save_plots = true;  
end
%% 
%%%%%%%%%%% Read Data %%%%%%%%%%%%%%%
% read diagnostic file
diag_file_name = 'diagnostics';
if exist([diag_file_name,'.mat'], 'file') == 2
    % load cleaned data, if it exists
    diag_file = [diag_file_name,'.mat'];
    diagnos = load(diag_file);
    % also load the txt file to check if mat file is up to date
    diag_txt = readtable([diag_file_name,'.txt']);
elseif exist([diag_file_name,'.txt'], 'file') == 2
    % otherwise read raw data (this will likely be fine, unless a restart caused overlapping data)
    diag_file = [diag_file_name,'.txt'];
    diagnos = table2struct(readtable(diag_file), 'ToScalar', true);
else
    error('Diagnostics file not found.')
end

% load parameters and grid (read separately so that it's not a global variable)
params = spins_params();
if strcmp(params.mapped_grid, 'true')
    %gd.x = xgrid_reader;
    
    [gd.x gd.z] = spinsgrid2d();
else
    gd = spins_grid('vector');
end
params.ndims = length(fieldnames(gd)); % Number of dimensions
gdpar.gd = gd;
gdpar.params = params;

% place into output structure
all_diagnos.diagnos = diagnos;

%%%%%%%%%%% Parse Data %%%%%%%%%%%%%%%
% shorten parameters
rho_0 = params.rho_0;
visco = params.visco;
Vol = find_volume(gdpar); % Volume of domain
kappa_min = find_min_diffu(params); % find minimum diffusivity
% get time variables
try
    clk_time = diagnos.Clock_time;
    sim_time = diagnos.Time;
    if exist('diag_txt', 'var')
        sim_txt_time = diag_txt.Time; % simulation time from raw diagnostics.txt file
    end
catch
    error(['Diagnostics file incorrectly labeled the time (clock time and simulation time).',...
        newline,'Should be "Clock_time" and "Time"'])
end
if sum(diff(sim_time) < 0) > 0
    error(['Simulation has been restarted and there are overlapping times,',...
        newline,'Run "clean_diagnostics" to remove the overlaps.'])
end

%% Parse time into step times and parse restarts
% find wall clock and simulation step time
clk_step_time = [clk_time(1); clk_time(2:end)-clk_time(1:end-1)];
sim_step_time = [sim_time(1); sim_time(2:end)-sim_time(1:end-1)];
clk_step_time(clk_step_time < 0) = 0;
% check for restarts
restart_ind = find(diagnos.Iter == 1);
restart_ind(1) = 1;
strt_ind = [restart_ind; length(diagnos.Iter)+1];
N_start = length(restart_ind);
n_steps = length(diagnos.Iter);
% total steps and clock time
tot_clk_time = sum(clk_time(strt_ind(2:end)-1));
% clock time minus write time
if exist('plot_times.txt', 'file')
    % if outputs exist, remove output write time
    warning('off','all'); % suppress warning for this only
    plottimes = readtable('plot_times.txt');
    warning('on','all');
    avg_write = nanmean(plottimes.WriteTime_s_);
    tot_write =  nansum(plottimes.WriteTime_s_);
else
    avg_write = 0;
    tot_write = 0;
end
avg_clk_step = (tot_clk_time - tot_write)/n_steps;
avg_sim_step = sim_time(end)/n_steps;
avg_clk_per_sim = avg_clk_step/avg_sim_step;
clk_per_sim = clk_step_time./sim_step_time;
% estimated clock time required to complete simulation
if ~isfield(params, 'final_time')
    params.final_time = diag_txt{end, 2};
end
sim_time_remain = params.final_time - sim_time(end);
if length(clk_per_sim) > 50
    clk_time_remain = sim_time_remain*mean(clk_per_sim(end-50:end));
else
    clk_time_remain = sim_time_remain*mean(clk_per_sim(2:end));
end

% legend labels for each run
run_label = cell(1,N_start);
for ii = 1:N_start
    run_label{ii} = ['Run ',num2str(ii)];
end

%%%%%%%%%%% Print Timing info %%%%%%%%%%%%%%%
% functions for converting seconds into other units
s2ms   = @(sec) datestr(datenum(0,0,0,0,0,sec),'MM:SS');
s2hms  = @(sec) datestr(datenum(0,0,0,0,0,sec),'HH:MM:SS');%#ok
s2dhm  = @(sec) datestr(datenum(0,0,0,0,0,sec),'DD:HH:MM');
s2dhms = @(sec) datestr(datenum(0,0,0,0,0,sec),'DD:HH:MM:SS');
% print info
fprintf('\n')
disp('---- Timings ----')
if exist('diag_txt', 'var')
    if sim_txt_time(end) - sim_time(end) > 0
        fprintf(2,'Sim. is %5.2g s ahead of mat file, re-run "clean_diagnostics"\n\n',sim_txt_time(end) - sim_time(end))
    end
end
fprintf('Most recent sim. time: %5.2f s\n',sim_time(end))
fprintf('Total clock time:      %s (D:H:M)\n',s2dhm(tot_clk_time))
%% print out average write time and step time
if avg_write > 0
    if avg_write < 10
        fprintf('Average write time: %8.3f s\n',avg_write)
    else
        fprintf('Average write time:    %s (M:S)\n',s2ms(avg_write))
    end
end
if avg_sim_step < 1e-2
    fprintf('Avg. sim. step time: %7.3f ms\n',avg_sim_step*1000)
else
    fprintf('Avg. sim. step time: %7.3f s\n',avg_sim_step)
end
if avg_clk_step < 10
    fprintf('Avg. clock step time:  %5.3f s\n',avg_clk_step)
else
    fprintf('Avg. clock step time:  %s (M:S)\n',s2ms(avg_clk_step))
end
fprintf('Avg. clock time per sim. sec.: %s (D:H:M:S)\n',s2dhms(avg_clk_per_sim))
if clk_time_remain < 60*60*24*32
    fprintf('Est. clock time remaining:     %s (D:H:M:S)\n',s2dhms(clk_time_remain))
else
    fprintf('Est. clock time remaining:     01:%s (M:D:H:M:S)\n',s2dhms(clk_time_remain-60*60*24*32))
end
all_diagnos.diagnos.AvgClockTimePerSimSec = (avg_clk_per_sim);
all_diagnos.diagnos.TotClockTime = (tot_clk_time);



%%%%%%%%%%% Which optional diagnostics were computed? %%%%%%%%%%%%%%%
if isfield(diagnos, 'Enst_y_tot')
    compute_enstrophy = true;
else
    compute_enstrophy = false; 
end
if isfield(diagnos, 'Diss_tot') 
    compute_dissipation = true;
else
    compute_dissipation = false; 
end
if isfield(diagnos, 'BPE_tot') 
    compute_BPE = true;
else
    compute_BPE = false; 
end
if isfield(diagnos, 'BPE_from_int') 
    compute_BPE_from_int = true;
else
    compute_BPE_from_int = false; 
end

%%%%%%%%%%% Shorten entstrophy variable and compare against dissipation %%%%%%%%%%%%%%%
% from Mike Waite's Turbulence class: 
%   total diss. = 2*mu*(total enstrophy)
%   when boundary conditions are not no-slip (ie. free-slip or periodic on all sides)
%   see also Yeung et al. 2012 in JFM
% check how close the ratio is to 1
if compute_enstrophy
    if params.ndims == 3
        enst_x_tot = diagnos.Enst_x_tot;
        enst_y_tot = diagnos.Enst_y_tot;
        enst_z_tot = diagnos.Enst_z_tot;
    else
        enst_x_tot = 0;
        enst_y_tot = diagnos.Enst_y_tot;
        enst_z_tot = 0;
    end
    enst_tot = enst_x_tot + enst_y_tot + enst_z_tot;
    if compute_dissipation
        enst_diss = diagnos.Diss_tot./(enst_tot*2*visco*rho_0);
        all_diagnos.diagnos.enst_diss = enst_diss; % place into output structure
    end
    all_diagnos.diagnos.Enst_tot  = enst_tot; % place into output structure
end

%%%%%%%%%%% Print Kolmogorov and Batchelor scales %%%%%%%%%%%%%%%
Scales = find_Kolm_Batch(diagnos, gdpar, kappa_min);
all_diagnos.Scales = Scales; % place into output structure

%%%%%%%%%%% Max density variation %%%%%%%%%%%%%%%
if isfield(diagnos, 'Max_density')
    rho_init = diagnos.Max_density(1);
    rho_var = diagnos.Max_density/rho_init - 1;
    rho_max_diff = diagnos.Max_density - rho_init;
    fprintf('\n')
    disp('---- Max density deviation ----')
    fprintf('Max {Max rho / Max rho(0)} - 1 =  %4.2g\n',max(rho_var));
    fprintf('Max {Max rho - Max rho(0)}     =  %4.2g\n',max(rho_max_diff));
    % place into output structure
    all_diagnos.diagnos.rho_var = rho_var;
end
if isfield(diagnos, 'Min_density')
    rho_min_init = diagnos.Min_density(1);
    rho_min_var = diagnos.Min_density/rho_min_init - 1;
    rho_min_diff = diagnos.Min_density - rho_min_init;
    fprintf('Min {Min rho / Min rho(0)} - 1 =  %4.2g\n',min(rho_var));
    fprintf('Min {Min rho - Min rho(0)}     =  %4.2g\n',min(rho_min_diff));
    % place into output structure
    all_diagnos.diagnos.rho_min_var = rho_min_var;
end
if isfield(diagnos, 'Max_temperature') || isfield(diagnos, 'Max_salinity') || ...
        isfield(diagnos, 'Min_temperature') || isfield(diagnos, 'Min_salinity')
    fprintf('\n')
    disp('---- Max and Min salt/temp change ----')
    if isfield(diagnos, 'Max_temperature')
        temp_max_init = diagnos.Max_temperature(1);
        temp_max_diff = diagnos.Max_temperature - temp_max_init;
        fprintf('Max {Max T - Max T(0)} =  %4.2f\n',max(temp_max_diff));
    end
    if isfield(diagnos, 'Min_temperature')
        temp_min_init = diagnos.Min_temperature(1);
        temp_min_diff = diagnos.Min_temperature - temp_min_init;
        fprintf('Max {Min T(0) - Min T} =  %4.2f\n',max(-temp_min_diff));
    end
    if isfield(diagnos, 'Max_salinity')
        salt_max_init = diagnos.Max_salinity(1);
        salt_max_diff = diagnos.Max_salinity - salt_max_init;
        fprintf('Max {Max s - Max s(0)} =  %4.2f\n',max(salt_max_diff));
    end
    if isfield(diagnos, 'Min_salinity')
        salt_min_init = diagnos.Min_salinity(1);
        salt_min_diff = diagnos.Min_salinity - salt_min_init;
        fprintf('Max {Min s(0) - Min s} =  %4.2f\n',max(-salt_min_diff));
    end
end

%%%%%%%%%%% Parse and Print Energy diagnostics %%%%%%%%%%%%%%%
% compute total KE
if ~isfield(diagnos, 'KE_y')
    diagnos.KE_y = diagnos.KE_x*0;
end
KE_tot = diagnos.KE_x + diagnos.KE_y + diagnos.KE_z;

% is there KE forcing?
if isfield(diagnos, 'KE_from_forcing')
    KE_forcing = true;
    F2KE_rate = diagnos.KE_from_forcing;
    F2KE_tot  = cumtrapz(diagnos.Time, F2KE_rate);
else
    KE_forcing = false;
    F2KE_rate = 0;
    F2KE_tot = 0;
end

% total energy
E_tot  = KE_tot + diagnos.PE_tot;
E_loss = E_tot(1) - E_tot;
% compute rates
time_rate = linspace(sim_time(1), sim_time(end),round(length(sim_time)/2));
Dmat = FiniteDiff(time_rate, 1, 2, true);
KE_rate = Dmat*interp1(sim_time, KE_tot, time_rate, 'pchip')';
E_rate  = Dmat*interp1(sim_time, E_tot,  time_rate, 'pchip')';

% compute APE, total Avail. Energy (AE_tot), and change in BPE, and their rates of change
if compute_BPE
    if do_filter
        BPE_tot = mlptdenoise(diagnos.BPE_tot, diagnos.Time, 5, 'DualMoments', 3);
        
    else
        BPE_tot = diagnos.BPE_tot;
    end
    APE_tot = diagnos.PE_tot - BPE_tot;
    AE_tot  = E_tot - BPE_tot;
    E0_and_W = AE_tot(1) + F2KE_tot(end);
    AE_loss = AE_tot(1) - AE_tot;
    BPE_change = BPE_tot - BPE_tot(1);
    BPE_rate = Dmat*interp1(sim_time, BPE_change, time_rate, 'pchip')';
    APE_rate = Dmat*interp1(sim_time, APE_tot,    time_rate, 'pchip')';
    AE_rate  = Dmat*interp1(sim_time,  AE_tot,    time_rate, 'pchip')';
end
% energy transfered from internal to BPE (phi_i)
if compute_BPE_from_int
    Int2BPE_rate = diagnos.BPE_from_int;
    Int2BPE_tot  = cumtrapz(diagnos.Time, Int2BPE_rate);
    % energy transfered from APE to BPE (phi_m = phi_d - phi_i)
    if compute_BPE
        APE2BPE_tot  = BPE_change - Int2BPE_tot;
        APE2BPE_rate = Dmat*interp1(sim_time, APE2BPE_tot, time_rate, 'pchip')';
    end
end
% energy dissipation rates and cumulative amounts
if compute_dissipation
    % dissipation (epsilon)
    KE2Int_rate = diagnos.Diss_tot;
    KE2Int_tot  = cumtrapz(diagnos.Time, KE2Int_rate);
    % KE to APE (phi_z)
    KE2APE_tot  = -(KE_tot - KE_tot(1)) - KE2Int_tot + F2KE_tot;
    KE2APE_rate = Dmat*interp1(sim_time, KE2APE_tot, time_rate, 'pchip')';
    % offset onto the grid to use with other rates which are off_set by 1st order FD derivative
    diss_offset = interp1(sim_time, diagnos.Diss_tot, time_rate, 'pchip')';
    % compute rate of change of internal energy, and the cumulative change 
    if compute_BPE_from_int
        Int_rate = diagnos.Diss_tot - diagnos.BPE_from_int;
        Int_change = cumtrapz(diagnos.Time, Int_rate);
    end
end

% compute energy created/removed by numerics (filter and numerical errors/scheme etc.)
% only if dissipation, BPE, and BPE_from_int were calculated
if compute_dissipation && compute_BPE && compute_BPE_from_int
    NumE_tot  = E_loss - Int_change + F2KE_tot;
    NumE_rate = Dmat*interp1(sim_time, NumE_tot, time_rate, 'pchip')';

    fprintf('\n')
    fprintf('---- Energy ----\n')
    if KE_forcing
        fprintf('As a percentage of the initial available energy + work done\n')
    else
        fprintf('As a percentage of the initial available energy\n')
    end
    if KE_forcing
        fprintf('  and total work applied by forcing\n')
    end
    fprintf('Available energy lost:   %6.2f %%\n',AE_loss(end)/E0_and_W*100)
    fprintf('  > through mixing:        %6.2f %%\n',APE2BPE_tot(end)/E0_and_W*100)
    fprintf('  > through dissipation:   %6.2f %%\n',KE2Int_tot(end)/E0_and_W*100)
    fprintf('  > through numerics:      %6.2f %%\n',NumE_tot(end)/E0_and_W*100)
    if KE_forcing
        fprintf('  > through forcing:       %6.2f %%\n',-F2KE_tot(end)/E0_and_W*100)
    end
    fprintf('Total Mech. energy lost: %6.2f %%\n',E_loss(end)/E0_and_W*100)
    fprintf('  > through diffusion:     %6.2f %%\n',-Int2BPE_tot(end)/E0_and_W*100)
    fprintf('  > through dissipation:   %6.2f %%\n',KE2Int_tot(end)/E0_and_W*100)
    fprintf('  > through numerics:      %6.2f %%\n',NumE_tot(end)/E0_and_W*100)
    fprintf('Change in BPE:           %6.2f %%\n',BPE_change(end)/E0_and_W*100)
end

% place into output structure
%  Energy budget
all_diagnos.EnergyBudget.Time       = diagnos.Time;
all_diagnos.EnergyBudget.E_tot      =   E_tot;
all_diagnos.EnergyBudget.KE_tot     =  KE_tot;
if compute_BPE
    all_diagnos.EnergyBudget.AE_tot     =  AE_tot;
    all_diagnos.EnergyBudget.APE_tot    = APE_tot;
    all_diagnos.EnergyBudget.BPE_change = BPE_change;
end
if compute_dissipation && compute_BPE_from_int
    all_diagnos.EnergyBudget.Int_change = Int_change; end
if compute_dissipation && compute_BPE && compute_BPE_from_int
    all_diagnos.EnergyBudget.NumE_tot   = NumE_tot; end
if KE_forcing
    all_diagnos.EnergyBudget.F2KE_tot = F2KE_tot; end
% Energy conversion budget
if compute_dissipation
    all_diagnos.EnergyBudget.KE2APE_tot  = KE2APE_tot;
    all_diagnos.EnergyBudget.KE2Int_tot  = KE2Int_tot;
end
if compute_BPE
    all_diagnos.EnergyBudget.APE2BPE_tot = APE2BPE_tot; end
if compute_BPE_from_int
    all_diagnos.EnergyBudget.Int2BPE_tot = Int2BPE_tot; end
% Energy rates
all_diagnos.EnergyRates.Time      = time_rate';
all_diagnos.EnergyRates.E_rate    =   E_rate;
all_diagnos.EnergyRates.KE_rate   =  KE_rate;
if compute_BPE
    all_diagnos.EnergyRates.AE_rate   =  AE_rate;
    all_diagnos.EnergyRates.APE_rate  = APE_rate;
    all_diagnos.EnergyRates.BPE_rate  = BPE_rate;
end
if compute_dissipation && compute_BPE && compute_BPE_from_int
    all_diagnos.EnergyRates.NumE_rate = NumE_rate;
end
if compute_dissipation && compute_BPE_from_int
    all_diagnos.EnergyRates.Int_rate  = interp1(sim_time, Int_rate,  time_rate, 'pchip')'; % note, should use original data
end
if KE_forcing
    all_diagnos.EnergyRates.F_rate = interp1(sim_time, F2KE_rate, time_rate, 'pchip')'; % note, should use original data
end
% Energy conversion rates
if compute_dissipation
    all_diagnos.EnergyRates.KE2APE_rate  = KE2APE_rate;
    all_diagnos.EnergyRates.KE2Int_rate  = interp1(sim_time, diagnos.Diss_tot,     time_rate, 'pchip')'; % note, should use original data
end
if compute_BPE
    all_diagnos.EnergyRates.APE2BPE_rate = APE2BPE_rate;
end
if compute_BPE_from_int
    all_diagnos.EnergyRates.Int2BPE_rate = interp1(sim_time, diagnos.BPE_from_int, time_rate, 'pchip')'; % note, should use original data
end

% compute mixing efficiencies (see Gregg 2018)
if compute_BPE
    fprintf('\n')
    fprintf('---- Mixing efficiency ----\n')
    fprintf('Change in BPE / |change in AE|: %6.2f \n', ...
        BPE_change(end)/abs(AE_tot(end) - AE_tot(1)))
    all_diagnos.Mixing.BPE_change_perc = BPE_change(end)/abs(AE_tot(end) - AE_tot(1)); % place into output structure
    if compute_dissipation
        % ignore first second in mixing efficiency numbers, this time is dominated by dissipation of the initial random perturbations
        inds = nearest_index(time_rate,1):length(BPE_rate);
        % We include the energy removed by the numerics as this is likely
        % the filter acting as numerical viscosity
        mix_cof = APE2BPE_rate./(diss_offset + NumE_rate);
        mix_eff = APE2BPE_rate./(diss_offset + NumE_rate + APE2BPE_rate);
        mix_eff_cum = APE2BPE_tot./(KE2Int_tot + APE2BPE_tot + NumE_tot);
        fprintf('Cumulative mixing eff.:         %6.2f \n', mix_eff_cum(end))
        fprintf('Max inst. mixing eff.:          %6.2f \n', max(mix_eff(inds)))
        fprintf('Max inst. mixing coeff.:        %6.2f \n', max(mix_cof(inds)))
        % place into output structure
        all_diagnos.Mixing.Time          = time_rate';
        warning('off','all'); % suppress warning for this only
        all_diagnos.Mixing.mix_eff_cum   = interp1(sim_time, mix_eff_cum, time_rate, 'pchip')'; % note, should use original data
        warning('on','all'); % suppress warning for this only
        all_diagnos.Mixing.mix_eff       = mix_eff;
        all_diagnos.Mixing.mix_coeff     = mix_cof;
        all_diagnos.Mixing.mix_eff_max   = max(mix_eff);
        all_diagnos.Mixing.mix_coeff_max = max(mix_cof);
    end
end

%% %%%%%%%%%%% make plots %%%%%%%%%%%%%%%
if make_plots
    fn = 20;    % first figure number
    fm = fn+15; % first figure number for unknown diagnostics
    cols = get(groot,'DefaultAxesColorOrder');
    coly = cols(2,:);
    tr_ind = 0; % tracer number

    for name = fieldnames(diagnos)'
        if strcmp(name, 'Iter')
            continue
            
            %%%% Clock Time per Step %%%%
        %% %% Diag_SimRunRate %% %%
        elseif strcmp(name, 'Clock_time')
            figure(fn), clf
            subplot(3,1,1), cla reset
            hold on
            for ii = 1:N_start
                if ii == 1
                    its = (strt_ind(ii)+1):strt_ind(ii+1)-1;
                else
                    its = strt_ind(ii):strt_ind(ii+1)-1;
                end
                plot(diagnos.Iter(its) + strt_ind(ii)-1, clk_step_time(its), '.')
            end
            ylabel('\Delta T_c (s)')
            title('Clock time step')
            if max(clk_step_time)/mean(clk_step_time) > 10
                yl = ylim();
                set(gca, 'YLim', [yl(1) prctile(clk_step_time,99.98)]);
            end
            legend(run_label)
            legend('location','best')
            %legend('boxoff')
            box on
            hold off

            %%%% Simulation Time per Step %%%%
        elseif strcmp(name, 'Time')
            figure(fn)
            subplot(3,1,2), cla reset
            hold on
            for ii = 1:N_start
                its = strt_ind(ii):strt_ind(ii+1)-1;
                plot(diagnos.Iter(its) + strt_ind(ii)-1, sim_step_time(its), '.')
            end
            ylabel('\Delta T_s (s)')
            title('Simulation time step')
            if length(sim_step_time) > 500
                set(gca, 'YLim', [0 1.2*max(sim_step_time(500:end))]);
            end
            box on
            hold off
            
            
            %%%% Clock Time per Simulation Time per Step %%%%
            subplot(3,1,3), cla reset
            hold on
            for ii = 1:N_start
                its = strt_ind(ii):strt_ind(ii+1)-1;
                plot(diagnos.Iter(its) + strt_ind(ii)-1, clk_per_sim(its), '.')
            end
            xlabel('Iteration')
            ylabel('\Delta{T_c} / \Delta{T_s}')
            title('Ratio of time steps')
            try
                max_clk_per_sim = max(clk_per_sim(~isoutlier(clk_per_sim,'ThresholdFactor',15)));
            catch
                max_clk_per_sim = max(clk_per_sim);
            end
            set(gca, 'YLim', [0 1.1*max_clk_per_sim]);
            box on
            hold off
    
            
            %% %% Maximum (absolute value of) Velocities %% %%
        elseif strcmp(name, 'Max_u')
            figure(fn+1)
            subplot(2,1,1), cla reset
            if ~isfield(diagnos, 'Max_v')
                diagnos.Max_v = diagnos.Max_u*0;
            end
            plot(diagnos.Time, [diagnos.Max_u, diagnos.Max_v, diagnos.Max_w, diagnos.Max_vel])
            if params.ndims == 3
                set(gca,'yscale','log');
            else
                set(gca,'yscale','linear');
            end
            ylabel('Velocity (m/s)')
            title('Max velocity')
            leg = legend({'Max $|u|$','Max $|v|$','Max $|w|$','Max $|\vec{u}|$'},...
                'Interpreter','Latex');
            leg.Location = 'best';
            leg.Box = 'off';
        elseif strcmp(name, 'Max_v') || strcmp(name, 'Max_w') || strcmp(name, 'Max_vel')
            continue

            %%%% Energy components %%%%
        elseif strcmp(name, 'KE_x')
            %% Plot just kinetic energy
            figure(fn+1)
            subplot(2,1,2), cla reset
            plot(diagnos.Time,[diagnos.KE_x, diagnos.KE_y, diagnos.KE_z, KE_tot])
            if params.ndims == 3
                set(gca,'yscale','log');
            else
                set(gca,'yscale','linear');
            end
            xlabel('time (s)')
            ylabel('Kinetic Energy (J)')
            title('KE components')
            leg = legend({'KE_x','KE_y','KE_z','KE_{tot}'});
            leg.Location = 'best';
            leg.Box = 'off';
            
            %% Plot all energy components
            %(KE, APE, Change in internal, Change in BPE, and numerical)
            figure(fn+2)
            subplot(2,1,1), cla reset
            hold on
            plot(diagnos.Time, KE_tot)
            energy_label = {'KE'};
            % Add APE, AE, and change in BPE, if calculated
            if compute_BPE
                plot(diagnos.Time, [APE_tot, AE_tot, BPE_change])
                energy_label = [energy_label, {'APE','AE','\Delta{BPE}'}]; %#ok
            end
            % Add change in internal energy if calculated
            if compute_dissipation && compute_BPE_from_int
                plot(diagnos.Time, Int_change,'k')
                energy_label = [energy_label, '\Delta{Int}'];%#ok
            end
            % add energy created/removed by numerics (filter and numerical errors/scheme etc.)
            % only if dissipation, BPE, and BPE_from_int were calculated
            if compute_dissipation && compute_BPE && compute_BPE_from_int
                plot(diagnos.Time, NumE_tot)
                %plot([0 diagnos.Time(end)],[1 1]*AE_tot(1),'k')
                %energy_label = [energy_label, {'Numerics','BPE from int','AE(0)'}];
                energy_label = [energy_label, {'Numerics'}]; %#ok
            end
            if KE_forcing
                plot(diagnos.Time, -F2KE_tot)
                energy_label = [energy_label, {'- Work'}]; %#ok
            end
            grid on
            box on
            ylabel('Energy (J)')
            title('Energy components')
            legend(energy_label)
            legend('location','best')
            legend('boxoff')

            %% plot the rate of change of Energy components
            subplot(2,1,2), cla reset
            hold on
            inds = 10:length(time_rate);
            plot(time_rate(inds), KE_rate(inds))
            energy_label = {'KE'};
            % Add APE, AE, and change in BPE, if calculated
            if compute_BPE
                plot(time_rate(inds), [APE_rate(inds), AE_rate(inds), BPE_rate(inds)])
                energy_label = [energy_label, {'APE','AE','BPE (\phi_d)'}]; %#ok
            end
            % Add change in internal energy if calculated
            if compute_dissipation && compute_BPE_from_int
                plot(diagnos.Time(10:end), Int_rate(10:end),'k')
                energy_label = [energy_label, 'Internal']; %#ok
            end
            % add energy created/removed by numerics (filter and numerical errors/scheme etc.)
            % only if dissipation, BPE, and BPE_from_int were calculated
            if compute_dissipation && compute_BPE && compute_BPE_from_int
                plot(time_rate(inds), NumE_rate(inds))
                energy_label = [energy_label, {'Numerics'}];%#ok
            end
            if KE_forcing
                plot(diagnos.Time, -F2KE_rate)
                energy_label = [energy_label, {'- Work'}]; %#ok
            end
            grid on
            box on
            xlabel('time (s)')
            ylabel('Rate (J/s)')
            title('Energy rates of change')
            legend(energy_label)
            legend('location','best')
            legend('boxoff')
            
            
            %% Plot all energy component rates (KE, APE, Diss, Change in BPE)
            figure(fn+3)
            subplot(2,1,1), cla reset
            hold on
            % interpolate onto a regular grid of similar size
            energy_label = {};
            if compute_dissipation
                plot(diagnos.Time, KE2APE_tot)
                energy_label = {'KE to APE (\int \phi_z dt)'};
            end
            % Add APE if calculated
            if compute_BPE && compute_BPE_from_int
                plot(diagnos.Time, APE2BPE_tot)
                energy_label = [energy_label, {'APE to BPE (\int \phi_m dt)'}];%#ok
            end
            % Add energy converted from internal energy
            if compute_BPE_from_int
                plot(diagnos.Time, Int2BPE_tot,'Color',cols(4,:))
                energy_label = [energy_label, 'Int to BPE (\int \phi_i dt)'];%#ok
            end
            % Add dissipation if calculated
            if compute_dissipation
                plot(diagnos.Time(8:end), KE2Int_tot(8:end),'k')
                energy_label = [energy_label, 'KE to Int (diss, \int \epsilon dt)'];%#ok
            end
            if compute_dissipation && compute_BPE && compute_BPE_from_int
                plot(diagnos.Time, NumE_tot,'Color',cols(5,:))
                energy_label = [energy_label, 'to Numerics'];%#ok
            end
            if KE_forcing
                plot(diagnos.Time, F2KE_tot)
                energy_label = [energy_label, {'Work to KE'}];%#ok
            end
            grid on
            ylabel('Energy (J)')
            title('Energy converted')
            legend(energy_label)
            legend('location','best')
            legend('boxoff')
            box on
            hold off

            subplot(2,1,2), cla reset
            hold on
            inds = 5:length(time_rate);
            % interpolate onto a regular grid of similar size
            energy_label = {};
            if compute_dissipation
                plot(time_rate(inds), KE2APE_rate(inds))
                energy_label = {'KE  to APE (\phi_z)'};
            end
            % Add APE if calculated
            if compute_BPE && compute_BPE_from_int
                plot(time_rate(inds), APE2BPE_rate(inds))
                energy_label = [energy_label, {'APE to BPE (\phi_m)'}];%#ok
            end
            % Add energy converted from internal energy
            if compute_BPE_from_int
                plot(diagnos.Time, diagnos.BPE_from_int,'Color',cols(4,:))
                energy_label = [energy_label, 'Int to BPE (\phi_i)']; %#ok
            end
            % Add dissipation if calculated
            if compute_dissipation
                plot(diagnos.Time(15:end), diagnos.Diss_tot(15:end),'k')
                energy_label = [energy_label, 'KE  to Int (diss, \epsilon)'];%#ok
            end
            if compute_dissipation && compute_BPE && compute_BPE_from_int
                plot(time_rate(inds), NumE_rate(inds),'Color',cols(5,:))
                energy_label = [energy_label, 'to Numerics'];%#ok
            end
            if KE_forcing
                plot(diagnos.Time, F2KE_rate)
                energy_label = [energy_label, {'Work to KE'}];%#ok
            end
            grid on
            xlabel('time (s)')
            ylabel('Rate (J/s)')
            title('Energy conversion rates')
            legend(energy_label)
            legend('location','best')
            legend('boxoff')
            box on
            hold off
            
        elseif strcmp(name, 'KE_y') || strcmp(name, 'KE_z') || strcmp(name, 'BPE_tot') ||...
                strcmp(name, 'PE_tot') || strcmp(name, 'BPE_from_int') || strcmp(name, 'KE_from_forcing')
            continue

            %%%% Total Mass %%%%
        elseif strcmp(name, 'Mass')
            figure(fn+4)
            subplot(3,1,1), cla reset
            first_ind = 50;
            % ignore the early adjustment caused by
            % the non-divergence free initial random perturbation
            if length(diagnos.Mass) > first_ind
                M0 = mean(diagnos.Mass(first_ind-20:first_ind));
                plot(diagnos.Time(first_ind-20:end),...
                    (diagnos.Mass(first_ind-20:end) - M0)/M0,'.')
                ylabel('M/M_0 - 1')
                title('Mass deviation')
            end

            %%%% Maximum Density or Salt/Temp %%%%
        elseif strcmp(name, 'Max_density')
            figure(fn+4)
            subplot(3,1,2), cla reset
            if isfield(diagnos,'Min_density')
                plot(diagnos.Time, [rho_max_diff -rho_min_diff]);
                ylabel('Density deviation')
                title('Max/Min density deviation')
                legend({'$\rho_{max} - \rho_{max}(0)$','$\rho_{min}(0) - \rho_{min}$'},...
                    'Interpreter','Latex','FontSize',12)
                legend('location','best')
                legend('boxoff')
                ax = gca;
                ax.YGrid = 'on';
            else
                plot(diagnos.Time, rho_var)
                ylabel('$\rho_\mathrm{max}/\rho_\mathrm{max}(0) - 1$',...
                    'Interpreter','Latex','FontSize',12)
                title('Max density deviation')
            end
        elseif strcmp(name, 'Min_density')
            continue
        elseif strcmp(name, 'Max_temperature') || strcmp(name, 'Max_salinity') || ...
                strcmp(name, 'Min_temperature') || strcmp(name, 'Min_salinity')
            figure(fn+10)
            if strcmp(name, 'Max_temperature') || strcmp(name, 'Min_temperature')
                subplot(2,1,1), cla reset, hold on
                leg_text = {};
                if strcmp(name, 'Max_temperature')
                    plot(diagnos.Time, temp_max_diff)
                    leg_text = [leg_text,'$T_{max} - T_{max}(0)$'];%#ok
                end
                if strcmp(name, 'Min_temperature')
                    plot(diagnos.Time, -temp_min_diff)
                    leg_text = [leg_text,'$T_{min}(0) - T_{min}$'];%#ok
                end
                ylabel('Temperature')
                legend(leg_text, 'Interpreter','Latex','FontSize',12)
                legend('location','best')
                legend('boxoff')
                ax = gca;
                ax.YGrid = 'on';
            end
            if strcmp(name, 'Max_salinity') || strcmp(name, 'Min_salinity')
                subplot(2,1,2), cla reset, hold on
                leg_text = {};
                if strcmp(name, 'Max_salinity')
                    plot(diagnos.Time, salt_max_diff)
                    leg_text = [leg_text,'$s_{max} - s_{max}(0)$'];%#ok
                end
                if strcmp(name, 'Min_salinity')
                    plot(diagnos.Time, -salt_min_diff)
                    leg_text = [leg_text,'$s_{min}(0) - s_{min}$'];%#ok
                end
                ylabel('Salinity')
                xlabel('t (s)')
                legend(leg_text, 'Interpreter','Latex','FontSize',12)
                legend('location','best')
                legend('boxoff')
                ax = gca;
                ax.YGrid = 'on';
            end


            %%%% Max Dissipation %%%%
        elseif strcmp(name, 'Max_diss')
            first_ind = 20;
            if length(diagnos.Max_diss) > first_ind
                figure(fn+4)
                subplot(3,1,3), cla reset
                inds = first_ind:length(diagnos.Max_diss);
                plot(diagnos.Time(inds), diagnos.Max_diss(inds))
                xlabel('time (s)')
                ylabel('$\epsilon_\mathrm{max}$  (J s$^{-1}$ m$^{-3}$)',...
                    'Interpreter','Latex','FontSize',12)
                title('Max dissipation')
            end

            %%%% Area where density is over/under initial extents %%%%
        elseif strcmp(name, 'Rho_over_vol')
            figure(fn+5), clf
            plot(diagnos.Time, [diagnos.Rho_over_vol diagnos.Rho_under_vol ...
                diagnos.Rho_over_extra_vol diagnos.Rho_under_extra_vol])
            xlabel('time (s)')
            ylabel('Volume where $\xi$ is true','Interpreter','Latex')
            %ylabel('$1/E_\mathrm{tot}\int_\Omega E(\xi) dV$','Interpreter','Latex')
            legend({'\xi = (\rho > \rho_{max})',...
                '\xi = (\rho < \rho_{min})',...
                '\xi = (\rho > 1.00001 \rho_{max})',...
                '\xi = (\rho < 1.00001 \rho_{min})'})
            legend('location','best')
            legend('boxoff')
        elseif strcmp(name, 'Rho_under_vol') || ...
                strcmp(name, 'Rho_over_extra_vol') || strcmp(name, 'Rho_under_extra_vol')
            continue

            %%%% Maximum of Dye or Tracer %%%%
        elseif strncmp(name, 'Max_dye',6) || strcmp(name, 'Max_tracer')
            % initialize figure
            figure(fn+6)
            tr_ind = tr_ind + 1;
            if tr_ind == 1
                clf
                hold on
                leg_text = {};
            end
            max_tr = diagnos.(name{1}) / diagnos.(name{1})(1) - 1;
            plot(diagnos.Time, max_tr)
            xlabel('time (s)')
            ylabel('Tr_{max}/Tr(0) - 1')
            title('Maximum tracer')

            % Make legend
            if length(name{1}) > 7
                tr_num = name{1}(8);
                if strcmp(name, 'Max_tracer')
                    tr_num = {};
                end
            else
                tr_num = {};
            end
            leg_text = [leg_text,tr_num];%#ok
            % add legend
            if ~isempty(leg_text)
                legend(leg_text)
                legend('location','best')
            end

            %%%% Max (absolute value of) vorticity %%%%
        elseif strcmp(name, 'Max_vort_y')
            first_ind = 20; % remove values due to initial perturbations
            if length(diagnos.Max_vort_y) > first_ind
                figure(fn+7)
                subplot(2,1,1), cla reset
                inds = first_ind:length(diagnos.Max_vort_y);
                if params.ndims == 3
                    plot(diagnos.Time(inds),...
                        [diagnos.Max_vort_x(inds), diagnos.Max_vort_y(inds), diagnos.Max_vort_z(inds)])
                    set(gca,'yscale','log');
                    leg = legend({'Max $|\omega_x|$','Max $|\omega_y|$','Max $|\omega_z|$'},...
                        'Interpreter','Latex','FontSize',12);
                else
                    plot(diagnos.Time(inds), diagnos.Max_vort_y(inds),'Color',coly)
                    set(gca,'yscale','linear');
                    leg = legend({'Max $|\omega_y|$'}, 'Interpreter','Latex','FontSize',12);
                end
                leg.Location = 'best';
                leg.Box = 'off';
                ylabel('Maximum vorticity (1/s)')
                title('Maximum vorticity components')
            end
        elseif strcmp(name, 'Max_vort_x') || strcmp(name, 'Max_vort_z')
            continue

            %%%% Total Enstrophy components %%%%
        elseif strcmp(name, 'Enst_y_tot')
            first_ind = 20; % remove values due to initial perturbations
            if length(diagnos.Enst_y_tot) > first_ind
                figure(fn+7)
                subplot(2,1,2), cla reset
                inds = first_ind:length(diagnos.Enst_y_tot);
                if params.ndims == 3
                    plot(diagnos.Time(inds),...
                        [enst_x_tot(inds), enst_y_tot(inds), enst_z_tot(inds), enst_tot(inds)]/Vol)
                    set(gca,'yscale','log');
                    leg = legend({'$\Omega_x$','$\Omega_y$','$\Omega_z$','$\Omega_\mathrm{tot}$'},...
                        'Interpreter', 'Latex','FontSize',12);
                else
                    plot(diagnos.Time(inds), enst_y_tot(inds)/Vol, 'Color', coly)
                    set(gca,'yscale','linear');
                    leg = legend({'$\Omega_\mathrm{tot}$'},...
                        'Interpreter', 'Latex','FontSize',12);
                end
                leg.Location = 'best';
                leg.Box = 'off';
                xlabel('time (s)')
                ylabel('$\Omega_\mathrm{tot}/V$ (1/s$^2$)',...
                    'Interpreter','Latex','FontSize',12)
                title('Enstrophy components')
            end
        elseif strcmp(name, 'Enst_x_tot') || strcmp(name, 'Enst_z_tot')
            continue

        elseif strcmp(name, 'Diss_tot')
            %%%% Total Dissipation %%%%
            first_ind = 50;
            if length(diagnos.Diss_tot) < first_ind
                first_ind = 1;
            end
            inds = first_ind:length(diagnos.Diss_tot);
            figure(fn+8)
            subplot(3,1,2), cla reset
            plot(diagnos.Time(inds), diagnos.Diss_tot(inds))
            if compute_enstrophy
                set(gca, 'xticklabels', []);
            else
                xlabel('time (s)')
            end
            ylabel('$\epsilon_\mathrm{tot}$  (J/s)',...
                'Interpreter','Latex','FontSize',12)
            title('Total dissipation')

            %%%% Enstrophy-Dissipation ratio %%%%
            if compute_enstrophy
                figure(fn+8)
                subplot(3,1,1), cla reset
                plot(diagnos.Time(inds), enst_tot(inds)/Vol);
                set(gca, 'xticklabels', []);
                ylabel('$\Omega_\mathrm{tot}/V$ (1/s$^2$)',...
                    'Interpreter','Latex','FontSize',12)
                title('Total enstrophy')

                subplot(3,1,3), cla reset
                plot(diagnos.Time(inds), enst_diss(inds)-1,'.')
                xlabel('time (s)')
                ylabel('$\epsilon_{tot}/(2 \mu \Omega_{tot})-1$',...
                    'Interpreter','Latex','FontSize',12)
                title('Enstrophy-dissipation ratio error')
            end

        elseif strcmp(name, 'Properties')
            continue

        else
            disp([name,' not configured'])
            figure(fm), clf
            plot(diagnos.Time, diagnos.(name{1}))
            xlabel('time (s)')
            ylabel(strrep(name,'_',' '))
            fm = fm+1;
        end
    end

    %%%% Mixing efficiency %%%%
    if compute_BPE && compute_dissipation
        figure(fn+9)
        clf
        plot(time_rate, [mix_eff mix_cof])
        xlabel('time (s)')
        ylabel('Mixing')
        title('Mixing')
        yl = ylim;
        set(gca, 'YLim', [0 yl(2)]);
        leg = legend({'Mixing eff.','Mixing coeff.'});
        leg.Location = 'best';
        leg.Box = 'off';
    end
    %%
    if save_plots
        FigNames = {'SimRunRate', 'MaxVelsKE_Components', 'AllEnergyComponents', ...
            'AllEnergyComponentRates', 'MassDensityMaxDiss', 'AreaDeltaRho', ...
            'DyeTracer', 'MaxVortyTotEnstrophy', 'TotDissipation_Enstrophy',...
            'MixEfficiency'};
        
        for figs = 1:9
            if ishandle(fn+figs)
                figure(fn+figs);
                print(['Diag_', FigNames{figs}, '.png'], '-dpng');
            end
        end
    end
                
end
fprintf('\n')


all_diagnos.EnergeticChange.APE_Lost = AE_loss(end)/E0_and_W*100;
all_diagnos.EnergeticChange.APE_toMix = APE2BPE_tot(end)/E0_and_W*100;
all_diagnos.EnergeticChange.APE_toDiss = KE2Int_tot(end)/E0_and_W*100;

save('all_diagnos.mat', 'all_diagnos')
end

%% Volume of domain
function Vol = find_volume(gdpar)
    gd = gdpar.gd;
    params = gdpar.params;
    gdvec = get_vector_grid(gd);
    Vol = params.Lx * params.Ly * params.Lz;
    if strcmp(params.mapped_grid, 'true')
        if params.ndims == 3
            bot = gd.z(:,1,1);
            top = gd.z(:,1,params.Nz);
        else
            bot = gd.z(:,1);
            top = gd.z(:,params.Nz);
        end
        if strcmp(params.type_x, 'NO_SLIP')
            warning('Volume calculation is not setup for Cheb grid in x.')
        else
            Vol = Vol - params.Ly * trapz(gdvec.x, bot+top);
        end
    end
end

%% Minimum diffusivity
function kappa_min = find_min_diffu(params)
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

%% Kolmogorov and Batchelor Scales
function Scales = find_Kolm_Batch(diagnos, gdpar, kappa_min)
    % Kolmogorov Scale: eta = (rho_0 nu^3/epsilon)^(1/4)
    % Batchelor Scale:  lambda_B = eta * sqrt(kappa/nu)
    % rho_0 is included to make epsilon the dissipation per unit mass (makes dimensions work)

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
