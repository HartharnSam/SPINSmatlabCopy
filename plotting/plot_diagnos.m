%% plot quantities from diagnostic file (diagnostics.txt) 
%  and prints information (timing, length scales, and Energy)
%  about the run
%
% the Energy diagnostics are only true for a closed system with no external forcing
% ie. no tidal forcing, or forced shear flow
% Those energy transfer terms must be computed on their own

%%%%%%%%%%% Read Data %%%%%%%%%%%%%%%
% .mat is a cleaned version of .txt
% (run clean_diagnostics to remove overlapping regions)
% read diagnostic file
diag_file_name = 'diagnostics';
if exist([diag_file_name,'.mat'], 'file') == 2
    diag_file = [diag_file_name,'.mat'];
    diagnos = load(diag_file);
elseif exist([diag_file_name,'.txt'], 'file') == 2
    diag_file = [diag_file_name,'.txt'];
    diagnos = readtable(diag_file);
else
    error('Diagnostics file not found.')
end

% read grid and parameters from spins.conf
gdpar_vec = spins_gridparams('Vector',false);
params = gdpar_vec.params;

%%%%%%%%%%% Parse Data %%%%%%%%%%%%%%%
% shorten parameters
try
    clk_time = diagnos.Clock_time;
    sim_time = diagnos.Time;
catch
    error(['Diagnostics file incorrectly labeled the time (clock time and simulation time).',...
        newline,'Should be "Clock_time" and "Time"'])
end
rho_0 = params.rho_0;
visco = params.visco;

% find all diffusivity values
diffu_types = {'kappa','kappa_rho','kappa_tracer','kappa_dye','kappa_dye1','kappa_dye2',...
'kappa_T','kappa_S','kappa_t','kappa_s'};
for ii = 1:length(diffu_types)
    if isfield(params, diffu_types{ii})
        kappa(ii) = params.(diffu_types{ii});
    else
        kappa(ii) = NaN;
    end
end
kappa_min = min(kappa(:));

%% Parse time into step times and parse restarts
% find wall clock step time
clk_step_time = [clk_time(1); clk_time(2:end)-clk_time(1:end-1)];
clk_step_time(clk_step_time<0) = 0; % fix places where a restart occured
% find sim step time
sim_step_time = [sim_time(1); sim_time(2:end)-sim_time(1:end-1)];
sim_step_time(sim_step_time<0) = 0; % fix places where a restart occured
% check for restarts
restart_ind = find(diagnos.Iter == 1);
strt_ind = [restart_ind; length(diagnos.Iter)+1];
N_start = length(restart_ind);
N_restart = N_start - 1;
n_steps = length(diagnos.Iter);
% total steps and clock time
tot_clk_time = sum(clk_step_time);
tot_sim_time = sum(sim_step_time);

% legend labels for each run
run_label = cell(1,N_start);
for ii = 1:N_start
    run_label{ii} = ['Run ',num2str(ii)];
end

%% max density variation
try
    rho_init = diagnos.Max_density(1);
    rho_var = diagnos.Max_density/rho_init - 1;
catch
    temp_init = diagnos.Max_temperature(1);
    salt_init = diagnos.Max_salinity(1);
    temp_var = diagnos.Max_temperature/temp_init;
    salt_var = diagnos.Max_salinity/salt_init;
    rho_init = eqn_of_state(temp_init, salt_init);
    rho_var = eqn_of_state(temp_var, salt_var)/rho_init - 1;
end


%%%%%%%%%%% Print Timing info %%%%%%%%%%%%%%%
% functions for converting seconds into other units
s2ms   = @(sec) datestr(datenum(0,0,0,0,0,sec),'MM:SS');
s2hms  = @(sec) datestr(datenum(0,0,0,0,0,sec),'HH:MM:SS');
s2dhm  = @(sec) datestr(datenum(0,0,0,0,0,sec),'DD:HH:MM');
s2dhms = @(sec) datestr(datenum(0,0,0,0,0,sec),'DD:HH:MM:SS');
% print info
fprintf('\n')
disp('---- Timings ----')
disp(['Total simulated time:  ',num2str(tot_sim_time),' s'])
disp(['Most recent sim. time: ',num2str(sim_time(end)),' s'])
disp(['Total clock time:      ',num2str(s2dhm(tot_clk_time)),' (D:H:M)'])
%% print out average write time and step time
warning('off','MATLAB:table:ModifiedAndSavedVarnames'); % suppress warning for this only
if exist('plot_times.txt', 'file')
    % if outputs exist, remove output write time
    plottimes = readtable('plot_times.txt');
    avg_write = mean(plottimes.WriteTime_s_);
    tot_write =  sum(plottimes.WriteTime_s_);
    avg_clk_step = (tot_clk_time - tot_write)/n_steps;
    avg_sim_step = tot_sim_time/n_steps;
    clk_per_sim = avg_clk_step/avg_sim_step;
    if avg_write < 1
        fprintf('Average write time: %8.3f s\n',avg_write)
    else
        fprintf('Average write time: %s (M:S)\n',s2ms(avg_write))
    end
    fprintf('Avg. sim. step time: %7.3f s\n',avg_sim_step)
    if avg_clk_step < 1
        fprintf('Avg. clock step time:  %5.3f s\n',avg_clk_step)
    else
        fprintf('Avg. clock step time:  %s (M:S)\n',s2ms(avg_clk_step))
    end
    fprintf('Avg. clock time per sim. sec.: %s (D:H:M:S)\n',s2dhms(clk_per_sim))
else
    clk_per_sim = tot_clk_time/tot_sim_time;
    fprintf('Avg. clock time per sim. sec.: %s (D:H:M:S)\n',s2dhms(clk_per_sim))
    disp('plot_times.txt file not found, or incorrectly configured.');
end
warning('on','MATLAB:table:ModifiedAndSavedVarnames');

%%%%%%%%%%% Shorten entstrophy variable and compare against dissipation %%%%%%%%%%%%%%%
% from Mike Waite's Turbulence class: 
%   total diss. = 2*mu*(total enstrophy)
%   when boundary conditions are not no-slip (ie. free-slip or periodic on all sides)
%   see also Yeung et al. 2012 in JFM
% check how close the ratio is to 1
if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot')) ...
        && any(ismember(diagnos.Properties.VariableNames, 'Enst_y_tot'))
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
    enst_diss = diagnos.Diss_tot./(enst_tot*2*visco*rho_0);
end

%%%%%%%%%%% Print Kolmogorov and Batchelor scales %%%%%%%%%%%%%%%
% eta = (rho_0 nu^3/epsilon)^(1/4)
% lambda_B = eta * sqrt(kappa/nu)
% rho_0 is included to make epsilon the dissipation per unit mass (makes dimensions work)
if any(ismember(diagnos.Properties.VariableNames, 'Max_diss'))
    % find max diss. ignoring the first 100 points containing the random perturbations
    if length(diagnos.Diss_tot) >= 100
        max_diss = max(diagnos.Max_diss(100:end));
    else
        max_diss = max(diagnos.Max_diss);
    end
    Kolm = (rho_0*visco^3/max_diss)^(1/4);
    Batch = Kolm*sqrt(kappa_min/visco);
    % compare min grid size to these scales
    if strcmp(params.type_z, 'NO_SLIP')
        params.dz = max(gdpar_vec.gd.z(2:end) - gdpar_vec.gd.z(1:end-1));
    end
    if params.ndims == 3
        max_dxyz = max([params.dx,params.dy,params.dz]);
    else
        max_dxyz = max([params.dx,params.dz]);
    end
    dx_Kolm  = max_dxyz/Kolm;
    dx_Batch = max_dxyz/Batch;
    fprintf('\n')
    disp('---- Kolmogorov and Batchelor Scales ----')
    disp(['dx/eta =      ',num2str(dx_Kolm)])
    disp(['dx/lambda_B = ',num2str(dx_Batch)])
end

%%%%%%%%%%% Parse and Print Energy diagnostics %%%%%%%%%%%%%%%
% compute total KE
if params.ndims == 2
    diagnos.KE_y = diagnos.KE_x*0;
end
KE_tot = diagnos.KE_x + diagnos.KE_y + diagnos.KE_z;
E_tot  = KE_tot + diagnos.PE_tot;
E_loss = E_tot(1) - E_tot;
% compute rates
time_rate = linspace(sim_time(1), sim_time(end),round(length(sim_time)/2));
Dmat = FiniteDiff(time_rate, 1, 2, true);
KE_rate = Dmat*interp1(sim_time, KE_tot, time_rate, 'pchip')';
E_rate  = Dmat*interp1(sim_time, E_tot,  time_rate, 'pchip')';

% compute APE, total Avail. Energy (AE_tot), and change in BPE, and their rates of change
if any(ismember(diagnos.Properties.VariableNames, 'BPE_tot'))
    APE_tot = diagnos.PE_tot - diagnos.BPE_tot;
    AE_tot  = KE_tot + APE_tot;
    AE_loss = AE_tot(1) - AE_tot;
    BPE_change = diagnos.BPE_tot - diagnos.BPE_tot(1);
    BPE_rate = Dmat*interp1(sim_time, BPE_change, time_rate, 'pchip')';
    APE_rate = Dmat*interp1(sim_time, APE_tot,    time_rate, 'pchip')';
    AE_rate  = Dmat*interp1(sim_time,  AE_tot,    time_rate, 'pchip')';
end
% energy transfered from internal to BPE (phi_i)
if any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
    Int2BPE_rate = diagnos.BPE_from_int;
    Int2BPE_tot  = cumtrapz(diagnos.Time, Int2BPE_rate);
    % energy transfered from APE to BPE (phi_a = phi_d - phi_i)
    if any(ismember(diagnos.Properties.VariableNames, 'BPE_tot'))
        APE2BPE_tot  = BPE_change - Int2BPE_tot;
        APE2BPE_rate = Dmat*interp1(sim_time, APE2BPE_tot, time_rate, 'pchip')';
    end
end
% energy dissipation rates and cumulative amounts
if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot'))
    % dissipation (epsilon)
    KE2Int_rate = diagnos.Diss_tot;
    KE2Int_tot  = cumtrapz(diagnos.Time, KE2Int_rate);
    % KE to APE (phi_z)
    KE2APE_tot  = -(KE_tot - KE_tot(1)) - KE2Int_tot;
    KE2APE_rate = Dmat*interp1(sim_time, KE2APE_tot, time_rate, 'pchip')';
    % offset onto the grid to use with other rates which are off_set by 1st order FD derivative
    diss_offset = interp1(sim_time, diagnos.Diss_tot, time_rate, 'pchip')';
    % compute rate of change of internal energy, and the cumulative change 
    if any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
        Int_rate = diagnos.Diss_tot - diagnos.BPE_from_int;
        Int_change = cumtrapz(diagnos.Time, Int_rate);
    end
end

% compute energy created/removed by numerics (filter and numerical errors/scheme etc.)
% only if dissipation, BPE, and BPE_from_int were calculated
if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot')) ...
        && any(ismember(diagnos.Properties.VariableNames, 'BPE_tot')) ...
        && any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
    NumE_tot  = E_loss - Int_change;
    NumE_rate = Dmat*interp1(sim_time, NumE_tot, time_rate, 'pchip')';

    fprintf('\n')
    fprintf('---- Energy ----\n')
    fprintf('As a percentage of the initial available energy\n')
    fprintf('Total avail. energy lost: %6.2f %%\n',AE_loss(end)/AE_tot(1)*100)
    fprintf('  - Total energy lost:    %6.2f %%\n',E_loss(end)/AE_tot(1)*100)
    fprintf('     > to dissipation:    %6.2f %%\n',KE2Int_tot(end)/AE_tot(1)*100)
    fprintf('     > to numerics:       %6.2f %%\n',NumE_tot(end)/AE_tot(1)*100)
    fprintf('     > to internal:       %6.2f %%\n',-Int2BPE_tot(end)/AE_tot(1)*100)
    fprintf('  - Total BPE gained:     %6.2f %%\n',BPE_change(end)/AE_tot(1)*100)
    fprintf('     > from APE:          %6.2f %%\n',APE2BPE_tot(end)/AE_tot(1)*100)
    fprintf('     > from internal:     %6.2f %%\n',Int2BPE_tot(end)/AE_tot(1)*100)
end

% compute mixing efficiencies
if any(ismember(diagnos.Properties.VariableNames, 'BPE_tot'))
    inds = 20:length(BPE_rate);
    fprintf('\n')
    fprintf('---- Mixing efficiency ----\n')
    fprintf('Change in BPE / |change in AE|: %6.2f \n', ...
        BPE_change(end)/abs(AE_tot(end) - AE_tot(1)))
    if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot'))
        mix_eff1 = BPE_rate(inds)./diss_offset(inds); 
        max_mix1 = max(mix_eff1);
        mix_eff2 = BPE_rate(inds)./(diss_offset(inds) + BPE_rate(inds));
        max_mix2 = max(mix_eff2);
        mix_eff3 = APE2BPE_rate(inds)./(diss_offset(inds) + APE2BPE_rate(inds));
        max_mix3 = max(mix_eff3);
        fprintf('Max {phi_d/epsilon}:            %6.2f \n', max_mix1)
        fprintf('Max {phi_d/(epsilon+phi_d)}:    %6.2f \n', max_mix2)
        fprintf('Max {phi_a/(epsilon+phi_a)}:    %6.2f \n', max_mix3)
    end
end


%%%%%%%%%%% make plots %%%%%%%%%%%%%%%
fn = 20;    % first figure number
fm = fn+15; % first figure number for unknown diagnostics
cols = get(groot,'DefaultAxesColorOrder');
coly = cols(2,:);

for name = diagnos.Properties.VariableNames
    if strcmp(name, 'Iter')
        continue

    %%%% Clock Time per Step %%%%
    elseif strcmp(name, 'Clock_time')
        figure(fn), clf
        subplot(3,1,1), cla
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
            ylim([yl(1) prctile(clk_step_time,99.98)])
        end
        legend(run_label)
        legend('location','best')
        %legend('boxoff')
        box on
        hold off

    %%%% Simulation Time per Step %%%%
    elseif strcmp(name, 'Time')
        figure(fn)
        subplot(3,1,2), cla
        hold on
        for ii = 1:N_start
            its = strt_ind(ii):strt_ind(ii+1)-1;
            plot(diagnos.Iter(its) + strt_ind(ii)-1, sim_step_time(its), '.')
        end
        ylabel('\Delta T_s (s)')
        title('Simulation time step')
        if length(sim_step_time) > 500
            ylim([0 1.2*max(sim_step_time(500:end))])
        end
        box on
        hold off

        %%%% Clock Time per Simulation Time per Step %%%%
        subplot(3,1,3), cla
        hold on
        clk_per_sim = clk_step_time./sim_step_time;
        for ii = 1:N_start
            its = strt_ind(ii):strt_ind(ii+1)-1;
            plot(diagnos.Iter(its) + strt_ind(ii)-1, clk_per_sim(its), '.')
        end
        xlabel('Iteration')
        ylabel('\Delta{T_c} / \Delta{T_s}')
        title('Ratio of time steps')
        clk_per_sim_num = clk_per_sim(~isinf(clk_per_sim) & ~isnan(clk_per_sim));
        ylim([0 2.5*mean(clk_per_sim_num)]); 
        box on
        hold off


    %%%% Maximum (absolute value of) Velocities %%%%
    elseif strcmp(name, 'Max_u')
        figure(fn+1)
        subplot(2,1,1), cla
        if params.ndims == 2
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
        subplot(2,1,2), cla
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
        subplot(2,1,1), cla
        hold on
        plot(diagnos.Time, KE_tot)
        energy_label = {'KE'};
        % Add APE, AE, and change in BPE, if calculated
        if any(ismember(diagnos.Properties.VariableNames, 'BPE_tot'))
            plot(diagnos.Time, [APE_tot, AE_tot, BPE_change])
            energy_label = [energy_label, {'APE','AE','\Delta{BPE}'}];
        end
        % Add change in internal energy if calculated
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot')) && ...
                any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(diagnos.Time, Int_change,'k')
            energy_label = [energy_label, '\Delta{Int}'];
        end
        % add energy created/removed by numerics (filter and numerical errors/scheme etc.)
        % only if dissipation, BPE, and BPE_from_int were calculated
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot')) && ...
                any(ismember(diagnos.Properties.VariableNames, 'BPE_tot')) && ...
                any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(diagnos.Time, NumE_tot)
            %plot([0 diagnos.Time(end)],[1 1]*AE_tot(1),'k')
            %energy_label = [energy_label, {'Numerics','BPE from int','AE(0)'}];
            energy_label = [energy_label, {'Numerics'}];
        end
        grid on
        box on
        ylabel('Energy (J)')
        title('Energy components')
        legend(energy_label)
        legend('location','best')
        legend('boxoff')

        %% plot the rate of change of Energy components
        subplot(2,1,2), cla
        hold on
        plot(time_rate, KE_rate)
        energy_label = {'KE'};
        % Add APE, AE, and change in BPE, if calculated
        if any(ismember(diagnos.Properties.VariableNames, 'BPE_tot'))
            plot(time_rate, [APE_rate, AE_rate, BPE_rate])
            energy_label = [energy_label, {'APE','AE','BPE (\phi_d)'}];
        end
        % Add change in internal energy if calculated
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot')) && ...
                any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(diagnos.Time, Int_rate,'k')
            energy_label = [energy_label, 'Internal'];
        end
        % add energy created/removed by numerics (filter and numerical errors/scheme etc.)
        % only if dissipation, BPE, and BPE_from_int were calculated
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot')) && ...
                any(ismember(diagnos.Properties.VariableNames, 'BPE_tot')) && ...
                any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(time_rate, NumE_rate)
            %plot([0 diagnos.Time(end)],[1 1]*AE_tot(1),'k')
            %energy_label = [energy_label, {'Numerics','BPE from int','AE(0)'}];
            energy_label = [energy_label, {'Numerics'}];
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
        subplot(2,1,1), cla
        hold on
        % interpolate onto a regular grid of similar size
        energy_label = {};
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot'))
            plot(diagnos.Time, KE2APE_tot)
            energy_label = {'KE to APE (\int \phi_z dt)'};
        end
        % Add APE if calculated
        if any(ismember(diagnos.Properties.VariableNames, 'BPE_tot')) ...
                && any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(diagnos.Time, APE2BPE_tot)
            energy_label = [energy_label, {'APE to BPE (\int \phi_a dt)'}];
        end
        % Add energy converted from internal energy
        if any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(diagnos.Time, Int2BPE_tot,'Color',cols(4,:))
            energy_label = [energy_label, 'Int to BPE (\int \phi_i dt)'];
        end
        % Add dissipation if calculated
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot'))
            plot(diagnos.Time(8:end), KE2Int_tot(8:end),'k')
            energy_label = [energy_label, 'KE to Int (diss, \int \epsilon dt)'];
        end
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot')) ...
                && any(ismember(diagnos.Properties.VariableNames, 'BPE_tot')) ...
                && any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(diagnos.Time, NumE_tot,'Color',cols(5,:))
            energy_label = [energy_label, 'to Numerics'];
        end
        grid on
        ylabel('Energy (J)')
        title('Energy converted')
        legend(energy_label)
        legend('location','best')
        legend('boxoff')
        box on
        hold off

        subplot(2,1,2), cla
        hold on
        % interpolate onto a regular grid of similar size
        energy_label = {};
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot'))
            plot(time_rate(4:end), KE2APE_rate(4:end))
            energy_label = {'KE  to APE (\phi_z)'};
        end
        % Add APE if calculated
        if any(ismember(diagnos.Properties.VariableNames, 'BPE_tot')) ...
                && any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(time_rate(4:end), APE2BPE_rate(4:end))
            energy_label = [energy_label, {'APE to BPE (\phi_a)'}];
        end
        % Add energy converted from internal energy
        if any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(diagnos.Time(4:end), diagnos.BPE_from_int(4:end),'Color',cols(4,:))
            energy_label = [energy_label, 'Int to BPE (\phi_i)'];
        end
        % Add dissipation if calculated
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot'))
            plot(diagnos.Time(8:end), diagnos.Diss_tot(8:end),'k')
            energy_label = [energy_label, 'KE  to Int (diss, \epsilon)'];
        end
        if any(ismember(diagnos.Properties.VariableNames, 'Diss_tot')) ...
                && any(ismember(diagnos.Properties.VariableNames, 'BPE_tot')) ...
                && any(ismember(diagnos.Properties.VariableNames, 'BPE_from_int'))
            plot(time_rate(4:end), NumE_rate(4:end),'Color',cols(5,:))
            energy_label = [energy_label, 'to Numerics'];
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
        strcmp(name, 'PE_tot') || strcmp(name, 'BPE_from_int')
        continue

    %%%% Total Mass %%%%
    elseif strcmp(name, 'Mass')
        figure(fn+4)
        subplot(3,1,1), cla
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
        subplot(3,1,2), cla
        plot(diagnos.Time, rho_var)
        ylabel('$\rho_\mathrm{max}/\rho_\mathrm{max}(0) - 1$',...
            'Interpreter','Latex','FontSize',14)
        title('Max tracers deviation (Density, etc.)')
    elseif strcmp(name, 'Max_temp') || strcmp(name, 'Max_salt')
        figure(fn+4)
        subplot(2,1,2), cla
        plot(diagnos.Time, [temp_var, salt_var])
        ylabel('Tracer ratio')
        legend({'$\rho_\mathrm{max}/\rho_\mathrm{max}(0)$',...
            'T_{max}/T_{max}(0)','S_{max}/S_{max}(0)'},...
            'Interpreter','Latex','FontSize',14)
        legend('location','best')
        legend('boxoff')

    %%%% Max Dissipation %%%%
    elseif strcmp(name, 'Max_diss')
        first_ind = 20;
        if length(diagnos.Max_diss) > first_ind
            figure(fn+4)
            subplot(3,1,3), cla
            inds = first_ind:length(diagnos.Max_diss);
            plot(diagnos.Time(inds), diagnos.Max_diss(inds))
            xlabel('time (s)')
            ylabel('$\epsilon_\mathrm{max}$  (J/s)',...
                'Interpreter','Latex','FontSize',14) 
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
    elseif strcmp(name, 'Max_dye') || strcmp(name, 'Max_tracer')
        figure(fn+6), clf
        plot(diagnos.Time, diagnos.(name{1}))
        xlabel('time (s)')
        ylabel('Maximum tracer') 

    %%%% Max (absolute value of) vorticity %%%%
    elseif strcmp(name, 'Max_vort_y')
        first_ind = 20; % remove values due to initial perturbations
        if length(diagnos.Max_vort_y) > first_ind
            figure(fn+7)
            subplot(2,1,1), cla
            inds = first_ind:length(diagnos.Max_vort_y);
            if params.ndims == 3
                % need to add in max vorticity component max{ om = sqrt(om_x^2 + om_y^2 + om_z^2)}
                plot(diagnos.Time(inds),...
                [diagnos.Max_vort_x(inds), diagnos.Max_vort_y(inds), diagnos.Max_vort_z(inds)])
                set(gca,'yscale','log');
                leg = legend({'Max $|\omega_x|$','Max $|\omega_y|$','Max $|\omega_z|$'},...
                    'Interpreter','Latex','FontSize',14);
            else
                plot(diagnos.Time(inds), diagnos.Max_vort_y(inds),'Color',coly)
                set(gca,'yscale','linear');
                leg = legend({'Max $|\omega_y|$'}, 'Interpreter','Latex','FontSize',14);
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
            subplot(2,1,2), cla
            inds = first_ind:length(diagnos.Enst_y_tot);
            if params.ndims == 3
                plot(diagnos.Time(inds),...
                [enst_x_tot(inds), enst_y_tot(inds), enst_z_tot(inds), enst_tot(inds)])
                set(gca,'yscale','log');
                leg = legend({'$\Omega_x$','$\Omega_y$','$\Omega_z$','$\Omega_\mathrm{tot}$'},...
                    'Interpreter', 'Latex','FontSize',14);
            else
                plot(diagnos.Time(inds), enst_y_tot(inds), 'Color', coly)
                set(gca,'yscale','linear');
                leg = legend({'$\Omega_\mathrm{tot}$'},...
                    'Interpreter', 'Latex','FontSize',14);
            end
            leg.Location = 'best';
            leg.Box = 'off';
            xlabel('time (s)')
            ylabel('$\Omega_\mathrm{tot}$ (1/s$^2$)',...
                'Interpreter','Latex','FontSize',14)
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
        subplot(3,1,2), cla
        plot(diagnos.Time(inds), diagnos.Diss_tot(inds))
        if any(ismember(diagnos.Properties.VariableNames, 'Enst_y_tot'))
            xticklabels([])
        else
            xlabel('time (s)')
        end
        ylabel('$\epsilon_\mathrm{tot}$  (J/s)',...
            'Interpreter','Latex','FontSize',14) 
        title('Total dissipation')

        %%%% Enstrophy-Dissipation ratio %%%%
        if any(ismember(diagnos.Properties.VariableNames, 'Enst_y_tot'))
            figure(fn+8)
            subplot(3,1,1), cla
            plot(diagnos.Time(inds), enst_tot(inds));
            xticklabels([])
            ylabel('$\Omega_\mathrm{tot}$ (1/s$^2$)',...
                'Interpreter','Latex','FontSize',14)
            title('Total enstrophy')

            subplot(3,1,3), cla
            plot(diagnos.Time(inds), enst_diss(inds)-1,'.')
            xlabel('time (s)')
            ylabel('$\epsilon_{tot}/(2 \mu \Omega_{tot})-1$',...
                'Interpreter','Latex','FontSize',14)
            title('Enstrophy-dissipation ratio error')
        end

    else
        disp([name,' not configured'])
        figure(fm), clf
        plot(diagnos.Time, diagnos.(name{1}))
        xlabel('time (s)')
        ylabel(strrep(name,'_',' '))
        fm = fm+1;
    end
end
