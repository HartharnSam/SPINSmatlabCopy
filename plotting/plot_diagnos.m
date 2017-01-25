%% plot quantities from diagnostic files (analysis.txt and enstrophy.txt) 
%% and prints timing information from plot_times.txt

%%%%%%%%%%% Read Data %%%%%%%%%%%%%%%
% .mat is a cleaned version of .txt (could have doubled times)
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
% read enstrophy file
enst_file_name = 'enstrophy';
if exist([enst_file_name,'.mat'], 'file') == 2
    enst_file = [enst_file_name,'.mat'];
    enst = load(enst_file);
    is_enst = true;
elseif exist([enst_file_name,'.txt'], 'file') == 2
    enst_file = [enst_file_name,'.txt'];
    enst = readtable(enst_file);
    is_enst = true;
else
    warning('Enstrophy file not found.')
    is_enst = false;
end

% read spins.conf
gdpar = spins_gridparams('Vector',false);
params = gdpar.params;
clear gdpar

%%%%%%%%%%% Parse Data %%%%%%%%%%%%%%%
%  shorten enstrophy components
try
    enst_x = enst.enst_x;
    enst_y = enst.enst_y;
    enst_z = enst.enst_z;
    enst_tot = enst.enst_tot;
catch
    disp('Enstrophy file incorrectly configured.')
end
% shorten other parameters
try
    clk_time = diagnos.Clock_time;
    sim_time = diagnos.Sim_time;
catch
    error('Diagnostics file incorrectly configured.')
end
if params.ndims == 2
    Vol = params.Lx*params.Lz;
else
    Vol = params.Lx*params.Ly*params.Lz;
end
rho_0 = params.rho_0;
visco = params.visco;
% find all diffusivity values
diffu_types = {'kappa','kappa_rho','kappa_tracer','kappa_dye','kappa_dye1','kappa_dye2',...
               'kappa_T','kappa_S'};
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
sim_step_time = [0; sim_time(2:end)-sim_time(1:end-1)];
sim_step_time(sim_step_time<0) = 0; % fix places where a restart occured
% check for restarts
restart_ind = find(diagnos.Iter == 1);
strt_ind = [restart_ind; length(diagnos.Iter)+1];
N_start = length(restart_ind);
N_restart = N_start - 1;
n_steps = length(diagnos.Iter);
% total steps and clock time
tot_clk_time = sum(clk_step_time) + clk_time(1);
tot_sim_time = sum(sim_step_time);

%% print out time information
% functions for times
s2ms   = @(sec) datestr(datenum(0,0,0,0,0,sec),'MM:SS');
s2hms  = @(sec) datestr(datenum(0,0,0,0,0,sec),'HH:MM:SS');
s2dhm  = @(sec) datestr(datenum(0,0,0,0,0,sec),'DD:HH:MM');
s2dhms = @(sec) datestr(datenum(0,0,0,0,0,sec),'DD:HH:MM:SS');
% print info
disp(['Total simulated time:  ',num2str(tot_sim_time),' s'])
disp(['Most recent sim. time: ',num2str(sim_time(end)),' s'])
disp(['Total physical time:   ',num2str(s2dhm(tot_clk_time)),' (D:H:M)'])
% print out average write time
warning('off','MATLAB:table:ModifiedVarnames'); % suppress warning for this only
try
    % print a table for each run (include the number of processors)
    plottimes = readtable('plot_times.txt');
    avg_write = mean(plottimes.WriteTime_s_);
    tot_write =  sum(plottimes.WriteTime_s_);
    avg_clk_step = (tot_clk_time - tot_write)/n_steps;
    avg_sim_step = tot_sim_time/n_steps;
    clk_per_sim = avg_clk_step/avg_sim_step;
    if avg_write < 1
        disp(['Average write time:    ',num2str(avg_write),' s'])
    else
        disp(['Average write time:    ',num2str(s2ms(avg_write)),' (M:S)'])
    end
    disp(['Avg. sim. step time:     ',num2str(avg_sim_step),' s'])
    if avg_clk_step < 1
        disp(['Avg. physical step time: ',num2str(avg_clk_step),' s'])
    else
        disp(['Avg. physical step time: ',num2str(s2ms(avg_clk_step)),' (M:S)'])
    end
    disp(['Avg. phys. time per sim. sec.: ',num2str(s2dhms(clk_per_sim)),' (D:H:M:S)'])
catch
    clk_per_sim = tot_clk_time/tot_sim_time;
    disp(['Avg. phys. time per sim. sec.: ',num2str(s2dhms(clk_per_sim)),' (D:H:M:S)'])
    disp('plot_times.txt file not found, or incorrectly configured.');
end
warning('on','MATLAB:table:ModifiedVarnames');

%% print Kolmogorov and Batchelor scales
% eta = (rho_0 nu^3/epsilon)^(1/4)
% rho_0 is included to make epsilon the dissipation per unit mass (makes dimensions work)
% find max diss. ignoring the first 100 points containing the random perturbations
if length(diagnos.Total_dissipation) >= 100
    max_diss = max(diagnos.Total_dissipation(100:end));
else
    max_diss = max(diagnos.Total_dissipation);
end
% if the dissipation was calculated, print info
if max_diss > 0
    Kolm = (rho_0*visco^3/max_diss)^(1/4);
    Batch = Kolm*sqrt(kappa_min/visco);
    if strcmp(params.mapped_grid, 'true')
        if params.ndims == 3
            max_dxyz = max([params.dx,params.dy]);
        else
            max_dxyz = max([params.dx]);
        end
    elseif strcmp(params.mapped_grid, 'false')
        if params.ndims == 3
            max_dxyz = max([params.dx,params.dy, params.dz]);
        else
            max_dxyz = max([params.dx,params.dz]);
        end
    end
    dx_Kolm = max_dxyz/Kolm;
    dx_Batch = max_dxyz/Batch;
    disp(['dx/eta =      ',num2str(dx_Kolm)])
    disp(['dx/lambda_B = ',num2str(dx_Batch)])
end

% max density variation
try
    rho_init = diagnos.Max_density(1);
    rho_var = diagnos.Max_density/rho_init;
catch
    temp_init = diagnos.Max_temperature(1);
    salt_init = diagnos.Max_salinity(1);
    temp_var = diagnos.Max_temperature/temp_init;
    salt_var = diagnos.Max_salinity/salt_init;
end

% legend labels for each run
lab = cell(1,N_start);
for ii = 1:N_start
    lab{ii} = ['Run ',num2str(ii)];
end


%%%%%%%%%%% make plots %%%%%%%%%%%%%%%
n = 20;    % first figure number

%%%% Clock Time per Step %%%%
figure(n), clf
hold on
for ii = 1:N_start
    if ii == 1
        its = (strt_ind(ii)+1):strt_ind(ii+1)-1;
    else
        its = strt_ind(ii):strt_ind(ii+1)-1;
    end
    plot(diagnos.Iter(its) + strt_ind(ii)-1, clk_step_time(its), '.')
end
xlabel('Iteration')
ylabel('Clock Time per step (s)')
legend(lab)
legend('location','best')
legend('boxoff')
%---- add simulation time axis on top ----%
% shift and create additional axis
ax1 = gca;
ax1_pos = ax1.Position;
ax1_pos(4) = ax1_pos(4) - 0.04;
ax1.Position = ax1_pos;
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(diagnos.Iter, clk_step_time,'Parent',ax2, 'LineStyle', 'none')
% get matching ticks
ax2.XLim = ax1.XLim;
xticks = ax2.XTick(1:end-1); % drop last one
ax2.XTick = xticks;         % make new ticks
%% make new tick labels
xticks(1) = 1;
xticks = sim_time(xticks);
if xticks(1) < 0.01
    xticks(1) = 0;
end
xlab = cell(1,length(xticks));
for ii = 1:length(xticks);
    xlab{ii} = num2str(xticks(ii),'%.3g');
end
% place on axis
ax2.XTickLabel = xlab;
ax2.YTick = [];
ax2.TickDir = 'Out';
xlabel('Sim time (s)')
axes(ax1)
box on
hold off
pause(0.25) % graphic errors occur if not paused - not sure why

%%%% Simulation Time per Step %%%%
n = n+1;
figure(n), clf
hold on
for ii = 1:N_start
    its = strt_ind(ii):strt_ind(ii+1)-1;
    plot(diagnos.Iter(its) + strt_ind(ii)-1, sim_step_time(its), '.')
end
xlabel('Iteration')
ylabel('Simulation Time per step (s)')
legend(lab)
legend('location','best')
legend('boxoff')
%---- add simulation time axis on top ----%
% shift and make new axis
ax1 = gca;
ax1_pos = ax1.Position;
ax1_pos(4) = ax1_pos(4) - 0.04;
ax1.Position = ax1_pos;
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(diagnos.Iter, sim_step_time,'Parent',ax2, 'LineStyle', 'none')
% get matching ticks
ax2.XLim = ax1.XLim;
xticks = ax2.XTick(1:end-1); % drop last one
ax2.XTick = xticks;         % make new ticks
% place on axis
ax2.XTickLabel = xlab;
ax2.YTick = [];
ax2.TickDir = 'Out';
xlabel('Sim time (s)')
axes(ax1)
box on
hold off

%%%% Maximum Velocities %%%%
n = n+1;
figure(n)
plot(diagnos.Sim_time, [diagnos.Max_U, diagnos.Max_V, diagnos.Max_W, diagnos.Max_vel])
if params.ndims == 3
    set(gca,'yscale','log');
else
    set(gca,'yscale','linear');
end
xlabel('time (s)')
ylabel('Maximum velocity (m/s)')
legend('u','v','w','u_{vec}')
legend('location','best')
legend('boxoff')

%%%% Maximum Density of Salt/Temp %%%%
n = n+1;
figure(n)
try
    plot(diagnos.Sim_time, rho_var)
    ylabel('$\rho_\mathrm{max}/\rho_\mathrm{max}(0)$','Interpreter','Latex')
catch
    plot(diagnos.Sim_time, [temp_var, salt_var])
    ylabel('Tracer ratio')
    legend({'T_{max}/T_{max}(0)','S_{max}/S_{max}(0)'})
    legend('location','best')
    legend('boxoff')
end
xlabel('time (s)')

%%%% KE components %%%%
n = n+1;
figure(n)
plot(diagnos.Sim_time, [diagnos.KE_x, diagnos.KE_y, diagnos.KE_z, diagnos.Total_KE]/Vol)
if params.ndims == 3
    set(gca,'yscale','log');
else
    set(gca,'yscale','linear');
end
xlabel('time (s)')
ylabel('KE$_\mathrm{tot}$ / V (J/m$^3$)','Interpreter','Latex')
legend('KE_x', 'KE_y', 'KE_z', 'KE_{tot}')
legend('location','best')
legend('boxoff')

%%%% Total PE %%%%
n = n+1;
figure(n)
plot(diagnos.Sim_time, diagnos.Total_PE/Vol)
xlabel('time (s)')
ylabel('PE$_\mathrm{tot} / V$ (J/m$^3$)','Interpreter','Latex')

%%%% Enstrophy components %%%%
n = n+1;
enst_start = 20;
if is_enst && length(enst_x) > enst_start
    figure(n)
    enst_inds = enst_start:length(enst_x);
    plot(enst.Time(enst_inds),...
    [enst_x(enst_inds), enst_y(enst_inds), enst_z(enst_inds), enst_tot(enst_inds)]/Vol)
    if params.ndims == 3
        set(gca,'yscale','log');
    else
        set(gca,'yscale','linear');
    end
    legend({'Enst_x','Enst_y','Enst_z','Enst_{tot}'})
    legend('location','best')
    legend('boxoff')
    xlabel('time (s)')
    ylabel('$\Omega_\mathrm{tot}$ / V (1/s$^2$)','Interpreter','Latex')
    title('Total Enstrophy')
else
    clf
    disp('Not enough outputs for enstrophy')
end

%%%% Dissipation %%%%
n = n+1;
figure(n)
if length(diagnos.Total_dissipation) > enst_start && ...
    max(diagnos.Total_dissipation) > 0
    plot(diagnos.Sim_time(enst_start:end), diagnos.Total_dissipation(enst_start:end))
    xlabel('time (s)')
    ylabel('$\epsilon_\mathrm{tot} / V$  (J/s /m$^3$)','Interpreter','Latex') 
    title('Total Dissipation')

%%%% Enstrophy-Dissipation ratio %%%%
% from Mike Waite's Turbulence class: total diss. = 2*mu*(total enstrophy)
% check how close the ratio is to 1
    len_diss = length(diagnos.Total_dissipation);
    len_enst = length(enst.enst_tot);
    len_d_e = min(len_diss, len_enst);
    enst_diss = diagnos.Total_dissipation(1:len_d_e)./(enst.enst_tot(1:len_d_e)*(2*visco*rho_0));
    n = n+1;
    figure(n)
    plot(diagnos.Sim_time(1:len_d_e), enst_diss-1,'.')
    xlabel('time (s)')
    ylabel('$\epsilon_{tot}/(2 \mu \Omega_{tot})-1$','Interpreter','Latex')
    title('Enstrophy-dissipation ratio error')
else
    clf
    disp('Not enough outputs for dissipation, or dissipation not computed')
end

%%%% Maximum of Dye or Tracer %%%%
if  ~isempty(strmatch('Max_dye'   ,diagnos.Properties.VariableNames, 'exact')) || ...
    ~isempty(strmatch('Max_tracer',diagnos.Properties.VariableNames, 'exact'))
    n = n+1;
    figure(n)
    try
        plot(diagnos.Sim_time, diagnos.Max_dye)
    catch
        plot(diagnos.Sim_time, diagnos.Max_tracer)
    end
    xlabel('time (s)')
    ylabel('Maximum tracer') 
else
    try
        close (n+1)
    end
end
