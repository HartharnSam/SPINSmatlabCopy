% plots things from the stresses diagnostic file

% read parameters
par = spins_params();

% read the file
if isfield(par, 'compute_stresses_bottom')
    stress_bot = readtable('stresses_bottom.txt');
end
if isfield(par, 'compute_stresses_top')
    stress_top = readtable('stresses_top.txt');
end
% error check
if ~isfield(par, 'compute_stresses_top') && ~isfield(par, 'compute_stresses_bottom')
    error('No stress files exist.')
end

% get items, all are surface integrals over the entire surface
strt = 20;
if isfield(par, 'compute_stresses_bottom')
    Bot_time   = stress_bot.Time(strt:end);
    Bot_tx_tot = stress_bot.Bottom_tx_tot(strt:end);
    Bot_tx_abs = stress_bot.Bottom_tx_abs(strt:end);
    Bot_ty_tot = stress_bot.Bottom_ty_tot(strt:end);
    Bot_ty_abs = stress_bot.Bottom_ty_abs(strt:end);
    Bot_ts     = stress_bot.Bottom_ts(strt:end);
end
if isfield(par, 'compute_stresses_top')
    Top_time   = stress_top.Time(strt:end);
    Top_tx_tot = stress_top.Top_tx_tot(strt:end);
    Top_tx_abs = stress_top.Top_tx_abs(strt:end);
    Top_ty_tot = stress_top.Top_ty_tot(strt:end);
    Top_ty_abs = stress_top.Top_ty_abs(strt:end);
    Top_ts     = stress_top.Top_ts(strt:end);
end

% create legends
if par.Ny == 1
    lab = {'$\int t_x dx$','$\int |t_x| dx$'};
else
    lab = {'$\int t_x dA$','$\int t_y dA$','$\int |t_x| dA$','$\int |t_y| dA$','$\int t_s dA$'};
end

% plot
figure(19), clf
if isfield(par, 'compute_stresses_bottom')
    subplot(2,1,1)
    if par.Ny == 1
        plot(Bot_time,[Bot_tx_tot, Bot_tx_abs]);
    else
        plot(Bot_time,[Bot_tx_tot, Bot_ty_tot, Bot_tx_abs, Bot_ty_abs, Bot_ts]);
    end
    xlabel('Time (s)')
    ylabel('Stress (N)')
    title('Bottom stresses')
    legend(lab,'Interpreter','Latex','FontSize',14)
    legend('location','best')
    legend('boxoff')
end

if isfield(par, 'compute_stresses_top')
    subplot(2,1,2)
    if par.Ny == 1
        plot(Top_time,[Top_tx_tot, Top_tx_abs]);
    else
        plot(Top_time,[Top_tx_tot, Top_ty_tot, Top_tx_abs, Top_ty_abs, Top_ts]);
    end
    xlabel('Time (s)')
    ylabel('Stress (N)')
    title('Top stresses')
    legend(lab,'Interpreter','Latex','FontSize',14)
    legend('location','best')
    legend('boxoff')
end
