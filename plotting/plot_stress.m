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
if par.perturb == 0
    strt = 1;
else
    strt = 20;
end
if isfield(par, 'compute_stresses_bottom')
    % integral values
    Bot_time   = stress_bot.Time(strt:end);
    Bot_tx_tot = stress_bot.Bottom_tx_tot(strt:end);
    Bot_tx_abs = stress_bot.Bottom_tx_abs(strt:end);
    Bot_ty_tot = stress_bot.Bottom_ty_tot(strt:end);
    Bot_ty_abs = stress_bot.Bottom_ty_abs(strt:end);
    Bot_ts     = stress_bot.Bottom_ts(strt:end);
    try
        % extrema values
        Bot_tx_max = stress_bot.Bottom_tx_max(strt:end);
        Bot_tx_min = stress_bot.Bottom_tx_min(strt:end);
        Bot_ty_max = stress_bot.Bottom_ty_max(strt:end);
        Bot_ty_min = stress_bot.Bottom_ty_min(strt:end);
        Bot_ts_max = stress_bot.Bottom_ts_max(strt:end);
        Bot_max = true;
    catch
        Bot_max = false;
    end
end
if isfield(par, 'compute_stresses_top')
    % integral values
    Top_time   = stress_top.Time(strt:end);
    Top_tx_tot = stress_top.Top_tx_tot(strt:end);
    Top_tx_abs = stress_top.Top_tx_abs(strt:end);
    Top_ty_tot = stress_top.Top_ty_tot(strt:end);
    Top_ty_abs = stress_top.Top_ty_abs(strt:end);
    Top_ts     = stress_top.Top_ts(strt:end);
    try
        % extrema values
        Top_tx_max = stress_top.Top_tx_max(strt:end);
        Top_tx_min = stress_top.Top_tx_min(strt:end);
        Top_ty_max = stress_top.Top_ty_max(strt:end);
        Top_ty_min = stress_top.Top_ty_min(strt:end);
        Top_ts_max = stress_top.Top_ts_max(strt:end);
        Top_max = true;
    catch
        Top_max = false;
    end
end

% integral values
% create legends
if par.Ny == 1
    lab = {'$\int t_x dx$','$\int |t_x| dx$'};
    ylab = 'Stress (N/m)';
else
    lab = {'$\int t_x dA$','$\int t_y dA$','$\int |t_x| dA$','$\int |t_y| dA$','$\int t_s dA$'};
    ylab = 'Stress (N/m^2)';
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
    ylabel('Force (N)')
    title('Bottom Forces')
    %legend(lab,'Interpreter','Latex','FontSize',14)
    %legend('location','best')
    %legend('boxoff')
    grid on
end

if isfield(par, 'compute_stresses_top')
    subplot(2,1,2)
    if par.Ny == 1
        plot(Top_time,[Top_tx_tot, Top_tx_abs]);
    else
        plot(Top_time,[Top_tx_tot, Top_ty_tot, Top_tx_abs, Top_ty_abs, Top_ts]);
    end
    xlabel('Time (s)')
    ylabel('Force (N)')
    title('Top Forces')
    legend(lab,'Interpreter','Latex','FontSize',12)
    legend('location','best')
    legend('boxoff')
    grid on
end


% extrema values
if Bot_max || Top_max
    % create legends
    if par.Ny == 1
        lab = {'max $t_x$','min $t_x$'};
    else
        lab = {'max $t_x$','max $t_y$','min $t_x$','min $t_y$','max $t_s$'};
    end

    % plot
    figure(18), clf
    if isfield(par, 'compute_stresses_bottom')
        subplot(2,1,1)
        if par.Ny == 1
            plot(Bot_time,[Bot_tx_max, Bot_tx_min]);
        else
            plot(Bot_time,[Bot_tx_max, Bot_ty_max, Bot_tx_min, Bot_ty_min, Bot_ts_max]);
        end
        xlabel('Time (s)')
        ylabel(ylab)
        title('Bottom extrema stresses')
        %legend(lab,'Interpreter','Latex','FontSize',14)
        %legend('location','best')
        %legend('boxoff')
        grid on
    end

    if isfield(par, 'compute_stresses_top')
        subplot(2,1,2)
        if par.Ny == 1
            plot(Top_time,[Top_tx_max, Top_tx_min]);
        else
            plot(Top_time,[Top_tx_max, Top_ty_max, Top_tx_min, Top_ty_min, Top_ts_max]);
        end
        xlabel('Time (s)')
        ylabel(ylab)
        title('Top extrema stresses')
        legend(lab,'Interpreter','Latex','FontSize',12)
        legend('location','best')
        legend('boxoff')
        grid on
    end
else
    figure(18)
    close
end
