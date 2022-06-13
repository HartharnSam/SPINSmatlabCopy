function all_stress = plot_stress(make_plots, save_plots)
%  PLOT_STRESS  Plot surface stresses and surface forces
%
%  David Deepwell, 2019

%% set default plotting option
if nargin == 0
    make_plots = true;
end

% read parameters
par = spins_params();

% which stresses were computed?
if isfield(par, 'compute_stresses_bottom')
    if strcmp(par.compute_stresses_bottom, 'true')
        compute_bottom = true;
    else
        compute_bottom = false;
    end
else
    compute_bottom = false;
end
if isfield(par, 'compute_stresses_top')
    if strcmp(par.compute_stresses_top, 'true')
        compute_top = true;
    else
        compute_top = false;
    end
else
    compute_top = false;
end


% read the file
if compute_bottom
    check_txt_file('stresses_bottom.txt', 'Time')
    stress_bot = readtable('stresses_bottom.txt');
    all_stress.bottom = stress_bot;
end
if compute_top
    check_txt_file('stresses_top.txt',    'Time')
    stress_top = readtable('stresses_top.txt');
    all_stress.top = stress_top;
end
% error check
if ~compute_bottom && ~compute_top
    error('No stress files exist.')
end

% get items, all are surface integrals over the entire surface
if par.perturb == 0
    strt = 1;
else
    strt = 20;
end
if compute_bottom
    % bottom surface forces
    Bot_time   = stress_bot.Time(strt:end);
    Bot_fx_tot = stress_bot.fx_tot(strt:end);
    Bot_fx_abs = stress_bot.fx_abs(strt:end);
    Bot_fy_tot = stress_bot.fy_tot(strt:end);
    Bot_fy_abs = stress_bot.fy_abs(strt:end);
    % Bottom stress extrema
    Bot_tx_max = stress_bot.tx_max(strt:end);
    Bot_tx_min = stress_bot.tx_min(strt:end);
    Bot_ty_max = stress_bot.ty_max(strt:end);
    Bot_ty_min = stress_bot.ty_min(strt:end);
    Bot_ts_max = stress_bot.ts_max(strt:end);
end
if compute_top
    % top surface forces
    Top_time   = stress_top.Time(strt:end);
    Top_fx_tot = stress_top.fx_tot(strt:end);
    Top_fx_abs = stress_top.fx_abs(strt:end);
    Top_fy_tot = stress_top.fy_tot(strt:end);
    Top_fy_abs = stress_top.fy_abs(strt:end);
    % extrema values
    Top_tx_max = stress_top.tx_max(strt:end);
    Top_tx_min = stress_top.tx_min(strt:end);
    Top_ty_max = stress_top.ty_max(strt:end);
    Top_ty_min = stress_top.ty_min(strt:end);
    Top_ts_max = stress_top.ts_max(strt:end);
end

if make_plots
    %% Surface forces
    % create legends
    if par.Ny == 1
        lab = {'$\int t_x dx$','$\int |t_x| dx$'};
        ylab = 'Stress (N/m)';
    else
        lab = {'$\int t_x dA$','$\int t_y dA$','$\int |t_x| dA$','$\int |t_y| dA$'};
        ylab = 'Stress (N/m^2)';
    end
    
    % plot
    if compute_top
        figure(19), clf
        
        subplot(2,1,1)
        if par.Ny == 1
            plot(Top_time,[Top_fx_tot, Top_fx_abs]);
        else
            plot(Top_time,[Top_fx_tot, Top_fy_tot, Top_fx_abs, Top_fy_abs]);
        end
        xlabel('t (s)')
        ylabel('Force (N)')
        title('Top Forces')
        %legend(lab,'Interpreter','Latex','FontSize',12)
        %legend('location','best')
        %legend('boxoff')
        grid on
        
        subplot(2,1,2)
        if par.Ny == 1
            plot(Top_time,[Top_tx_max, Top_tx_min]);
        else
            plot(Top_time,[Top_tx_max, Top_ty_max, Top_tx_min, Top_ty_min, Top_ts_max]);
        end
        xlabel('t (s)')
        ylabel(ylab)
        title('Top extrema stresses')
        %legend(lab,'Interpreter','Latex','FontSize',12)
        %legend('location','best')
        %legend('boxoff')
        grid on
        if save_plots
            print('TopStresses.png', '-dpng');
        end
    end
    
    %% extrema values
    % create legends
    if par.Ny == 1
        lab = {'max $t_x$','min $t_x$'};
    else
        lab = {'max $t_x$','max $t_y$','min $t_x$','min $t_y$','max $t_s$'};
    end
    
    % plot
    figure(18), clf
    if compute_bottom
        subplot(2,1,1)
        if par.Ny == 1
            plot(Bot_time,[Bot_fx_tot, Bot_fx_abs]);
        else
            plot(Bot_time,[Bot_fx_tot, Bot_fy_tot, Bot_fx_abs, Bot_fy_abs]);
        end
        xlabel('t (s)')
        ylabel('Force (N)')
        title('Bottom Forces')
        legend(lab,'Interpreter','Latex','FontSize',12)
        legend('location','best')
        legend('boxoff')
        grid on
        
        subplot(2,1,2)
        if par.Ny == 1
            plot(Bot_time,[Bot_tx_max, Bot_tx_min]);
        else
            plot(Bot_time,[Bot_tx_max, Bot_ty_max, Bot_tx_min, Bot_ty_min, Bot_ts_max]);
        end
        xlabel('t (s)')
        ylabel(ylab)
        title('Bottom extrema stresses')
        legend(lab,'Interpreter','Latex','FontSize',12)
        legend('location','best')
        legend('boxoff')
        grid on
        
        if save_plots
            print('BottomStresses.png', '-dpng');
        end
    end
end

function check_txt_file(filename, header_ptrn)
    % error check if multiple runs have been completed in this directory
    file_txt = fileread(filename);
    file_lines = regexp(file_txt, '\r\n|\r|\n', 'split');
    if sum(contains(file_lines, header_ptrn)) > 1
        error('%s file has too many headers\n%s',...
            filename,...
            'Remove headers and outputs from previous runs before continuing')
    end
end
