% plots things from the stresses diagnostic file

% read the file
try
    stressfile = 'stresses.txt';
    stress = readtable(stressfile);
catch
    error('stresses.txt not found.')
end

% find rho_0 to adjust if necessary
par = spins_params();
if isfield(par,'rho_0')
    if par.rho_0 == 1
        mult = 1000;    % use physical reference density
    else
        mult = 1;       % otherwise don't change anything
    end
else
    mult = 1000;
    disp('rho_0 not found in spins.conf. Setting rho_0 = 1000 kg/m^3.')
end

% get items, all are surface integrals over the entire surface
strt = 20;
time = stress.Time(strt:end);
Bot_tx_tot = stress.Bottom_tx_tot(strt:end)*mult;
Bot_tx_abs = stress.Bottom_tx_abs(strt:end)*mult;
Bot_ty_tot = stress.Bottom_ty_tot(strt:end)*mult;
Bot_ty_abs = stress.Bottom_ty_abs(strt:end)*mult;
Bot_ts     = stress.Bottom_ts(strt:end)*mult;
Top_tx_tot = stress.Top_tx_tot(strt:end)*mult;
Top_tx_abs = stress.Top_tx_abs(strt:end)*mult;
Top_ty_tot = stress.Top_ty_tot(strt:end)*mult;
Top_ty_abs = stress.Top_ty_abs(strt:end)*mult;
Top_ts     = stress.Top_ts(strt:end)*mult;
% legends
lab1 = {'Bottom t_x','Bottom t_y','Top t_x','Top t_y','Bottom t_s','Top t_s'};
lab2 = {'Bottom |t_x|','Bottom |t_y|','Top |t_x|','Top |t_y|','Bottom t_s','Top t_s'};

% plot
n = 30;
figure(n), clf
plot(time,[Bot_tx_tot, Bot_ty_tot, Top_tx_tot, Top_ty_tot, Bot_ts, Top_ts]);
xlabel('Time (s)')
ylabel('Stress (N)')
legend(lab1)
legend('location','best')
legend('boxoff')

n = n+1;
figure(n), clf
plot(time,[Bot_tx_abs, Bot_ty_abs, Top_tx_abs, Top_ty_abs, Bot_ts, Top_ts]);
xlabel('Time (s)')
ylabel('Stress (N)')
legend(lab2)
legend('location','best')
legend('boxoff')
