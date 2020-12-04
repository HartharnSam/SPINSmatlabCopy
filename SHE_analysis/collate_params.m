%clearvars
close all;
runnames = run_directory_names('2 Layer');

%%
for i = 1:length(runnames)
    cd(['../', runnames{i}]);
    params = spins_params;
    Lx(i) = params.Lx;
    Ly(i) = params.Ly;
    Lz(i) = params.Lz;
    pyc_thickness(i) = params.h_halfwidth;
    h1(i) = Lz(i) + (params.pyc_loc - pyc_thickness(i));
    
    all_diagnos = plot_diagnos(false, false, false);
    APE_Loss(i, 1) = all_diagnos.EnergeticChange.APE_Lost;
    APE_toMix(i, 1) = all_diagnos.EnergeticChange.APE_toMix;
    APE_toDiss(i, 1) = all_diagnos.EnergeticChange.APE_toDiss;
    
    all_diagnos = find_diss_peaks;
    MaxMaxDiss_time(i, 1) = all_diagnos.Energetics.MaxDiss_time; % Time when maximum dissipation (point) reached - should indicate break point
    MaxMaxDiss(i, 1) = all_diagnos.Energetics.MaxDiss_tot; % Maximum instantaneous dissipation measurement
    plot_diss_uc_Ri(MaxMaxDiss_time(i,1));
    print('BreakPoint_Diagnostic.png', '-dpng')
    
    %close all;
end
Lx = Lx';
Ly = Ly';
Lz = Lz';
pyc_thickness = pyc_thickness';
h1 = h1';
