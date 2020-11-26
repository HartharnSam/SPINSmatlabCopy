%clearvars
close all;
%% Two Layer
runnames = {'200520_02', '250520_03', '270520_04', '280520_05', '040620_06', ...
    '050620_07', '070620_08', '080620_09', '080620_10', '090620_11', '110620_12', ...
    '120620_13', '130620_14', '140620_15', '150620_16', '160620_17', '170620_18'...
    '230620_19', '240620_20', '250620_21', '260620_22', '270620_23', '280620_24',...
    '010720_25', '020720_26', '030720_27', '040720_28', '050720_29', '250720_31',...
    '260720_32', '270720_33', '300720_36', '081020_44'};
runnames = {'250720_31', '270520_04'};
%% Three Layer
%runnames = {'02_090720', '03_100720', '04_110720', '05_120720',...
%    '07_310720', '24_071020', '25_221020', '26_091120', '27_111120', '28_121120'};
%runnames = {'27_111120', '26_091120', '24_071020', '08_010820'};
%% One Layer
%runnames = {'091020_45', '101020_46', '111020_47', '121020_48', '101120_49'}

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
