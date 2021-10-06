clearvars; close all; clc;
runnames = {'091020_45', '101020_46', '111020_47', '121020_48',...
                '101120_49', '120220_50', '130220_51', '140220_52', '150220_53', ...
                '160220_54', '170220_55', '180220_56'}
for i = 1:length(runnames)
    cd(['../', runnames{i}]);
    
    load('wave_characteristics', 'WaveStats')
    params = spins_params;
    d = params.Lz;
    c = WaveStats.meanWaveSpeed;
    nu = params.visco;
    
    Re(i) = c*d/nu;
end
Re = Re';
