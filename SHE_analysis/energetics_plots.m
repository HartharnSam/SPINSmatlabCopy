%% Script to plot diagnostics we're actually interested in for energy changes
clearvars; close all; clc

all_diagnos = plot_diagnos(false, false, false);
EnergyBudget_inds = find(all_diagnos.EnergyBudget.Time > 10);
Mixing_inds = find(all_diagnos.Mixing.Time > 10);
diagnos_inds = find(all_diagnos.diagnos.Time > 10);


%subplot(2, 1, 1)
plot(all_diagnos.EnergyBudget.Time(EnergyBudget_inds), all_diagnos.EnergyBudget.KE2Int_tot(EnergyBudget_inds), 'b-');
hold on
plot(all_diagnos.EnergyBudget.Time(EnergyBudget_inds), all_diagnos.EnergyBudget.APE2BPE_tot(EnergyBudget_inds), 'r-');
 xlabel('time')
 ylabel('Energy Converted (J)');
 set(gca, 'XDir', 'normal');
 legend('Dissipated', 'Mixed', 'Location', 'northwest');
 %%
%  subplot(2, 1, 2)
%  plot(all_diagnos.Mixing.Time(Mixing_inds), all_diagnos.Mixing.mix_eff_cum(Mixing_inds), 'r-');
%  hold on
%  plot(all_diagnos.diagnos.Time(diagnos_inds), cumsum(all_diagnos.diagnos.Diss_tot(diagnos_inds)), 'b-');
%  
%  xlabel('time')
%  ylabel('');
%  
%  set(gca, 'XDir', 'normal')
%  legend('Mixing Efficiency', 'Dissipation')
