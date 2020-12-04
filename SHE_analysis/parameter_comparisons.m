clc; clearvars; close all; 

DATA = SPINS_parameter_collator(list_filename);
set(groot, 'DefaultAxesXDir', 'normal');
%% Look at generation vs amplitude
figure(1)
subaxis(1, 3, 1)
displacement = DATA.pyc_loc - DATA.pyc_adj_loc;
DATA=DATA(displacement>.06, :);
displacement = DATA.pyc_loc - DATA.pyc_adj_loc;

[m, ~, R2] = plot_slope(displacement./DATA.Lz, DATA.meanAmp./DATA.Lz, 1);
xlabel('$(H_t - H_0)/H$', 'interpreter', 'latex');
ylabel('$A/H$', 'interpreter', 'latex');
text(.5, .05, {['$A_{sw} = $'], ['$',num2str(m), ' (H_t - H_0)$']}, 'HorizontalAlignment', 'center', 'interpreter', 'latex');

%% Look at generation vs wavelength
subaxis(1, 3, 2)

[m, c, r2] = plot_slope(displacement./DATA.Lz, DATA.meanWavelength./DATA.Lz, 1);
xlabel('$(H_t - H_0)/H$', 'interpreter', 'latex');
ylabel('$L_{sw}/H$', 'interpreter', 'latex');
text(.5, 1.7, {['$L_{sw} = ', num2str(m), '(H_t - H_0)$']}, 'HorizontalAlignment', 'center', 'interpreter', 'latex');

%% Look at generation vs Wave Slope
 figure(1)
 subaxis(1, 3, 3)
 [m, c, r2] = plot_slope(displacement./DATA.Lz, DATA.meanAmp./DATA.meanWavelength, 1);
 %[~, ~, ~] = plot_slope(displacement./DATA.Lz, DATA.meanAmp./DATA.meanWavelength, 2);
text(.5, .02, {['$A/L = ', num2str(m), '(H_t - H_0)$']}, 'HorizontalAlignment', 'center', 'interpreter', 'latex');

 xlabel('$(H_t - H_0)/H$', 'interpreter', 'latex');
 ylabel('$A/L_{sw}$', 'interpreter', 'latex');

 %% Format figure 1
 set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [2 2 30 17]); set(gcf, 'PaperPosition', get(gcf, 'Position'));
print('C:\Users\samha\OneDrive - Newcastle University\Project\Shoal_Core\figures\step_size_results.png', '-dpng');

 
%% Plot Sw and S for our waves
WaveSteepness = DATA.meanAmp./DATA.meanWavelength;
figure
plot(WaveSteepness, DATA.hill_slope, 'bx')
set(gca, 'xlim', [0 0.2]);
xlabel('$S_w$', 'interpreter', 'latex')
set(gca, 'ylim', [0 0.35]);
ylabel('$S$', 'interpreter', 'latex');

load 'Aghsaee2010.mat'
indSurges = find(strcmp(Aghsaee2010.ShoalingType, 'S'));
indColSurges = find(strcmp(Aghsaee2010.ShoalingType, 'CS'));
indFission = find(strcmp(Aghsaee2010.ShoalingType, 'F'));
indCollapse = find(strcmp(Aghsaee2010.ShoalingType, 'C'));
indSurges = find(strcmp(Aghsaee2010.ShoalingType, 'S'));
indColPlung = find(strcmp(Aghsaee2010.ShoalingType, 'CP'));
indPlunging = find(strcmp(Aghsaee2010.ShoalingType, 'P'));

hold on
plot(Aghsaee2010.Amp./Aghsaee2010.Wavelength, Aghsaee2010.Slope, 'x', 'color', [.2 .2 .2]);


%% Plot Sw and S for our waves
WaveIr = DATA.hill_slope./sqrt(DATA.meanAmp./DATA.meanWavelength);
figure
plot(WaveIr, DATA.hill_slope, 'bx')
set(gca, 'xlim', [0 0.2]);
xlabel('$S_w$', 'interpreter', 'latex')
set(gca, 'ylim', [0 0.35]);
ylabel('$S$', 'interpreter', 'latex');

load 'Aghsaee2010.mat'
hold on
WaveIr_AG2010 = Aghsaee2010.Slope./sqrt(Aghsaee2010.Amp./Aghsaee2010.Wavelength);
plot(WaveIr_AG2010, Aghsaee2010.Slope, 'x', 'color', [.2 .2 .2]);
print('C:\Users\samha\OneDrive - Newcastle University\Project\Shoal_Core\figures\parameter_space.png', '-dpng');

%%


function [m, c,  r2] = plot_slope(x, y, nums)

plot(x, y, 'kx')
range_y = range(y);
ylim([min(y)-.05*range_y max(y)+.05*range_y]);

%% Regression part
[a, b] = size(x);
if a<b
    x = x';
    y = y';
end
if nums==1
    X = x;
elseif nums ==2
    X = [ones(length(x), 1) x];
end
b1 = X\y;
yCalc1 = X*b1;
hold on
plot(x, yCalc1, 'k-')

r2 = 1 - sum((y-yCalc1).^2)/sum((y - mean(y)).^2);
if nums == 1
    m = b1(1);
    c = 0;
else
    m = b1(2);
    c = b1(1);
end
title(['R^2 = ', num2str(r2)]);
end

